import random
import heapq
import re
from typing import List, Dict, Tuple, Optional
from .utils import (
    make_hairpin, generate_random_nucleotides, hairpin_ok_at_42C,
    calculate_tm, calculate_gc_content, _dg, reverse_complement,
    reverse, complement, vienna_fold_single, heterodimer_ascii_onepic,
    vienna_cofold_dotbracket, VIENNA_AVAILABLE
)

def check_cross_reactivity(rt_primer: str, variant_seq: str) -> Dict:
    """
    Check if the RT primer cross-reacts with a variant sequence.
    Returns binding energy (dG) and structure if ViennaRNA is available.
    """
    if not VIENNA_AVAILABLE:
        return {"error": "ViennaRNA not available"}
    
    # Calculate interaction between variant (RNA) and RT primer (DNA)
    struct, mfe = vienna_cofold_dotbracket(variant_seq, rt_primer)
    return {
        "mfe": mfe,
        "struct": struct
    }

def build_rt_primer(rt5: str,
                    rt3: str,
                    *,
                    stem_len=(7, 8),
                    loop_len=(5, 5),
                    rand5_len=(8, 10),
                    rand3_len=(8, 10),
                    n_rep=300,
                    hairpin_ok_fn=None):
    """
    Try up to n_rep times to build a primer:
      5'- rt5 + random_5 + hairpin(stem_len,loop_len) + random_3 + rt3 -3'
    Lengths can be fixed ints, inclusive ranges, or discrete lists.
    hairpin_ok_fn: optional callable(rt_primer)->bool (e.g., hairpin_ok_at_42C)
    """
    for _ in range(n_rep):
        hairpin = make_hairpin(stem_len=stem_len, loop_len=loop_len)
        rand5 = generate_random_nucleotides(rand5_len)
        rand3 = generate_random_nucleotides(rand3_len)
        rt_primer = rt5 + rand5 + hairpin + rand3 + rt3
        if hairpin_ok_fn is None or hairpin_ok_fn(rt_primer):
            return rt_primer, rand5, rand3, hairpin
    return None

def tune_primer(
    primer_seq,
    target_tm_range,
    oligo_conc, Na_conc, Mg_conc, dNTPs_conc,
    max_length=22, min_length=16,
    context_5prime=None,       # optional: upstream sequence in 5'→3', with the base immediately upstream at the RIGHT end (…UUU|[HERE])
    return_nearest=False,      # if True, return the closest Tm if exact window cannot be met
    max_steps=100,             # safety cap
):
    """
    Adjust primer Tm into target_tm_range by ONLY modifying the 5′ end.
    - If Tm too LOW: prepend bases to 5′ end, first from context_5prime (if provided),
      then GC/A 'clamps' (G/C preferred, occasional A as gentle bump).
    - If Tm too HIGH: trim from the 5′ end.
    The 3′ end is never altered.

    Parameters
    ----------
    primer_seq : str
        Initial primer (5'→3').
    target_tm_range : (float, float)
        Desired Tm window in °C, inclusive (e.g., (58, 65)).
    oligo_conc, Na_conc, Mg_conc, dNTPs_conc : float
        Concentrations for Tm calculation.
    max_length, min_length : int
        Allowed final primer length bounds.
    context_5prime : str or None
        Upstream sequence (5'→3'); the base immediately upstream of `primer_seq`
        MUST be the RIGHTMOST char of this string. Example:
            template:  ... A  G  [P r i m e r ...]
                                 ^ rightmost context base goes here
    return_nearest : bool
        If True, returns the closest (seq, Tm) even if outside the target range.
        If False, raises ValueError if window can’t be met.
    max_steps : int
        Safety cap on total edit steps.

    Returns
    -------
    (str, float)
        Tuned primer sequence (5'→3') and its Tm (°C).
    """
    lo, hi = target_tm_range
    if min_length > max_length:
        raise ValueError("min_length cannot be greater than max_length")

    def tm(seq): return calculate_tm(seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc)

    # internal helper: choose a clamp base with GC preference
    def choose_clamp():
        # weights: G:5, C:5, A:1
        pool = ['G']*5 + ['C']*5 + ['A']
        return random.choice(pool)

    # current state
    seq = primer_seq
    best = (abs((tm(seq)) - (lo+hi)/2.0), seq, tm(seq))  # (distance, seq, tm)

    # for context pulls: consume from RIGHT end (closest to primer)
    ctx = list(context_5prime) if context_5prime else None

    steps = 0
    while steps < max_steps:
        steps += 1
        cur_tm = tm(seq)
        # record best-so-far
        d = abs(cur_tm - (lo+hi)/2.0)
        if d < best[0]:
            best = (d, seq, cur_tm)

        # already good and within length bounds
        if lo <= cur_tm <= hi and (min_length <= len(seq) <= max_length):
            return seq, cur_tm

        # too low: try to prepend from context, else a clamp
        if cur_tm < lo:
            if len(seq) >= max_length:
                break  # cannot extend
            if ctx and len(ctx) > 0:
                base = ctx.pop()  # take rightmost base (immediately upstream)
                seq = base + seq
            else:
                # clamp (avoid creating long homopolymers at 5′ if possible)
                base = choose_clamp()
                if len(seq) >= 4 and seq[:4] == base*4:
                    # nudge away from long runs
                    base = 'G' if base != 'G' else 'C'
                seq = base + seq
            continue

        # too high: trim from 5′ if possible
        if cur_tm > hi:
            if len(seq) <= min_length:
                break  # cannot trim
            seq = seq[1:]
            continue

    # If we reach here, we couldn’t hit the exact window
    if return_nearest and (min_length <= len(best[1]) <= max_length):
        return best[1], best[2]
    raise ValueError("Could not adjust primer to meet Tm and length criteria")

# ---------------- Scoring & search ----------------
def score_design(rt_tm, f_tm, r_tm, gc, dG3, dG5, length):
    # Heuristic scoring consistent with your targets (Tm ~60, GC 40–60, small length bonus)
    score = 0.0
    score += max(0, 20 - abs(gc - 50))
    score += max(0, 20 - abs(f_tm - 60))
    score += max(0, 20 - abs(r_tm - 60))
    score += max(0, 15 - abs(f_tm - r_tm))
    # enforce windows already, add small reward for being closer to −10
    score += max(0, 5 - abs((dG3 - (-10))))
    score += max(0, 5 - abs((dG5 - (-10))))
    # avoid long homopolymers
    score -= 15 if re.search(r'(A{5,}|T{5,}|G{5,}|C{5,})', 'X'*0) else 0  # placeholder; RT already filtered
    # shorter is slightly better beyond ~80 nt
    score += max(0, 10 - 0.1*(length - 80))
    return score

def enumerate_strict_candidates(miDNA, Na, Mg, dNTPs, C_rt, C_qpcr,
                                dg_min=-14.0, dg_max=-8.0,
                                top_n=5, seed=13,
                                # iterative knobs:
                                tries_per_pair=120,
                                beam_per_pair=3,
                                stem_len=(7, 8),
                                loop_len=(5, 5),
                                rand5_len=(8, 10),
                                rand3_len=(8, 10),
                                tm_window=(58, 65),
                                max_primer_len=25,
                                min_primer_len=16,
                                # NEW:
                                l3_min=4, l3_max=8,
                                l5_min=4, l5_max=8,
                                dg_sum_max=-20.0,
                                dg_temp_C=37.0,
                                variant_sequences: List[str] = None,
                                p5_range: Tuple[Optional[int], Optional[int]] = (None, None),
                                p3_range: Tuple[Optional[int], Optional[int]] = (None, None)): # Added p3_range
    """
    Iterative search for best designs; returns top_n by score.
    """
    random.seed(seed)
    n = len(miDNA)

    # Precompute cDNA (DNA alphabet) -- complement (not reverse) per your snippet
    cDNA_seq = complement(miDNA)
    
    # Convert temp to Kelvin for _dg
    temp_k = 273.15 + dg_temp_C

    p3_min, p3_max_arg = p3_range

    # ---- 3′ hemiprobe: length l3_min–l3_max, ΔG in [dg_min, dg_max]
    # Now supports floating position if p3_range is specified or default behavior
    valid3 = []
    for L3 in range(l3_min, l3_max + 1):
        if L3 > n: break
        # Iterate over all possible start positions for 3' probe
        for s3 in range(0, n - L3 + 1):
            pos_1based = s3 + 1
            
            # Check p3 constraints
            if p3_min is not None and pos_1based < p3_min: continue
            if p3_max_arg is not None and pos_1based > p3_max_arg: continue
            
            # If no p3 constraints are provided, default to anchoring at the 3' end?
            # The original behavior was anchored at 3' end (s3 = n - L3).
            # If the user leaves p3 inputs empty (None), should we allow floating or enforce anchor?
            # To preserve original behavior when inputs are empty: enforce anchor.
            if p3_min is None and p3_max_arg is None:
                if s3 != n - L3: continue

            comp3 = miDNA[s3:s3+L3]
            dG3   = _dg(comp3, temp_k=temp_k) * 2
            if dg_min <= dG3 <= dg_max:
                rt3 = reverse_complement(comp3)
                valid3.append((s3, L3, comp3, rt3, dG3))

    p5_min, p5_max_arg = p5_range

    # ---- 5′ hemiprobe: 5′ end OR middle, length l5_min–l5_max, ΔG in [dg_min, dg_max]
    valid5_all = []
    for L5 in range(l5_min, l5_max + 1):
        if L5 >= n: break
        for s in range(0, n - L5 + 1):
            # Check position constraints
            pos_1based = s + 1
            if p5_min is not None and pos_1based < p5_min: continue
            if p5_max_arg is not None and pos_1based > p5_max_arg: continue

            seg = miDNA[s:s+L5]
            dG5 = _dg(seg, temp_k=temp_k) * 2
            if dg_min <= dG5 <= dg_max:
                rt5 = reverse_complement(seg)
                valid5_all.append((s, L5, seg, rt5, dG5))

    global_heap = []  # min-heap of (-score, unique_id, cand_dict)
    seen_rt_primers = set()
    unique_id_counter = 0

    # Iterate pairs
    for (s3, L3, comp3, rt3, dG3) in valid3:
        # per-pair beam
        pair_heap = []
        pair_seq_set = set()

        for (s5, L5, comp5, rt5, dG5) in valid5_all:
            # no overlap; 3′ probe must be downstream of 5' probe
            # s5 is start of 5' probe, s3 is start of 3' probe
            # We need 5' probe to end before 3' probe starts: s5 + L5 <= s3
            if s5 + L5 > s3: 
                continue
            if not (dG3 <= dG5 and (dG3 + dG5) <= dg_sum_max):
                continue

            # Iterate randomizations for this pair
            for _ in range(tries_per_pair):
                built = build_rt_primer(
                    rt5=rt5, rt3=rt3,
                    stem_len=stem_len, loop_len=loop_len,
                    rand5_len=rand5_len, rand3_len=rand3_len,
                    n_rep=1,  # one build attempt per inner iteration
                    hairpin_ok_fn=hairpin_ok_at_42C
                )
                if not built:
                    continue
                rt_primer, rand5, rand3, hairpin = built

                # RT reaction product template (what reverse primer sees): RT primer + reverse(cDNA)
                rt_product = rt5 + rand5 + hairpin + rand3 + reverse(cDNA_seq)

                # Forward: take first 18–20 nt, then tune into window
                seed_fwd = rt_primer[:18]
                try:
                    fwd, f_tm = tune_primer(seed_fwd, tm_window, C_qpcr, Na, Mg, dNTPs,
                                            max_length=max_primer_len, min_length=min_primer_len)
                except ValueError:
                    continue

                # Reverse: take last 18 nt of RT product (template), RC, then tune
                seed_rev_template = rt_product[-18:]
                seed_rev = reverse_complement(seed_rev_template)
                try:
                    rev, r_tm = tune_primer(seed_rev, tm_window, C_qpcr, Na, Mg, dNTPs,
                                            max_length=max_primer_len, min_length=min_primer_len)
                except ValueError:
                    continue

                # RT primer properties
                rt_tm = calculate_tm(rt_primer, C_rt, Na, Mg, dNTPs)
                rt_gc = calculate_gc_content(rt_primer)

                # (Optional) structure of RT at 42°C (should show only the hairpin)
                rt_struct, rt_mfe = vienna_fold_single(rt_primer, 42.0)

                # Colored HTML parts
                colored = (
                    f'<span style="color:#6A1B9A;">{rt5}</span>'
                    f'<span style="color:#F9A825;">{rand5}</span>'
                    f'<span style="color:#1E88E5;">{hairpin}</span>'
                    f'<span style="color:#1B5E20;">{rand3}</span>'
                    f'<span style="color:#D81B60;">{rt3}</span>'
                )

                # One-picture ASCII
                ascii_fig, basepairs = heterodimer_ascii_onepic(
                    miRNA_rna=miDNA.replace('T','U'),
                    s5=s5, L5=L5,
                    rt5=rt5,  # Pass DNA directly
                    L3=L3,
                    rt3=rt3   # Pass DNA directly
                )

                # Score candidate
                score = score_design(
                    rt_tm=rt_tm, f_tm=f_tm, r_tm=r_tm, gc=rt_gc,
                    dG3=dG3, dG5=dG5, length=len(rt_primer)
                )

                # Deduplicate within pair
                if rt_primer in pair_seq_set:
                    continue
                pair_seq_set.add(rt_primer)

                cand = {
                    "score": score,
                    "rt_primer": rt_primer,
                    "rt_product": rt_product,
                    "colored_html": colored,
                    "rt_primer_5": rt5, "random_nt_5": rand5, "hairpin": hairpin,
                    "random_nt_3": rand3, "rt_primer_3": rt3,
                    "RT_primer_Tm": f"{rt_tm:.1f}°C", "RT_primer_GC": f"{rt_gc:.1f}%",
                    "Forward_primer": fwd, "Forward_primer_Tm": f"{f_tm:.1f}°C",
                    "Forward_primer_GC": f"{calculate_gc_content(fwd):.1f}%",
                    "Reverse_primer": rev, "Reverse_primer_Tm": f"{r_tm:.1f}°C",
                    "Reverse_primer_GC": f"{calculate_gc_content(rev):.1f}%",
                    "dg3": f"{dG3:.2f} kcal/mol", "dg5": f"{dG5:.2f} kcal/mol",
                    "dg_sum": f"{(dG3 + dG5):.2f} kcal/mol",
                    "rt_struct_42C": rt_struct if rt_struct else "N/A",
                    "rt_mfe_42C": f"{rt_mfe:.2f} kcal/mol" if rt_mfe is not None else "N/A",
                    "heterodimer_ascii": ascii_fig, "heterodimer_pairs": basepairs,
                    "binding_5_start": s5 + 1,
                    "binding_5_end": s5 + L5,
                    "binding_3_start": n - L3 + 1,
                    "binding_3_end": n,
                    "binding_5": f"{s5 + 1}-{s5 + L5}",
                    "binding_3": f"{n - L3 + 1}-{n}",
                }

                if variant_sequences:
                    cand["cross_reactivity"] = []
                    for var_seq in variant_sequences:
                        if not var_seq.strip(): continue
                        xr = check_cross_reactivity(rt_primer, var_seq.strip())
                        xr['variant'] = var_seq.strip()
                        if 'mfe' in xr and isinstance(xr['mfe'], (int, float)):
                             xr['mfe'] = f"{xr['mfe']:.2f}"
                        cand["cross_reactivity"].append(xr)

                # push to per-pair beam (min-heap of size beam_per_pair)
                # push to per-pair beam (min-heap of size beam_per_pair)
                unique_id_counter += 1
                heapq.heappush(pair_heap, (score, unique_id_counter, cand))
                if len(pair_heap) > beam_per_pair:
                    heapq.heappop(pair_heap)

        # Merge pair beam into global heap, dedup on full RT primer
        while pair_heap:
            _, _, cand = heapq.heappop(pair_heap)
            rt_seq = cand["rt_primer"]
            if rt_seq in seen_rt_primers:
                continue
            seen_rt_primers.add(rt_seq)
            unique_id_counter += 1
            heapq.heappush(global_heap, (cand["score"], unique_id_counter, cand))
            if len(global_heap) > (top_n * 4):  # keep a small global pool
                heapq.heappop(global_heap)

    # Final top_n by score
    final = sorted((c for _, _, c in global_heap), key=lambda x: x["score"], reverse=True)
    return final[:top_n]
