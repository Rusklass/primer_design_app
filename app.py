# app.py

from flask import Flask, render_template, request, jsonify
import numpy as np
import random, os, re, heapq
from typing import Optional, List, Dict, Tuple

app = Flask(__name__)

# Optional: ViennaRNA (skip gracefully if not present)
try:
    import RNA  # pip install ViennaRNA (may need system package)
    VIENNA_AVAILABLE = True
except Exception:
    VIENNA_AVAILABLE = False

# ---------------- Thermo constants ----------------
R = 1.987    # cal/(K*mol)
# T = 298.15   # K (25°C) for ΔG window checks
T = 310.15   # K (37°C) for ΔG window checks

# DNA↔RNA maps
DNA2RNA = str.maketrans({'A':'A','T':'U','G':'G','C':'C'})
RNA2DNA = str.maketrans({'A':'A','U':'T','G':'G','C':'C'})

# Nearest-neighbor parameters (DNA/DNA)
nn_params = {
    'AA': {'dH': -7.9, 'dS': -22.2}, 'TT': {'dH': -7.9, 'dS': -22.2},
    'AT': {'dH': -7.2, 'dS': -20.4}, 'TA': {'dH': -7.2, 'dS': -21.3},
    'CA': {'dH': -8.5, 'dS': -22.7}, 'TG': {'dH': -8.5, 'dS': -22.7},
    'GT': {'dH': -8.4, 'dS': -22.4}, 'AC': {'dH': -8.4, 'dS': -22.4},
    'CT': {'dH': -7.8, 'dS': -21.0}, 'AG': {'dH': -7.8, 'dS': -21.0},
    'GA': {'dH': -8.2, 'dS': -22.2}, 'TC': {'dH': -8.2, 'dS': -22.2},
    'CG': {'dH': -10.6, 'dS': -27.2}, 'GC': {'dH': -9.8, 'dS': -24.4},
    'GG': {'dH': -8.0, 'dS': -19.9},  'CC': {'dH': -8.0, 'dS': -19.9},
}

# ---------------- Utilities ----------------
def reverse_complement(seq):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(comp.get(b, b) for b in seq[::-1])

def reverse(seq):
    return seq[::-1]

def complement(seq):
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement_map.get(base, base) for base in seq])

def calculate_gc_content(seq):
    return 0.0 if not seq else 100.0*(seq.count('G') + seq.count('C'))/len(seq)

def _nn_thermo(seq):
    """Return (ΔH kcal/mol, ΔS cal/K/mol) for DNA/DNA with simple init & terminal corrections."""
    if not seq: return 0.0, 0.0
    dH, dS = 0.2, -5.7  # initiation
    if seq[0] in ('A','T'): dH, dS = dH+2.2, dS+6.9
    if seq[-1] in ('A','T'): dH, dS = dH+2.2, dS+6.9
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        v = nn_params.get(pair) or nn_params.get(pair[::-1])
        if v: dH, dS = dH + v['dH'], dS + v['dS']
    return dH, dS # kcal/mol, cal/K/mol

def _dg(seq):
    dH, dS = _nn_thermo(seq)
    return dH - (T * dS / 1000.0)

def calculate_tm(seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc):
    dH, dS = _nn_thermo(seq)
    dH *= 1000.0
    oligo_conc_adjusted = max(oligo_conc/4.0, 1e-12)  # avoid log(0)
    tm_K = dH / (dS + R*np.log(oligo_conc_adjusted))
    tm_C = tm_K - 273.15
    # simplified salt correction
    salt_effect = 16.6 * np.log10(max(Na_conc + 4.0*np.sqrt(max(Mg_conc,0.0)), 1e-6))
    return tm_C + salt_effect

def generate_random_nucleotides(len_spec):
    """
    len_spec: int, (min,max), or list of allowed lengths.
    """
    L = _choose_len(len_spec, "random_nt_len", min_allowed=1)
    return ''.join(random.choice('AGCT') for _ in range(L))

def _choose_len(spec, name: str, *, min_allowed: int = 1) -> int:
    """
    Accepts:
      - int -> use as-is
      - (min,max) tuple/list -> pick randint(min,max) inclusive
      - list/tuple of discrete values (len != 2) -> pick random choice
    """
    if isinstance(spec, int):
        L = spec
    elif isinstance(spec, (tuple, list)):
        if len(spec) == 2 and all(isinstance(x, (int, float)) for x in spec):
            a, b = int(min(spec)), int(max(spec))
            L = random.randint(a, b)
        else:
            # discrete set of candidates, e.g. [8,9,10]
            if not spec:
                raise ValueError(f"{name}: empty list of candidates.")
            L = int(random.choice(spec))
    else:
        raise TypeError(f"{name}: expected int, (min,max), or list of candidates.")

    if L < min_allowed:
        raise ValueError(f"{name}: length {L} < minimum allowed {min_allowed}.")
    return L

def _has_forbidden_complementarity(loop: str, stem: str, min_k: int = 5) -> bool:
    """
    Return True if the loop contains reverse-complement matches of length >= min_k
    to either stem arm, or if the loop is self-complementary over >= min_k.
    """
    left_arm  = stem
    right_arm = reverse_complement(stem)  # the other arm sequence

    # Check loop vs BOTH stem arms
    for arm in (left_arm, right_arm):
        for k in range(min_k, min(len(loop), len(arm)) + 1):
            # scan all k-mers of the arm; if rc(kmer) appears in loop -> forbidden
            for i in range(len(arm) - k + 1):
                kmer = arm[i:i+k]
                if reverse_complement(kmer) in loop:
                    return True

    # Check loop self-complementarity (avoid internal mini-hairpins)
    for k in range(min_k, len(loop) + 1):
        for i in range(len(loop) - k + 1):
            kmer = loop[i:i+k]
            rc_kmer = reverse_complement(kmer)
            # ensure we’re not trivial same-position match; any occurrence is suspicious
            if rc_kmer in loop:
                # optionally, require a *separate* location:
                j = loop.find(rc_kmer)
                if j != -1 and not (j == i and rc_kmer == kmer[::-1]):  # allow exact palindromic overlap only if you want
                    return True
    return False

def make_hairpin(stem_len=(8, 8), loop_len=(5, 5), *, min_k: int = 4, max_tries: int = 1000) -> str:
    """
    Generate a hairpin:  [stem] + [loop] + rc(stem)
    while ensuring the loop does NOT contain reverse-complement segments
    (length >= min_k) to either stem arm, and has no strong self-complementarity.
    stem_len, loop_len can be:
      - int (fixed)
      - (min,max) inclusive tuple/list
      - list of discrete allowed values, e.g. [7,8,9]
    """
    def balanced_random(n: int) -> str:
        s = ''
        while len(s) < n:
            b = random.choice('AAGCTGCT')  # light GC bias
            if len(s) >= 3 and s[-3:] == b * 3:
                continue  # avoid >=3 homopolymers
            s += b
        return s

    for _ in range(max_tries):
        Ls = _choose_len(stem_len, "stem_len", min_allowed=3)
        Ll = _choose_len(loop_len, "loop_len", min_allowed=3)

        stem = balanced_random(Ls)

        # ensure stem GC is reasonable (e.g., 35–70%)
        gc = (stem.count('G') + stem.count('C')) / Ls
        if not (0.35 <= gc <= 0.70):
            continue

        # Try several loop candidates for this stem
        for _ in range(200):
            loop = balanced_random(Ll)

            # Quick filters: avoid loop homopolymers of 4, and extreme GC
            if re.search(r'(A{4,}|T{4,}|G{4,}|C{4,})', loop):
                continue
            gc_loop = (loop.count('G') + loop.count('C')) / Ll
            if not (0.2 <= gc_loop <= 0.8):
                continue

            if not _has_forbidden_complementarity(loop, stem, min_k=min_k):
                return stem + loop + reverse_complement(stem)

    raise RuntimeError("Failed to generate a valid hairpin without loop–stem complementarity.")

def hairpin_ok_at_42C(rt_primer: str) -> bool:
    if not VIENNA_AVAILABLE:
        return True  # skip check if library absent
    try:
        md = RNA.md(); md.temperature = 42.0
        fc = RNA.fold_compound(rt_primer.replace('T','U'), md)
        struct, _ = fc.mfe()
        # Expect a single hairpin: roughly "((((....))))" around the loop region.
        # You can add stricter checks (e.g., exact loop length) if desired.
        return '(' in struct and ')' in struct  # minimal sanity; customize as needed
    except Exception:
        return True

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

def get_bool(form, name, default=False):
    v = form.get(name)
    if v is None:
        return default
    return str(v).lower() in ('1','true','yes','on')

def colorize(rt3, r3rand, hairpin, r5rand, rt5):
    palette = {
        'rt3': '#D81B60',  # magenta
        'r3':  '#1B5E20',  # dark green
        'hp':  '#1E88E5',  # blue
        'r5':  '#F9A825',  # amber
        'rt5': '#6A1B9A',  # purple
    }
    return (
        f'<span style="color:{palette["rt3"]};">{rt3}</span>'
        f'<span style="color:{palette["r3"]};">{r3rand}</span>'
        f'<span style="color:{palette["hp"]};">{hairpin}</span>'
        f'<span style="color:{palette["r5"]};">{r5rand}</span>'
        f'<span style="color:{palette["rt5"]};">{rt5}</span>'
    )

def vienna_cofold_dotbracket(mirna_rna, rtprimer_dna):
    """Return (dotbracket, ΔG kcal/mol) or (None,None). Treat DNA primer as RNA by T->U."""
    if not VIENNA_AVAILABLE: return None, None
    s1 = mirna_rna.replace('T','U').upper()
    s2 = rtprimer_dna.replace('T','U').upper()
    try:
        struct, mfe = RNA.cofold(s1 + '&' + s2)
        return struct, float(mfe)
    except Exception:
        return None, None
    
def vienna_fold_single(seq_dna, tempC=42.0):
    if not VIENNA_AVAILABLE: return None, None
    try:
        md = RNA.md()
        md.temperature = float(tempC)
        fc = RNA.fold_compound(seq_dna.replace('T','U'), md)
        struct, mfe = fc.mfe()
        return struct, float(mfe)
    except Exception:
        return None, None

def _pairs_from_dotbracket(dot: str) -> Dict[int,int]:
    st, pairs = [], {}
    for i, ch in enumerate(dot):
        if ch == '(':
            st.append(i)
        elif ch == ')':
            if st:
                j = st.pop(); pairs[i]=j; pairs[j]=i
    return pairs


def fold_rt_and_check_hairpin_only(rt_primer_dna: str, hp_start: int, stem_len: int, loop_len: int,
                                   tempC: float = 42.0) -> Tuple[bool, Optional[str], Optional[float], Optional[str]]:
    if not VIENNA_AVAILABLE:
        return True, None, None, None
    try:
        hp_end = hp_start + stem_len + loop_len + stem_len - 1
        md = RNA.md(); md.temperature = float(tempC)
        fc = RNA.fold_compound(rt_primer_dna.translate(DNA2RNA), md)
        struct, mfe = fc.mfe()
        pairs = _pairs_from_dotbracket(struct)

        for i, ch in enumerate(struct):
            if (i < hp_start or i > hp_end) and ch in '()':
                return False, struct, float(mfe), f"pair outside hairpin at {i}"

        ok_pairs = 0
        for k in range(stem_len):
            i = hp_start + k
            j_exp = hp_end - k
            if pairs.get(i, -1) == j_exp:
                ok_pairs += 1
        if ok_pairs < stem_len - 1:
            return False, struct, float(mfe), f"weak hairpin ({ok_pairs}/{stem_len})"

        return True, struct, float(mfe), None
    except Exception as e:
        return False, None, None, f"fold error: {e}"
    
def parse_len_spec_field(val: str, field_name: str, default):
    """
    Parse a length spec from the form.
    Accepts:
      - "8"            -> 8 (int)
      - "7-10"         -> (7,10) inclusive range
      - "7,8,9"        -> [7,8,9] discrete set
    Returns default if val is empty.
    """
    if val is None or str(val).strip() == "":
        return default
    s = str(val).strip()
    try:
        if '-' in s and ',' not in s:
            a, b = s.split('-', 1)
            a, b = int(a), int(b)
            return (min(a,b), max(a,b))
        elif ',' in s:
            arr = [int(x) for x in s.split(',') if x.strip() != ""]
            if not arr:
                raise ValueError
            return arr
        else:
            return int(s)
    except Exception:
        raise ValueError(f"Invalid length spec for '{field_name}': {val}. Use forms like 8  or  7-10  or  7,8,9")

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
                                dg_sum_max=-20.0):
    """
    Iterative search for best designs; returns top_n by score.
    """
    random.seed(seed)
    n = len(miDNA)

    # Precompute cDNA (DNA alphabet) -- complement (not reverse) per your snippet
    cDNA_seq = complement(miDNA)

    # ---- 3′ hemiprobe: must sit at miRNA 3′ end, length l3_min–l3_max, ΔG in [dg_min, dg_max]
    valid3 = []
    for L3 in range(l3_min, l3_max + 1):
        if L3 > n: break
        comp3 = miDNA[n-L3:]
        dG3   = _dg(comp3) * 2
        if dg_min <= dG3 <= dg_max:
            rt3 = reverse_complement(comp3)
            valid3.append((L3, comp3, rt3, dG3))

    # ---- 5′ hemiprobe: 5′ end OR middle, length l5_min–l5_max, ΔG in [dg_min, dg_max]
    valid5_all = []
    for L5 in range(l5_min, l5_max + 1):
        if L5 >= n: break
        for s in range(0, n - L5 + 1):
            seg = miDNA[s:s+L5]
            dG5 = _dg(seg) * 2
            if dg_min <= dG5 <= dg_max:
                rt5 = reverse_complement(seg)
                valid5_all.append((s, L5, seg, rt5, dG5))

    global_heap = []  # min-heap of (-score, unique_id, cand_dict)
    seen_rt_primers = set()

    # Iterate pairs
    for (L3, comp3, rt3, dG3) in valid3:
        start_3 = n - L3
        # per-pair beam
        pair_heap = []
        pair_seq_set = set()

        for (s5, L5, comp5, rt5, dG5) in valid5_all:
            # no overlap; 3′ more negative than 5′; sum ≤ -20
            if s5 + L5 > start_3: 
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
                    rc_rt5=reverse_complement(rt5),  # equals comp5
                    L3=L3,
                    rc_rt3=reverse_complement(rt3)   # equals comp3
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
                }

                # push to per-pair beam (min-heap of size beam_per_pair)
                heapq.heappush(pair_heap, (score, cand))
                if len(pair_heap) > beam_per_pair:
                    heapq.heappop(pair_heap)

        # Merge pair beam into global heap, dedup on full RT primer
        while pair_heap:
            _, cand = heapq.heappop(pair_heap)
            rt_seq = cand["rt_primer"]
            if rt_seq in seen_rt_primers:
                continue
            seen_rt_primers.add(rt_seq)
            heapq.heappush(global_heap, (cand["score"], cand))
            if len(global_heap) > (top_n * 4):  # keep a small global pool
                heapq.heappop(global_heap)

    # Final top_n by score
    final = sorted((c for _, c in global_heap), key=lambda x: x["score"], reverse=True)
    return final[:top_n]

# ASCII hetero-dimer figure for two separated binding sites (5′ and 3′)
def heterodimer_ascii(miRNA_rna, rt5, s5, L5, rt3, L3, n):
    # miRNA (RNA letters) 5'->3'
    top = "5' " + miRNA_rna + "\n"
    # First block: 5′ hemiprobe binds positions [s5, s5+L5)
    bar5 = "    " + " " * s5 + "|" * L5 + "\n"
    bot5_seq = " " * s5 + reverse_complement(rt5).replace('T','U')  # equals miRNA segment
    bot5 = "3' " + bot5_seq + "\n\n"
    # Second block: 3′ hemiprobe binds the tail (last L3 bases)
    spaces = n - L3
    bar3 = "    " + " " * spaces + "|" * L3 + "\n"
    bot3_seq = " " * spaces + reverse_complement(rt3).replace('T','U')
    bot3 = "3' " + bot3_seq + "\n"
    # base pairs = L5 + L3 (approx)
    return top + bar5 + bot5 + top + bar3 + bot3, (L5 + L3)

def heterodimer_ascii_onepic(miRNA_rna, s5, L5, rc_rt5, L3, rc_rt3):
    """
    Single ASCII picture:
      Top:   5' <miRNA (RNA letters)>
      Mid:       pipes '|' at bound positions (5' segment and last L3)
      Bottom: 3' <DNA letters only where bound (rc_rt5 and rc_rt3), spaces elsewhere>

    rc_rt5 = reverse_complement(rt5) == comp5 (DNA)
    rc_rt3 = reverse_complement(rt3) == comp3 (DNA)
    """
    n = len(miRNA_rna)

    top_prefix = "5' "
    mid_prefix = " " * len(top_prefix)   # <- match '5' ' prefix length (3)
    bot_prefix = "3' "

    # Top (RNA)
    top = top_prefix + miRNA_rna + "\n"

    # Middle (pipes)
    bars = [" "] * n
    for i in range(s5, s5+L5): bars[i] = "|"
    for i in range(n - L3, n): bars[i] = "|"
    mid = mid_prefix + "".join(bars) + "\n"

    # Bottom (DNA only at bound positions)
    dna = [" "] * n
    dna[s5:s5+L5] = list(rc_rt5)       # DNA bases aligned to 5′ binding
    dna[n-L3:n]   = list(rc_rt3)       # DNA bases aligned to 3′ binding
    bot = bot_prefix + "".join(dna) + "\n"

    return top + mid + bot, (L5 + L3)

# ---------------- miRNA DB ----------------
def load_mirna_database():
    mirna_db = {}
    filepath = os.path.join(app.root_path, 'mirna_database.txt')
    if not os.path.exists(filepath):
        return mirna_db
    with open(filepath, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    current = None
    for line in lines:
        if line.startswith('>'):
            parts = line[1:].split()
            conventional = parts[0]
            mimat = parts[1] if len(parts) > 1 else ''
            if '-' in conventional:
                organism_code, mirna_short = conventional.split('-', 1)
            else:
                organism_code, mirna_short = '', conventional
            organism = ' '.join(parts[2:-1]) if len(parts) >= 4 else ''
            full_name = parts[-1] if len(parts) >= 3 else ''
            current = {
                'conventional_name': conventional,
                'organism_code': organism_code,
                'miRNA_name': mirna_short,
                'mimat_id': mimat,
                'organism': organism,
                'full_name': full_name,
                'sequence': ''
            }
        else:
            if current:
                current['sequence'] = line.upper()  # keep RNA (U) in DB
                mirna_db[current['conventional_name']] = current
                current = None
    return mirna_db

miRNA_DATABASE = load_mirna_database()

# ---------------- Flask routes (single page) ----------------
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html', vienna_available=VIENNA_AVAILABLE)

@app.route('/design', methods=['POST'])
def design():
    name = (request.form.get('miRNA_name') or '').strip()
    seq_rna = (request.form.get('miRNA_sequence') or '').strip()

    # concentrations
    try:
        Na  = float(request.form.get('Na_conc','5'))/1000.0
        Mg  = float(request.form.get('Mg_conc','3'))/1000.0
        dNTPs = float(request.form.get('dNTPs_conc','0.5'))/1000.0
        C_rt   = float(request.form.get('oligo_conc_rt','0.05'))/1e6
        C_qpcr = float(request.form.get('oligo_conc_qpcr','0.4'))/1e6
        if any(x<0 for x in [Na,Mg,dNTPs,C_rt,C_qpcr]):
            return jsonify({"error":"Concentrations must be non-negative."}), 400
    except ValueError:
        return jsonify({"error":"Please enter valid numeric concentrations."}), 400

    # resolve name
    if name:
        entry = miRNA_DATABASE.get(name.replace(' ',''))
        if not entry:
            return jsonify({"error":f"miRNA '{name}' not found in database."}), 404
        seq_rna = entry['sequence']

    if not seq_rna:
        return jsonify({"error":"Provide a miRNA sequence or conventional name."}), 400
    if not re.fullmatch(r'[ACGTUacgtu]+', seq_rna):
        return jsonify({"error":"Sequence must contain only A,C,G,U/T."}), 400

    # Advanced: ΔG window, seed, outputs
    try:
        dg_min = float(request.form.get('dg_min','-14'))
        dg_max = float(request.form.get('dg_max','-8'))
    except ValueError:
        return jsonify({"error":"ΔG window must be numeric."}), 400
    if dg_min > dg_max:
        dg_min, dg_max = dg_max, dg_min

        # --- Advanced numeric knobs (safe defaults)
    try:
        # ΔG temperature (for _dg): default 37 °C (your T is 310.15 K above)
        dg_temp_C = float(request.form.get('dg_temp_C', '37'))
        # Set global T used by _dg (kcal model uses T in Kelvin)
        global T
        T = 273.15 + dg_temp_C

        # Iterations / search size / output
        tries_per_pair = int(request.form.get('tries_per_pair', '120'))
        beam_per_pair  = int(request.form.get('beam_per_pair', '3'))
        top_n          = int(request.form.get('top_n', '5'))

        # Hemiprobe length ranges
        l3_min = int(request.form.get('l3_min', '4'))
        l3_max = int(request.form.get('l3_max', '8'))
        l5_min = int(request.form.get('l5_min', '4'))
        l5_max = int(request.form.get('l5_max', '8'))
        if l3_min > l3_max: l3_min, l3_max = l3_max, l3_min
        if l5_min > l5_max: l5_min, l5_max = l5_max, l5_min

        # ΔG sum cutoff
        dg_sum_max = float(request.form.get('dg_sum_max', '-20'))

        # qPCR Tm window and length bounds
        tm_lo = float(request.form.get('tm_lo', '58'))
        tm_hi = float(request.form.get('tm_hi', '65'))
        tm_window = (min(tm_lo, tm_hi), max(tm_lo, tm_hi))
        max_primer_len = int(request.form.get('max_primer_len', '25'))
        min_primer_len = int(request.form.get('min_primer_len', '16'))

        # Hairpin and random spacer length specs (accept 8  or  7-9  or  8,9)
        stem_len  = parse_len_spec_field(request.form.get('stem_len',  '7-8'),  'stem_len',  (7,8))
        loop_len  = parse_len_spec_field(request.form.get('loop_len',  '5'),    'loop_len',  (5,5))
        rand5_len = parse_len_spec_field(request.form.get('rand5_len', '8-10'), 'rand5_len', (8,10))
        rand3_len = parse_len_spec_field(request.form.get('rand3_len', '8-10'), 'rand3_len', (8,10))
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    # seed: if "randomize_seed" checked, create cryptographically-strong random seed
    randomize_seed = get_bool(request.form, 'randomize_seed', False)
    if randomize_seed:
        import secrets
        seed = secrets.randbits(8)
    else:
        try:
            seed = int(request.form.get('seed','1'))
        except ValueError:
            seed = 1

    include_rt_product = get_bool(request.form, 'include_rt_product', False)
    include_score      = get_bool(request.form, 'include_score', False)
    include_ascii      = get_bool(request.form, 'include_ascii', False)
    include_struct     = get_bool(request.form, 'include_struct', False)

    randomize_seed = get_bool(request.form, 'randomize_seed', False)
    auto_relax     = get_bool(request.form, 'auto_relax', False)

    miDNA = seq_rna.replace('U','T').upper()
    cDNA  = reverse_complement(miDNA)

    # search (with optional auto-relax)
    candidates = enumerate_strict_candidates(
        miDNA, Na, Mg, dNTPs, C_rt, C_qpcr,
        dg_min=dg_min, dg_max=dg_max,
        top_n=top_n, seed=seed,
        tries_per_pair=tries_per_pair,
        beam_per_pair=beam_per_pair,
        stem_len=stem_len, loop_len=loop_len,
        rand5_len=rand5_len, rand3_len=rand3_len,
        tm_window=tm_window,
        max_primer_len=max_primer_len,
        min_primer_len=min_primer_len,
        l3_min=l3_min, l3_max=l3_max,
        l5_min=l5_min, l5_max=l5_max,
        dg_sum_max=dg_sum_max
    )

    if not candidates and auto_relax:
        for slack in (0.5, 1.0, 1.5, 2.0):
            candidates = enumerate_strict_candidates(
                miDNA, Na, Mg, dNTPs, C_rt, C_qpcr,
                dg_min=dg_min - slack, dg_max=dg_max + slack,
                top_n=top_n, seed=seed,
                tries_per_pair=tries_per_pair,
                beam_per_pair=beam_per_pair,
                stem_len=stem_len, loop_len=loop_len,
                rand5_len=rand5_len, rand3_len=rand3_len,
                tm_window=tm_window,
                max_primer_len=max_primer_len,
                min_primer_len=min_primer_len,
                l3_min=l3_min, l3_max=l3_max,
                l5_min=l5_min, l5_max=l5_max,
                dg_sum_max=dg_sum_max
            )
            if candidates:
                break

    # Filter fields
    filtered = []
    for c in candidates:
        keep = {
            # always keep these core fields
            "colored_html": c.get("colored_html",""),
            "RT_primer_Tm": c.get("RT_primer_Tm","N/A"),
            "RT_primer_GC": c.get("RT_primer_GC","N/A"),
            "Forward_primer": c.get("Forward_primer",""),
            "Forward_primer_Tm": c.get("Forward_primer_Tm","N/A"),
            "Reverse_primer": c.get("Reverse_primer",""),
            "Reverse_primer_Tm": c.get("Reverse_primer_Tm","N/A"),
            "dg3": c.get("dg3","N/A"),
            "dg5": c.get("dg5","N/A"),
            "dg_sum": c.get("dg_sum","N/A"),
        }
        if include_rt_product:
            keep["rt_product"] = c.get("rt_product","")
        if include_score:
            sc = c.get("score", None)
            keep["score"] = round(float(sc), 2) if isinstance(sc,(int,float)) else sc
        if include_ascii:
            keep["heterodimer_ascii"] = c.get("heterodimer_ascii","")
            keep["heterodimer_pairs"] = c.get("heterodimer_pairs","")
        if include_struct:
            keep["rt_struct_42C"] = c.get("rt_struct_42C","N/A")
            keep["rt_mfe_42C"]    = c.get("rt_mfe_42C","N/A")
        filtered.append(keep)

    return jsonify({
        "miRNA_name": name or "Custom miRNA",
        "miRNA_sequence": seq_rna,
        "cDNA_sequence": cDNA,
        "vienna_available": VIENNA_AVAILABLE,
        "dg_window": [dg_min, dg_max],
        "seed_used": seed,
        "candidates": filtered,
        "advanced_used": {
            "dg_temp_C": dg_temp_C,
            "dg_sum_max": dg_sum_max,
            "tries_per_pair": tries_per_pair,
            "beam_per_pair": beam_per_pair,
            "top_n": top_n,
            "l3_range": [l3_min, l3_max],
            "l5_range": [l5_min, l5_max],
            "tm_window": [tm_window[0], tm_window[1]],
            "min_primer_len": min_primer_len,
            "max_primer_len": max_primer_len,
            "stem_len": str(stem_len),
            "loop_len": str(loop_len),
            "rand5_len": str(rand5_len),
            "rand3_len": str(rand3_len),
        }
    })

@app.route('/autocomplete', methods=['GET'])
def autocomplete():
    q = (request.args.get('q') or '').lower()
    suggestions = []
    if q:
        for key in miRNA_DATABASE:
            if q in key.lower():
                suggestions.append(key)
    return jsonify({"suggestions": suggestions})

if __name__ == '__main__':
    app.run(debug=True)
