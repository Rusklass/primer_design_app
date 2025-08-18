# app.py

from flask import Flask, render_template, request, jsonify
import numpy as np
import random
import os
import re

app = Flask(__name__)

# Optional: ViennaRNA (skip gracefully if not present)
try:
    import RNA  # pip install ViennaRNA (may need system package)
    VIENNA_AVAILABLE = True
except Exception:
    VIENNA_AVAILABLE = False

# ---------------- Thermo constants ----------------
R = 1.987    # cal/(K*mol)
T = 298.15   # K (25°C) for ΔG window checks

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

def dg_at_25C(seq):
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

def generate_random_nucleotides(n):
    return ''.join(random.choice('AGCT') for _ in range(n))

def make_hairpin(stem_len=8, loop_len=5):
    def balanced_random(n):
        s=''
        while len(s)<n:
            b=random.choice('AAGCTGCT')  # light GC bias, avoids too many homopolymers
            if len(s)>=3 and s[-3:]==b*3: continue
            s+=b
        return s
    stem = balanced_random(stem_len)
    loop = balanced_random(loop_len)
    return stem + loop + reverse_complement(stem)

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

def tune_primer(primer_seq, target_tm_range, oligo_conc, Na_conc, Mg_conc, dNTPs_conc, max_length=22):
    """
    Given an initial primer sequence, adjust its 5' end (by adding a GC clamp nucleotide)
    if the melting temperature is too low. Also, if Tm is too high, you might remove nucleotides.
    The final primer length must be between 20 and max_length.
    
    Returns the tuned primer sequence and its Tm.
    """
    current_seq = primer_seq
    current_tm = calculate_tm(current_seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc)
    
    # If Tm is too low, try prepending a clamp nucleotide at 5'
    while current_tm < target_tm_range[0] and len(current_seq) < max_length:
        clamp = random.choice(['G','C','A'])  # G/C preferred, A as gentle bump
        current_seq = clamp + current_seq
        current_tm = calculate_tm(current_seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc)
    
    # Alternatively, if Tm is too high, you might remove nucleotides (if length > 20)
    while current_tm > target_tm_range[1] and len(current_seq) > 20:
        # Remove the first nucleotide and recalc Tm
        current_seq = current_seq[1:]
        current_tm = calculate_tm(current_seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc)
    
    # Check that final length is within the allowed range
    if not (20 <= len(current_seq) <= max_length) or not (target_tm_range[0] <= current_tm <= target_tm_range[1]):
        raise ValueError("Could not adjust primer to meet Tm and length criteria")
    
    return current_seq, current_tm

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

# Windows exactly per spec: around -10 (accept -14..-9)
DG_MIN, DG_MAX = -15.0, -8.0

def enumerate_strict_candidates(miDNA, Na, Mg, dNTPs, C_rt, C_qpcr,
                                dg_min=-14.0, dg_max=-9.0,
                                top_n=5, seed=13):
    random.seed(seed)
    n = len(miDNA)

    # Precompute cDNA (DNA alphabet)
    cDNA_seq = complement(miDNA)

    # ---- 3′ hemiprobe: must sit at miRNA 3′ end, length 6–10, ΔG in [dg_min, dg_max]
    valid3 = []
    for L3 in range(6, 11):
        if L3 > n: break
        comp3 = miDNA[n-L3:]
        dG3   = dg_at_25C(comp3)
        if dg_min <= dG3 <= dg_max:
            rt3 = reverse_complement(comp3)  # RT primer 3′ end
            valid3.append((L3, comp3, rt3, dG3))

    # ---- 5′ hemiprobe: 5′ end OR middle, length 6–10, ΔG in [dg_min, dg_max]
    valid5_all = []
    for L5 in range(6, 11):
        if L5 >= n: break
        for s in range(0, n - L5 + 1):
            seg = miDNA[s:s+L5]
            dG5 = dg_at_25C(seg)
            if dg_min <= dG5 <= dg_max:
                rt5 = reverse_complement(seg)   # RT primer 5′ end
                valid5_all.append((s, L5, seg, rt5, dG5))

    candidates = []
    for (L3, comp3, rt3, dG3) in valid3:
        start_3 = n - L3
        for (s5, L5, comp5, rt5, dG5) in valid5_all:
            # no overlap; 3′ more negative than 5′; sum ≤ -20
            if s5 + L5 > start_3: continue
            if not (dG3 <= dG5 and (dG3 + dG5) <= -20.0): continue

            # Build RT primer: 5′ [rt5] + rand5 + hairpin + rand3 + [rt3] 3′
            hairpin = make_hairpin(8, 5)
            rand5 = generate_random_nucleotides(random.choice([8, 9]))
            rand3 = generate_random_nucleotides(random.choice([8, 9]))
            rt_primer = rt5 + rand5 + hairpin + rand3 + rt3

            # RT reaction product template (what reverse primer sees): RT primer + reverse(cDNA)
            rt_product = rt5 + rand5 + hairpin + rand3 + reverse(cDNA_seq)

            # Forward primer: start at 5′ end of RT primer; tune to 58–65 °C
            candidate_forward = rt_primer[:18]
            try:
                fwd, f_tm = tune_primer(candidate_forward, (58, 65), C_qpcr, Na, Mg, dNTPs, max_length=25)
            except ValueError:
                continue  # try next candidate

            # Reverse primer: last 18 nt of RT product, RC, then tune
            candidate_rev_template = rt_product[-18:]
            candidate_reverse = reverse_complement(candidate_rev_template)
            try:
                rev, r_tm = tune_primer(candidate_reverse, (58, 65), C_qpcr, Na, Mg, dNTPs, max_length=25)
            except ValueError:
                continue  # try next candidate

            # RT primer properties
            rt_tm = calculate_tm(rt_primer, C_rt, Na, Mg, dNTPs)
            rt_gc = calculate_gc_content(rt_primer)

            # Hairpin-only check at 42 °C (optional, requires ViennaRNA)
            rt_struct, rt_mfe = vienna_fold_single(rt_primer, 42.0)

            # Colored HTML parts (same palette)
            colored = (
                f'<span style="color:#6A1B9A;">{rt5}</span>'
                f'<span style="color:#F9A825;">{rand5}</span>'
                f'<span style="color:#1E88E5;">{hairpin}</span>'
                f'<span style="color:#1B5E20;">{rand3}</span>'
                f'<span style="color:#D81B60;">{rt3}</span>'
            )

            # One-picture ASCII (RNA top, DNA bottom) with aligned pipes
            ascii_fig, basepairs = heterodimer_ascii_onepic(
                miRNA_rna=miDNA.replace('T','U'),     # top line RNA
                s5=s5, L5=L5,
                rc_rt5=reverse_complement(rt5),       # DNA, equals comp5
                L3=L3,
                rc_rt3=reverse_complement(rt3)        # DNA, equals comp3
            )

            # Score
            score = score_design(rt_tm=rt_tm, f_tm=f_tm, r_tm=r_tm, gc=rt_gc,
                                dG3=dG3, dG5=dG5, length=len(rt_primer))

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
            candidates.append(cand)

    candidates.sort(key=lambda c: c["score"], reverse=True)
    return candidates[:top_n]

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
        dg_max = float(request.form.get('dg_max','-9'))
    except ValueError:
        return jsonify({"error":"ΔG window must be numeric."}), 400
    if dg_min > dg_max:
        dg_min, dg_max = dg_max, dg_min

    # seed: if "randomize_seed" checked, create cryptographically-strong random seed
    randomize_seed = get_bool(request.form, 'randomize_seed', False)
    if randomize_seed:
        import secrets
        seed = secrets.randbits(32)
    else:
        try:
            seed = int(request.form.get('seed','13'))
        except ValueError:
            seed = 13

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
        top_n=5, seed=seed
    )

    if not candidates and auto_relax:
        for slack in (0.5, 1.0, 1.5, 2.0):
            candidates = enumerate_strict_candidates(
                miDNA, Na, Mg, dNTPs, C_rt, C_qpcr,
                dg_min=dg_min - slack, dg_max=dg_max + slack,
                top_n=5, seed=seed
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
        "candidates": filtered
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
