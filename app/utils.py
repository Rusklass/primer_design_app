import random
import re
import os
import numpy as np
from typing import Optional, List, Dict, Tuple

# Optional: ViennaRNA (skip gracefully if not present)
try:
    import RNA  # pip install ViennaRNA (may need system package)
    VIENNA_AVAILABLE = True
except Exception:
    VIENNA_AVAILABLE = False

# ---------------- Thermo constants ----------------
R = 1.987    # cal/(K*mol)
DEFAULT_T = 310.15   # K (37°C) for ΔG window checks

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

def _dg(seq, temp_k=DEFAULT_T):
    dH, dS = _nn_thermo(seq)
    return dH - (temp_k * dS / 1000.0)

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

def heterodimer_ascii_onepic(miRNA_rna, s5, L5, rt5, L3, rt3):
    """
    Single ASCII picture:
      Top:   5' <miRNA (RNA letters)>
      Mid:       pipes '|' at bound positions (5' segment and last L3)
      Bottom: 3' <DNA in antiparallel orientation (reversed)>

    rt5 = DNA primer segment binding to 5' region (stored 5'→3')
    rt3 = DNA primer segment binding to 3' region (stored 5'→3')
    
    DNA is shown in 3'→5' orientation (reversed) to represent antiparallel binding.
    """
    n = len(miRNA_rna)

    top_prefix = "5' "
    mid_prefix = " " * len(top_prefix)
    bot_prefix = "3' "

    # Top (RNA)
    top = top_prefix + miRNA_rna + "\n"

    # Middle (pipes)
    bars = [" "] * n
    for i in range(s5, s5+L5): bars[i] = "|"
    for i in range(n - L3, n): bars[i] = "|"
    mid = mid_prefix + "".join(bars) + "\n"

    # Bottom (DNA in antiparallel orientation - reversed)
    dna = [" "] * n
    dna[s5:s5+L5] = list(reverse(rt5))  # Show rt5 reversed (3'→5')
    dna[n-L3:n] = list(reverse(rt3))    # Show rt3 reversed (3'→5')
    bot = bot_prefix + "".join(dna) + "\n"

    return top + mid + bot, (L5 + L3)

def load_mirna_database(root_path):
    mirna_db = {}
    filepath = os.path.join(root_path, 'mirna_database.txt')
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
