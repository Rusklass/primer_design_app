# app.py

from flask import Flask, render_template, request, jsonify, session, Response
import csv
import io
import numpy as np
import random
import os
from datetime import datetime, timedelta, UTC
import secrets

app = Flask(__name__)
app.secret_key = secrets.token_hex(32)

user_history = {}  # { ip: [ (timestamp, result_dict), ... ] }
HISTORY_DURATION = timedelta(hours=6)

def get_client_ip():
    return request.headers.get('X-Forwarded-For', request.remote_addr)

def clean_old_history():
    now = datetime.now(UTC)
    for ip in list(user_history.keys()):
        user_history[ip] = [(ts, data) for ts, data in user_history[ip] if now - ts < HISTORY_DURATION]
        if not user_history[ip]:
            del user_history[ip]

# Thermodynamic constants
R = 1.987  # cal/(K*mol)
Kx = 3e4
T = 298.15  # K (25°C)

# DNA nearest-neighbor parameters
nn_params = {
    'AA': {'dH': -7.9, 'dS': -22.2}, 'TT': {'dH': -7.9, 'dS': -22.2},
    'AT': {'dH': -7.2, 'dS': -20.4}, 'TA': {'dH': -7.2, 'dS': -21.3},
    'CA': {'dH': -8.5, 'dS': -22.7}, 'TG': {'dH': -8.5, 'dS': -22.7},
    'GT': {'dH': -8.4, 'dS': -22.4}, 'AC': {'dH': -8.4, 'dS': -22.4},
    'CT': {'dH': -7.8, 'dS': -21.0}, 'AG': {'dH': -7.8, 'dS': -21.0},
    'GA': {'dH': -8.2, 'dS': -22.2}, 'TC': {'dH': -8.2, 'dS': -22.2},
    'CG': {'dH': -10.6, 'dS': -27.2}, 'GC': {'dH': -9.8, 'dS': -24.4},
    'GG': {'dH': -8.0, 'dS': -19.9}, 'CC': {'dH': -8.0, 'dS': -19.9}
}

def load_mirna_database():
    mirna_db = {}
    filepath = os.path.join(app.root_path, 'mirna_database.txt')
    with open(filepath, 'r') as file:
        lines = file.readlines()
    current_entry = None
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            parts = line[1:].split()
            conventional_name = parts[0]
            mimat_id = parts[1]
            organism_code, miRNA_name = (conventional_name.split('-', 1) if '-' in conventional_name else ('', conventional_name))
            organism = ' '.join(parts[2:-1])
            full_name = parts[-1]
            current_entry = {
                'conventional_name': conventional_name,
                'organism_code': organism_code,
                'miRNA_name': miRNA_name,
                'mimat_id': mimat_id,
                'organism': organism,
                'full_name': full_name,
                'sequence': ''
            }
        elif current_entry:
            current_entry['sequence'] = line.upper()
            mirna_db[current_entry['conventional_name']] = current_entry
            current_entry = None
    return mirna_db

miRNA_DATABASE = load_mirna_database()

# --- Utility Functions ---
def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'U':'A'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])

def reverse(seq):
    return seq[::-1]

def complement(seq):
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement_map.get(base, base) for base in seq])

def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if seq else 0

def calculate_thermodynamics(seq):
    dH, dS = 0.2, -5.7
    if seq:
        if seq[0] in ('A', 'T'): dH += 2.2; dS += 6.9
        if seq[-1] in ('A', 'T'): dH += 2.2; dS += 6.9
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        param = nn_params.get(pair) or nn_params.get(pair[::-1])
        if param:
            dH += param['dH']
            dS += param['dS']
        else:
            print(f"Warning: Unknown pair {pair}")
    return dH, dS

def calculate_tm(seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc):
    dH, dS = calculate_thermodynamics(seq)
    dH *= 1000
    oligo_conc_adjusted = oligo_conc / 4
    tm_k = dH / (dS + R * np.log(oligo_conc_adjusted))
    tm_c = tm_k - 273.15
    salt_effect = 16.6 * np.log10(Na_conc + 4 * np.sqrt(Mg_conc))
    return tm_c + salt_effect

def generate_random_nucleotides(length):
    return ''.join(random.choices('AGCT', k=length))

def has_complementarity(seq1, seq2, max_bp=4):
    seq2_rc = reverse_complement(seq2)
    for i in range(len(seq1)):
        for j in range(len(seq2_rc)):
            match_len = 0
            while i + match_len < len(seq1) and j + match_len < len(seq2_rc) and seq1[i + match_len] == seq2_rc[j + match_len]:
                match_len += 1
                if match_len > max_bp:
                    return True
    return False

def tune_primer(primer_seq, target_tm_range, oligo_conc, Na_conc, Mg_conc, dNTPs_conc, max_length=22):
    """
    Given an initial primer sequence, adjust its 5' end (by adding a GC clamp nucleotide)
    if the melting temperature is too low. Also, if Tm is too high, you might remove nucleotides.
    The final primer length must be between 20 and max_length.
    
    Returns the tuned primer sequence and its Tm.
    """
    current_seq = primer_seq
    current_tm = calculate_tm(current_seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc)
    
    # If Tm is too low, try adding a nucleotide at 5' (GC clamp)
    while current_tm < target_tm_range[0] and len(current_seq) < max_length:
        # Prepend a 'G' (or choose a nucleotide based on your design rules)
        current_seq = "G" + current_seq
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

def design_primers(miRNA_name, miRNA_seq_rna, params):
    # Extract concentrations from params
    Na_conc = params['Na_conc']
    Mg_conc = params['Mg_conc']
    dNTPs_conc = params['dNTPs_conc']
    oligo_conc_rt = params['oligo_conc_rt']
    oligo_conc_qpcr = params['oligo_conc_qpcr']
    
    miRNA_seq = miRNA_seq_rna.replace('U', 'T').upper()
    # Generate cDNA sequence
    cDNA_seq = complement(miRNA_seq)
    
    for length in range(5, len(miRNA_seq) + 1):
        comp_region_5 = miRNA_seq[:length]
        rt_primer_5 = reverse_complement(comp_region_5)
        dH, dS = calculate_thermodynamics(comp_region_5)
        dG5 = dH - (T * dS / 1000)  # Gibbs free energy at 25°C
        if -14 <= dG5 <= -10:
            break
    else:
        return {"error": "No suitable 5' complementary region found."}
    
    for length in range(5, len(miRNA_seq) - len(comp_region_5) + 1):
        comp_region_3 = miRNA_seq[-length:]
        rt_primer_3 = reverse_complement(comp_region_3)
        dH, dS = calculate_thermodynamics(comp_region_3)
        dG3 = dH - (T * dS / 1000)
        if -14 <= dG3 <= -10:
            break
    else:
        comp_region_3 = ''
        rt_primer_3 = ''
        dG3 = None

    random_nt_length_5 = random.choice([8, 9])
    random_nt_length_3 = random.choice([8, 9])
    
    for _ in range(100):
        random_nt_5 = generate_random_nucleotides(random_nt_length_5)
        if not any([
            has_complementarity(random_nt_5, rt_primer_5),
            has_complementarity(random_nt_5, rt_primer_3)
        ]):
            break
    else:
        return {"error": "Failed to generate non-complementary random_nt_5."}
    
    for _ in range(100):
        random_nt_3 = generate_random_nucleotides(random_nt_length_3)
        if not any([
            has_complementarity(random_nt_3, rt_primer_5),
            has_complementarity(random_nt_3, rt_primer_3),
            has_complementarity(random_nt_3, random_nt_5)
        ]):
            break
    else:
        return {"error": "Failed to generate non-complementary random_nt_3."}
    
    hairpin_stem = generate_random_nucleotides(8)
    hairpin_loop = generate_random_nucleotides(5)
    hairpin = hairpin_stem + hairpin_loop + reverse_complement(hairpin_stem)
    
    rt_primer = rt_primer_5 + random_nt_5 + hairpin + random_nt_3 + rt_primer_3
    rt_product = rt_primer_5 + random_nt_5 + hairpin + random_nt_3 + reverse(cDNA_seq)
    
    # --------------------------
    # New forward primer design:
    # Take candidate of 18 nt from the 5' end of rt_primer and then adjust if needed.
    candidate_forward = rt_primer[:18]
    try:
        forward_primer, fwd_tm = tune_primer(candidate_forward, (58, 65), oligo_conc_qpcr,
                                             Na_conc, Mg_conc, dNTPs_conc, max_length=25)
    except ValueError as e:
        return {"error": "Forward primer design failed: " + str(e)}
    
    # --------------------------
    # New reverse primer design:
    # Take candidate of 18 nt from the 3' end of rt_product,
    # then obtain its reverse complement and adjust if needed.
    candidate_rev_template = rt_product[-18:]
    candidate_reverse = reverse_complement(candidate_rev_template)
    try:
        reverse_primer, rev_tm = tune_primer(candidate_reverse, (58, 65), oligo_conc_qpcr,
                                             Na_conc, Mg_conc, dNTPs_conc, max_length=25)
    except ValueError as e:
        return {"error": "Reverse primer design failed: " + str(e)}
    
    # Calculate melting temperatures and GC content for the other primers as needed
    rt_tm = calculate_tm(rt_primer, oligo_conc_rt, Na_conc, Mg_conc, dNTPs_conc)
    rt_gc = calculate_gc_content(rt_primer)
    fwd_gc = calculate_gc_content(forward_primer)
    rev_gc = calculate_gc_content(reverse_primer)
    
    results = {
        "miRNA_name": miRNA_name,
        "miRNA_sequence": miRNA_seq_rna,
        "cDNA_sequence": cDNA_seq,
        "RT_primer": rt_primer,
        "RT_primer_Tm": f"{rt_tm:.2f}°C",
        "RT_primer_GC": f"{rt_gc:.2f}%",
        "RT_primer_5_region": rt_primer_5,
        "RT_primer_5_dG": f"{dG5:.2f} kcal/mol" if dG5 else "N/A",
        "RT_primer_3_region": rt_primer_3,
        "RT_primer_3_dG": f"{dG3:.2f} kcal/mol",
        "RT_product": rt_product,
        "Forward_primer": forward_primer,
        "Forward_primer_Tm": f"{fwd_tm:.2f}°C",
        "Forward_primer_GC": f"{fwd_gc:.2f}%",
        "Reverse_primer": reverse_primer,
        "Reverse_primer_Tm": f"{rev_tm:.2f}°C",
        "Reverse_primer_GC": f"{rev_gc:.2f}%",
        # Including additional RT primer components for troubleshooting if needed:
        "rt_primer_3": rt_primer_3,
        "random_nt_3": random_nt_3,
        "hairpin": hairpin,
        "random_nt_5": random_nt_5,
        "rt_primer_5": rt_primer_5,
    }
    return results


# Routes and views
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        miRNA_name = request.form.get('miRNA_name', '').strip()
        miRNA_sequence = request.form.get('miRNA_sequence', '').strip()

        try:
            Na_conc = float(request.form.get('Na_conc', '5')) / 1000
            Mg_conc = float(request.form.get('Mg_conc', '3')) / 1000
            dNTPs_conc = float(request.form.get('dNTPs_conc', '0.5')) / 1000
            oligo_conc_rt = float(request.form.get('oligo_conc_rt', '0.05')) / 1e6
            oligo_conc_qpcr = float(request.form.get('oligo_conc_qpcr', '0.4')) / 1e6
            if any(c < 0 for c in [Na_conc, Mg_conc, dNTPs_conc, oligo_conc_rt, oligo_conc_qpcr]):
                return render_template('index.html', error="Concentration values must be non-negative.")
        except ValueError:
            return render_template('index.html', error="Please enter valid numerical values for concentrations.")

        if miRNA_name:
            miRNA_name = miRNA_name.replace(' ', '')
            if miRNA_name in miRNA_DATABASE:
                miRNA_sequence = miRNA_DATABASE[miRNA_name]['sequence']
            else:
                return render_template('index.html', error=f"miRNA '{miRNA_name}' not found in the database.")

        if not miRNA_sequence:
            return render_template('index.html', error="Please provide a miRNA sequence.")

        params = {
            'Na_conc': Na_conc,
            'Mg_conc': Mg_conc,
            'dNTPs_conc': dNTPs_conc,
            'oligo_conc_rt': oligo_conc_rt,
            'oligo_conc_qpcr': oligo_conc_qpcr
        }

        results = design_primers(miRNA_name or 'Custom miRNA', miRNA_sequence, params)

        if "error" in results:
            return render_template('index.html', error=results["error"])

        ip = get_client_ip()
        now = datetime.now(UTC)
        clean_old_history()
        user_history.setdefault(ip, []).append((now, results))

        session['history'] = session.get('history', [])[-4:] + [results]
        # return render_template('index.html', 
        #                          results=results, 
        #                          history=[h[1] for h in user_history[ip]], 
        #                          current_year=datetime.now().year)
        return render_template('index.html', 
                       results=results,
                       current_year=datetime.now().year, 
                       history=[h[1] for h in user_history.get(get_client_ip(), [])])
    return render_template('index.html', current_year=datetime.now().year)

# Autocomplete route
@app.route('/autocomplete', methods=['GET'])
def autocomplete():
    query = request.args.get('q', '').lower()
    suggestions = []
    if query:
        for key in miRNA_DATABASE:
            if query in key.lower():
                suggestions.append(key)
    return jsonify({'suggestions': suggestions})

@app.route('/download_csv')
def download_csv():
    ip = get_client_ip()
    clean_old_history()
    history = user_history.get(ip, [])[-10:]  # Last 10 entries

    # If no history, show empty CSV
    if not history:
        output = "No data available.\n"
        return Response(output, mimetype="text/csv", headers={"Content-Disposition": "attachment;filename=primer_history.csv"})

    # Create CSV in memory
    output = io.StringIO()
    writer = csv.writer(output)

    # Write header
    writer.writerow([
        "miRNA Name", "Forward Primer", "Reverse Primer",
        "RT Primer", "RT Primer Tm", "RT Primer GC",
        "Forward Primer Tm", "Forward Primer GC",
        "Reverse Primer Tm", "Reverse Primer GC"
    ])

    # Write data rows
    for _, result in history:
        writer.writerow([
            result.get("miRNA_name", ""),
            result.get("Forward_primer", ""),
            result.get("Reverse_primer", ""),
            result.get("RT_primer", ""),
            result.get("RT_primer_Tm", ""),
            result.get("RT_primer_GC", ""),
            result.get("Forward_primer_Tm", ""),
            result.get("Forward_primer_GC", ""),
            result.get("Reverse_primer_Tm", ""),
            result.get("Reverse_primer_GC", "")
        ])

    output.seek(0)
    return Response(output, mimetype="text/csv", headers={"Content-Disposition": "attachment;filename=primer_history.csv"})

if __name__ == '__main__':
    app.run(debug=True)
