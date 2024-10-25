# app.py

from flask import Flask, render_template, request, jsonify
import numpy as np
import random
import os

app = Flask(__name__)

# Define constants
R = 1.987  # Ideal gas constant in cal/(K*mol)
Kx = 3e4   # Equilibrium constant for Mg2+ and dNTPs binding
T = 298.15 # Temperature for Gibbs energy calculation (25°C + 273.15°K)

# Nearest-neighbor thermodynamic parameters for DNA/DNA duplexes
nn_params = {
    'AA': {'dH': -7.9, 'dS': -22.2},
    'TT': {'dH': -7.9, 'dS': -22.2},
    'AT': {'dH': -7.2, 'dS': -20.4},
    'TA': {'dH': -7.2, 'dS': -21.3},
    'CA': {'dH': -8.5, 'dS': -22.7},
    'TG': {'dH': -8.5, 'dS': -22.7},
    'GT': {'dH': -8.4, 'dS': -22.4},
    'AC': {'dH': -8.4, 'dS': -22.4},
    'CT': {'dH': -7.8, 'dS': -21.0},
    'AG': {'dH': -7.8, 'dS': -21.0},
    'GA': {'dH': -8.2, 'dS': -22.2},
    'TC': {'dH': -8.2, 'dS': -22.2},
    'CG': {'dH': -10.6, 'dS': -27.2},
    'GC': {'dH': -9.8, 'dS': -24.4},
    'GG': {'dH': -8.0, 'dS': -19.9},
    'CC': {'dH': -8.0, 'dS': -19.9},
}

# Function to load miRNA database
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
            # Extract organism code and miRNA name
            if '-' in conventional_name:
                organism_code, miRNA_name = conventional_name.split('-', 1)
            else:
                organism_code = ''
                miRNA_name = conventional_name
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
            current_entry['sequence'] = line.upper()  # Keep RNA format
            key = current_entry['conventional_name']
            mirna_db[key] = current_entry
            current_entry = None
    return mirna_db

# Load miRNA database at startup
miRNA_DATABASE = load_mirna_database()

# Function definitions (same as provided code)
def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_comp = ''.join([complement.get(base, base) for base in seq[::-1]])
    return rev_comp

def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def calculate_thermodynamics(seq):
    dH = 0.0  # Enthalpy in kcal/mol
    dS = 0.0  # Entropy in cal/(K*mol)
    # Initiation parameters
    init_dH = 0.2
    init_dS = -5.7
    dH += init_dH
    dS += init_dS
    # Terminal corrections for AT base pairs
    if seq and seq[0] in ('A', 'T'):
        dH += 2.2
        dS += 6.9
    if seq and seq[-1] in ('A', 'T'):
        dH += 2.2
        dS += 6.9
    # Sum nearest-neighbor parameters
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        if pair in nn_params:
            dH += nn_params[pair]['dH']
            dS += nn_params[pair]['dS']
        else:
            rev_pair = pair[::-1]
            if rev_pair in nn_params:
                dH += nn_params[rev_pair]['dH']
                dS += nn_params[rev_pair]['dS']
            else:
                print(f"Warning: Thermodynamic parameters not found for pair {pair}")
    return dH, dS

def calculate_tm(seq, oligo_conc, Na_conc, Mg_conc, dNTPs_conc):
    dH, dS = calculate_thermodynamics(seq)
    dH *= 1000  # Convert from kcal/mol to cal/mol
    # Adjust oligo concentration for non-self-complementary sequences
    oligo_conc_adjusted = oligo_conc / 4
    # Calculate Tm at standard conditions (1 M Na+)
    tm_kelvin = dH / (dS + R * (np.log(oligo_conc_adjusted)))
    tm_celsius = tm_kelvin - 273.15
    # Simplified salt corrections
    salt_effect = 16.6 * np.log10(Na_conc + 4 * np.sqrt(Mg_conc))
    tm_corrected = tm_celsius + salt_effect
    return tm_corrected

def generate_random_nucleotides(length):
    return ''.join(random.choice('AGCT') for _ in range(length))

def design_primers(miRNA_name, miRNA_seq_rna, params):
    # Extract concentrations from params
    Na_conc = params['Na_conc']
    Mg_conc = params['Mg_conc']
    dNTPs_conc = params['dNTPs_conc']
    oligo_conc_rt = params['oligo_conc_rt']
    oligo_conc_qpcr = params['oligo_conc_qpcr']

    miRNA_seq = miRNA_seq_rna.replace('U', 'T').upper()
    # Generate cDNA sequence
    cDNA_seq = reverse_complement(miRNA_seq)
    # Design RT primer
    # 3' end complementary region
    for length in range(6, len(miRNA_seq)+1):
        comp_region_3 = miRNA_seq[-length:]
        rt_primer_3 = reverse_complement(comp_region_3)
        dH, dS = calculate_thermodynamics(comp_region_3)
        dG = dH - (T * dS / 1000)  # Gibbs free energy at 25°C
        if -14 <= dG <= -10:
            break
    else:
        return {"error": "No suitable 3' complementary region found."}
    # 5' end complementary region (non-overlapping)
    for length in range(6, len(miRNA_seq) - len(comp_region_3) + 1):
        comp_region_5 = miRNA_seq[:length]
        if miRNA_seq.find(comp_region_3) < miRNA_seq.find(comp_region_5) + len(comp_region_5):
            continue
        rt_primer_5 = reverse_complement(comp_region_5)
        dH, dS = calculate_thermodynamics(comp_region_5)
        dG5 = dH - (T * dS / 1000)  # Gibbs free energy at 25°C
        if -14 <= dG5 <= -10:
            break
    else:
        comp_region_5 = ''
        rt_primer_5 = ''
        dG5 = None
    # Generate random nucleotides (8-9 bases) on both sides
    random_nt_length_5 = random.choice([8, 9])
    random_nt_length_3 = random.choice([8, 9])
    random_nt_5 = generate_random_nucleotides(random_nt_length_5)
    random_nt_3 = generate_random_nucleotides(random_nt_length_3)
    # Generate hairpin structure with 8 bp stem and 5-base loop
    hairpin_stem = generate_random_nucleotides(8)
    hairpin_loop = generate_random_nucleotides(5)
    hairpin = hairpin_stem + hairpin_loop + reverse_complement(hairpin_stem)
    # Construct RT primer
    rt_primer = rt_primer_3 + random_nt_3 + hairpin + random_nt_5 + rt_primer_5
    rt_product = rt_primer_3 + random_nt_3 + hairpin + random_nt_5 + cDNA_seq
    # Design forward primer (complementary to 5' end of RT primer)
    forward_primer = rt_primer[:20]
    # Design reverse primer (complementary to miRNA) #CHANGE/UPDATE
    reverse_primer_seq = miRNA_seq[:len(forward_primer)]
    reverse_primer = reverse_complement(reverse_primer_seq)
    # Calculate melting temperatures and GC content
    rt_tm = calculate_tm(rt_primer, oligo_conc_rt, Na_conc, Mg_conc, dNTPs_conc)
    rt_gc = calculate_gc_content(rt_primer)
    fwd_tm = calculate_tm(forward_primer, oligo_conc_qpcr, Na_conc, Mg_conc, dNTPs_conc)
    fwd_gc = calculate_gc_content(forward_primer)
    rev_tm = calculate_tm(reverse_primer, oligo_conc_qpcr, Na_conc, Mg_conc, dNTPs_conc)
    rev_gc = calculate_gc_content(reverse_primer)
    # Prepare results
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
        "RT_primer_3_dG": f"{dG:.2f} kcal/mol",
        "RT_product": rt_product,
        "Forward_primer": forward_primer,
        "Forward_primer_Tm": f"{fwd_tm:.2f}°C",
        "Forward_primer_GC": f"{fwd_gc:.2f}%",
        "Reverse_primer": reverse_primer,
        "Reverse_primer_Tm": f"{rev_tm:.2f}°C",
        "Reverse_primer_GC": f"{rev_gc:.2f}%",
        # Include the RT primer components
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

        # Retrieve concentration inputs
        try:
            Na_conc = float(request.form.get('Na_conc', '5')) / 1000  # Convert mM to M
            Mg_conc = float(request.form.get('Mg_conc', '3')) / 1000  # Convert mM to M
            dNTPs_conc = float(request.form.get('dNTPs_conc', '0.5')) / 1000  # Convert mM to M
            oligo_conc_rt = float(request.form.get('oligo_conc_rt', '0.05')) / 1e6  # Convert µM to M
            oligo_conc_qpcr = float(request.form.get('oligo_conc_qpcr', '0.4')) / 1e6  # Convert µM to M

            # Additional validation for non-negative values
            if any(conc < 0 for conc in [Na_conc, Mg_conc, dNTPs_conc, oligo_conc_rt, oligo_conc_qpcr]):
                error = "Concentration values must be non-negative."
                return render_template('index.html', error=error)
        except ValueError:
            error = "Please enter valid numerical values for concentrations."
            return render_template('index.html', error=error)

        # Existing miRNA sequence retrieval logic
        if miRNA_name:
            miRNA_name = miRNA_name.replace(' ', '')
            if miRNA_name in miRNA_DATABASE:
                miRNA_entry = miRNA_DATABASE[miRNA_name]
                miRNA_sequence = miRNA_entry['sequence']
            else:
                error = f"miRNA '{miRNA_name}' not found in the database."
                return render_template('index.html', error=error)
        if not miRNA_sequence:
            error = "Please provide a miRNA sequence."
            return render_template('index.html', error=error)

        # Prepare parameters dictionary
        params = {
            'Na_conc': Na_conc,
            'Mg_conc': Mg_conc,
            'dNTPs_conc': dNTPs_conc,
            'oligo_conc_rt': oligo_conc_rt,
            'oligo_conc_qpcr': oligo_conc_qpcr
        }

        results = design_primers(miRNA_name or 'Custom miRNA', miRNA_sequence, params)
        if "error" in results:
            error = results["error"]
            return render_template('index.html', error=error)
        return render_template('result.html', results=results)
    return render_template('index.html')

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

# Date in home route
@app.route('/')
def home():
    current_year = datetime.now().year  # Get the current year
    return render_template('index.html', current_year=current_year)

if __name__ == '__main__':
    app.run(debug=True)
