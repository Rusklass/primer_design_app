import os
import csv
import io
from flask import Blueprint, render_template, request, jsonify, current_app, Response
from .utils import (
    load_mirna_database, get_bool, parse_len_spec_field, 
    reverse_complement, VIENNA_AVAILABLE
)
from .primer_design import enumerate_strict_candidates

main = Blueprint('main', __name__)

# Global DB cache
miRNA_DATABASE = {}

def get_mirna_db():
    global miRNA_DATABASE
    if not miRNA_DATABASE:
        # Assuming mirna_database.txt is in the project root (one level up from app/)
        root_path = os.path.dirname(current_app.root_path)
        miRNA_DATABASE = load_mirna_database(root_path)
    return miRNA_DATABASE

@main.route('/', methods=['GET'])
def index():
    return render_template('index.html', vienna_available=VIENNA_AVAILABLE)

@main.route('/autocomplete', methods=['GET'])
def autocomplete():
    q = (request.args.get('q') or '').lower()
    db = get_mirna_db()
    suggestions = []
    if q:
        for key in db:
            if q in key.lower():
                suggestions.append(key)
    return jsonify({"suggestions": suggestions})

def run_design(form_data):
    name = (form_data.get('miRNA_name') or '').strip()
    seq_rna = (form_data.get('miRNA_sequence') or '').strip()

    # concentrations
    try:
        Na  = float(form_data.get('Na_conc','5'))/1000.0
        Mg  = float(form_data.get('Mg_conc','3'))/1000.0
        dNTPs = float(form_data.get('dNTPs_conc','0.5'))/1000.0
        C_rt   = float(form_data.get('oligo_conc_rt','0.05'))/1e6
        C_qpcr = float(form_data.get('oligo_conc_qpcr','0.4'))/1e6
        if any(x<0 for x in [Na,Mg,dNTPs,C_rt,C_qpcr]):
            return {"error":"Concentrations must be non-negative."}, 400
    except ValueError:
        return {"error":"Please enter valid numeric concentrations."}, 400

    # resolve name
    if name:
        db = get_mirna_db()
        entry = db.get(name.replace(' ',''))
        if not entry:
            return {"error":f"miRNA '{name}' not found in database."}, 404
        seq_rna = entry['sequence']

    if not seq_rna:
        return {"error":"Provide a miRNA sequence or conventional name."}, 400
    # Basic validation
    import re
    if not re.fullmatch(r'[ACGTUacgtu]+', seq_rna):
        return {"error":"Sequence must contain only A,C,G,U/T."}, 400

    # Advanced: ΔG window, seed, outputs
    try:
        dg_min = float(form_data.get('dg_min','-14'))
        dg_max = float(form_data.get('dg_max','-8'))
    except ValueError:
        return {"error":"ΔG window must be numeric."}, 400
    if dg_min > dg_max:
        dg_min, dg_max = dg_max, dg_min

    try:
        dg_temp_C = float(form_data.get('dg_temp_C', '37'))
        
        tries_per_pair = int(form_data.get('tries_per_pair', '120'))
        beam_per_pair  = int(form_data.get('beam_per_pair', '3'))
        top_n          = int(form_data.get('top_n', '5'))

        l3_min = int(form_data.get('l3_min', '4'))
        l3_max = int(form_data.get('l3_max', '8'))
        l5_min = int(form_data.get('l5_min', '4'))
        l5_max = int(form_data.get('l5_max', '8'))
        if l3_min > l3_max: l3_min, l3_max = l3_max, l3_min
        if l5_min > l5_max: l5_min, l5_max = l5_max, l5_min

        p5_min = form_data.get('p5_min')
        p5_min = int(p5_min) if p5_min and p5_min.strip() else None
        p5_max = form_data.get('p5_max')
        p5_max = int(p5_max) if p5_max and p5_max.strip() else None

        p3_min = form_data.get('p3_min')
        p3_min = int(p3_min) if p3_min and p3_min.strip() else None
        p3_max = form_data.get('p3_max')
        p3_max = int(p3_max) if p3_max and p3_max.strip() else None

        dg_sum_max = float(form_data.get('dg_sum_max', '-20'))

        tm_lo = float(form_data.get('tm_lo', '58'))
        tm_hi = float(form_data.get('tm_hi', '65'))
        tm_window = (min(tm_lo, tm_hi), max(tm_lo, tm_hi))
        max_primer_len = int(form_data.get('max_primer_len', '25'))
        min_primer_len = int(form_data.get('min_primer_len', '16'))

        stem_len  = parse_len_spec_field(form_data.get('stem_len',  '7-8'),  'stem_len',  (7,8))
        loop_len  = parse_len_spec_field(form_data.get('loop_len',  '5'),    'loop_len',  (5,5))
        rand5_len = parse_len_spec_field(form_data.get('rand5_len', '8-10'), 'rand5_len', (8,10))
        rand3_len = parse_len_spec_field(form_data.get('rand3_len', '8-10'), 'rand3_len', (8,10))
    except ValueError as e:
        return {"error": str(e)}, 400

    randomize_seed = get_bool(form_data, 'randomize_seed', False)
    if randomize_seed:
        import secrets
        seed = secrets.randbits(8)
    else:
        try:
            seed = int(form_data.get('seed','1'))
        except ValueError:
            seed = 1

    auto_relax     = get_bool(form_data, 'auto_relax', False)
    
    # Parse variant sequences (newline separated)
    variants_raw = (form_data.get('variant_sequences') or '').strip()
    variant_sequences = [v.strip().replace('U','T').upper() for v in variants_raw.split('\n') if v.strip()]

    miDNA = seq_rna.replace('U','T').upper()
    cDNA  = reverse_complement(miDNA)

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
        dg_sum_max=dg_sum_max,
        dg_temp_C=dg_temp_C,
        variant_sequences=variant_sequences,
        p5_range=(p5_min, p5_max),
        p3_range=(p3_min, p3_max)
    )

    # Auto-relax if we don't have enough candidates
    if auto_relax and len(candidates) < top_n:
        for slack in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
            # If we have some candidates, keep them? 
            # Actually, enumerate_strict_candidates returns a fresh list. 
            # We might want to merge or just take the relaxed set if it's better.
            # Simple approach: just try to get a full set with relaxed constraints.
            
            new_candidates = enumerate_strict_candidates(
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
                dg_sum_max=dg_sum_max + slack, # also relax sum
                dg_temp_C=dg_temp_C,
                variant_sequences=variant_sequences,
                p5_range=(p5_min, p5_max),
                p3_range=(p3_min, p3_max)
            )
            
            if len(new_candidates) > len(candidates):
                candidates = new_candidates
            
            if len(candidates) >= top_n:
                break
    
    return {
        "miRNA_name": name or "Custom miRNA",
        "miRNA_sequence": seq_rna,
        "cDNA_sequence": cDNA,
        "vienna_available": VIENNA_AVAILABLE,
        "dg_window": [dg_min, dg_max],
        "seed_used": seed,
        "candidates": candidates,
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
    }, 200

@main.route('/design', methods=['POST'])
def design():
    result, status = run_design(request.form)
    if status != 200:
        return jsonify(result), status
    
    # Filter fields for JSON response
    include_rt_product = get_bool(request.form, 'include_rt_product', False)
    include_score      = get_bool(request.form, 'include_score', False)
    include_ascii      = get_bool(request.form, 'include_ascii', False)
    include_struct     = get_bool(request.form, 'include_struct', False)

    filtered = []
    for c in result['candidates']:
        keep = {
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
            "binding_5": f"{c.get('binding_5_start')}-{c.get('binding_5_end')}",
            "binding_3": f"{c.get('binding_3_start')}-{c.get('binding_3_end')}",
            "cross_reactivity": c.get("cross_reactivity", [])
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
    
    result['candidates'] = filtered
    return jsonify(result)

@main.route('/export', methods=['POST'])
def export_csv():
    result, status = run_design(request.form)
    if status != 200:
        return jsonify(result), status
    
    candidates = result['candidates']
    if not candidates:
        return "No candidates to export", 400

    # Create CSV
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Headers
    headers = [
        "Rank", "Score", "RT_Primer", "RT_Tm", "RT_GC", 
        "Forward_Primer", "Fwd_Tm", "Fwd_GC", 
        "Reverse_Primer", "Rev_Tm", "Rev_GC",
        "dG3", "dG5", "dG_Sum"
    ]
    writer.writerow(headers)

    for i, c in enumerate(candidates, 1):
        writer.writerow([
            i,
            f"{c.get('score', 0):.2f}",
            c.get('rt_primer', ''),
            c.get('RT_primer_Tm', ''),
            c.get('RT_primer_GC', ''),
            c.get('Forward_primer', ''),
            c.get('Forward_primer_Tm', ''),
            c.get('Forward_primer_GC', ''),
            c.get('Reverse_primer', ''),
            c.get('Reverse_primer_Tm', ''),
            c.get('Reverse_primer_GC', ''),
            c.get('dg3', ''),
            c.get('dg5', ''),
            c.get('dg_sum', '')
        ])
    
    output.seek(0)
    return Response(
        output.getvalue(),
        mimetype="text/csv",
        headers={"Content-disposition": "attachment; filename=primers.csv"}
    )
