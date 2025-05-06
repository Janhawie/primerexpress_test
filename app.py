from flask import Flask, render_template, request, redirect, url_for, flash, send_file
from werkzeug.utils import secure_filename
import primer3
import io
import csv
import os
import re
from Bio import SeqIO
from tempfile import NamedTemporaryFile
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import time

app = Flask(__name__)
# In production, use environment variable: os.environ.get('SECRET_KEY')
app.secret_key = 'primer_design_secret_key'
app.config['MAX_CONTENT_LENGTH'] = 1 * 1024 * 1024  # 1MB max upload size
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['JSON_AS_ASCII'] = False  # Support for non-ASCII characters in JSON responses

# Create upload folder if it doesn't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Allowed file extensions
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'txt'}

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def parse_fasta(fasta_content):
    """Parse FASTA string content into a dictionary of {id: sequence}"""
    sequences = {}
    
    # Use StringIO to create a file-like object from the string
    fasta_handle = io.StringIO(fasta_content)
    
    # Parse with BioPython
    for record in SeqIO.parse(fasta_handle, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    
    return sequences

def design_primers(sequence, params):
    """Design primers using primer3-py with the provided parameters"""
    try:
        # Prepare primer3 parameters
        seq_args = {
            'SEQUENCE_ID': 'target_sequence',
            'SEQUENCE_TEMPLATE': sequence,
        }
        
        global_args = {
            'PRIMER_OPT_SIZE': int(params.get('primer_opt_size', 20)),
            'PRIMER_MIN_SIZE': int(params.get('primer_min_size', 18)),
            'PRIMER_MAX_SIZE': int(params.get('primer_max_size', 27)),
            'PRIMER_OPT_TM': float(params.get('primer_opt_tm', 60.0)),
            'PRIMER_MIN_TM': float(params.get('primer_min_tm', 57.0)),
            'PRIMER_MAX_TM': float(params.get('primer_max_tm', 63.0)),
            'PRIMER_MIN_GC': float(params.get('primer_min_gc', 20.0)),
            'PRIMER_MAX_GC': float(params.get('primer_max_gc', 80.0)),
            'PRIMER_PRODUCT_SIZE_RANGE': [[int(params.get('product_min_size', 100)), 
                                           int(params.get('product_max_size', 1000))]],
            'PRIMER_NUM_RETURN': int(params.get('num_primer_pairs', 5))
        }
        
        # Add internal oligo (probe) parameters if requested
        if params.get('include_probe', False):
            global_args.update({
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_INTERNAL_OPT_SIZE': int(params.get('probe_opt_size', 20)),
                'PRIMER_INTERNAL_MIN_SIZE': int(params.get('probe_min_size', 18)),
                'PRIMER_INTERNAL_MAX_SIZE': int(params.get('probe_max_size', 27)),
                'PRIMER_INTERNAL_OPT_TM': float(params.get('probe_opt_tm', 60.0)),
                'PRIMER_INTERNAL_MIN_TM': float(params.get('probe_min_tm', 57.0)),
                'PRIMER_INTERNAL_MAX_TM': float(params.get('probe_max_tm', 63.0)),
                'PRIMER_INTERNAL_MIN_GC': float(params.get('probe_min_gc', 20.0)),
                'PRIMER_INTERNAL_MAX_GC': float(params.get('probe_max_gc', 80.0)),
            })
        
        # Using the updated function name (design_primers instead of designPrimers)
        primer_results = primer3.bindings.design_primers(seq_args, global_args)
        return primer_results
    
    except Exception as e:
        print(f"Error in primer design: {e}")
        import traceback
        traceback.print_exc()  # Print full traceback for debugging
        return None

def process_primer_results(primer_results, include_probe=False):
    """Process primer3 results into a more usable format"""
    processed_results = []
    
    if not primer_results:
        print("No primer results to process")
        return processed_results
    
    # Debug the structure of primer_results
    print(f"Keys in primer_results: {list(primer_results.keys())}")
    
    # Check if the format is as expected
    num_primer_pairs = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
    if num_primer_pairs == 0:
        print("No primer pairs returned")
        return processed_results
    
    print(f"Processing {num_primer_pairs} primer pairs")
    
    for i in range(num_primer_pairs):
        try:
            # Check if we have the expected keys
            left_seq_key = f'PRIMER_LEFT_{i}_SEQUENCE'
            right_seq_key = f'PRIMER_RIGHT_{i}_SEQUENCE'
            
            if left_seq_key not in primer_results or right_seq_key not in primer_results:
                print(f"Missing primer sequence keys for pair {i}")
                continue
            
            # Get primer pair data
            pair = {
                'pair_id': i + 1,
                'forward_primer': {
                    'sequence': primer_results.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
                    'position': primer_results.get(f'PRIMER_LEFT_{i}', [0, 0])[0] + 1,  # 1-based position
                    'length': primer_results.get(f'PRIMER_LEFT_{i}', [0, 0])[1],
                    'tm': primer_results.get(f'PRIMER_LEFT_{i}_TM', 0.0),
                    'gc_percent': primer_results.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0.0)
                },
                'reverse_primer': {
                    'sequence': primer_results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                    'position': primer_results.get(f'PRIMER_RIGHT_{i}', [0, 0])[0] + 1,  # 1-based position
                    'length': primer_results.get(f'PRIMER_RIGHT_{i}', [0, 0])[1],
                    'tm': primer_results.get(f'PRIMER_RIGHT_{i}_TM', 0.0),
                    'gc_percent': primer_results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0.0)
                },
                'product_size': primer_results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
                'pair_penalty': primer_results.get(f'PRIMER_PAIR_{i}_PENALTY', 0.0)
            }
            
            # Add probe information if available
            if include_probe and f'PRIMER_INTERNAL_{i}_SEQUENCE' in primer_results:
                pair['probe'] = {
                    'sequence': primer_results.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', ''),
                    'position': primer_results.get(f'PRIMER_INTERNAL_{i}', [0, 0])[0] + 1,  # 1-based position
                    'length': primer_results.get(f'PRIMER_INTERNAL_{i}', [0, 0])[1],
                    'tm': primer_results.get(f'PRIMER_INTERNAL_{i}_TM', 0.0),
                    'gc_percent': primer_results.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT', 0.0)
                }
            
            processed_results.append(pair)
        except Exception as e:
            print(f"Error processing primer pair {i}: {e}")
            continue
    
    print(f"Returning {len(processed_results)} processed primer pairs")
    return processed_results

def generate_csv(results, sequence_id, include_probe=False):
    """Generate CSV content from primer design results"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Write headers
    headers = ['Pair ID', 'Sequence ID',
               'Forward Primer', 'Forward Position', 'Forward Length', 'Forward Tm', 'Forward GC%',
               'Reverse Primer', 'Reverse Position', 'Reverse Length', 'Reverse Tm', 'Reverse GC%',
               'Product Size', 'Pair Penalty']
    
    if include_probe:
        headers.extend(['Probe Sequence', 'Probe Position', 'Probe Length', 'Probe Tm', 'Probe GC%'])
    
    writer.writerow(headers)
    
    # Write data rows
    for pair in results:
        row = [
            pair['pair_id'], sequence_id,
            pair['forward_primer']['sequence'], pair['forward_primer']['position'], pair['forward_primer']['length'],
            f"{pair['forward_primer']['tm']:.2f}", f"{pair['forward_primer']['gc_percent']:.2f}",
            pair['reverse_primer']['sequence'], pair['reverse_primer']['position'], pair['reverse_primer']['length'],
            f"{pair['reverse_primer']['tm']:.2f}", f"{pair['reverse_primer']['gc_percent']:.2f}",
            pair['product_size'], f"{pair['pair_penalty']:.2f}"
        ]
        
        if include_probe and 'probe' in pair:
            row.extend([
                pair['probe']['sequence'], pair['probe']['position'], pair['probe']['length'],
                f"{pair['probe']['tm']:.2f}", f"{pair['probe']['gc_percent']:.2f}"
            ])
        
        writer.writerow(row)
    
    return output.getvalue()

def perform_blast_search(sequence, program="blastn", database="nt", hitlist_size=10):
    """Perform BLAST search using NCBI's web service"""
    try:
        # Perform BLAST search
        result_handle = NCBIWWW.qblast(program, database, sequence, hitlist_size=hitlist_size)
        
        # Parse the results
        from Bio.Blast import NCBIXML
        blast_records = NCBIXML.parse(result_handle)
        
        # Process results
        blast_results = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    blast_results.append({
                        'title': alignment.title,
                        'length': alignment.length,
                        'e_value': hsp.expect,
                        'score': hsp.score,
                        'identities': hsp.identities,
                        'gaps': hsp.gaps,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'sbjct_start': hsp.sbjct_start,
                        'sbjct_end': hsp.sbjct_end,
                        'alignment_length': hsp.align_length
                    })
        
        return blast_results
    except Exception as e:
        print(f"BLAST search error: {str(e)}")
        return None

@app.route('/')
def index():
    # Get a sample sequence for demo purposes
    sample_sequence = """>BRCA1_fragment
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGA
GTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAAT
TTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATAAA
GGAAATACAGCTTGCTCACTTCGTTACTTCAGGCCCCTCAAATCTGTGGGAGCCACACCTGATGACAC
"""
    return render_template('index.html', sample_sequence=sample_sequence)

@app.route('/design', methods=['POST'])
def design():
    try:
        print("=== Starting primer design process ===")
        
        # Get sequence input method
        sequence_input_method = request.form.get('sequence_input_method')
        print(f"Input method: {sequence_input_method}")
        
        # Get sequence data based on input method
        sequences = {}
        
        if sequence_input_method == 'paste':
            fasta_content = request.form.get('sequence_text', '')
            if not fasta_content:
                flash('No sequence provided!', 'error')
                return redirect(url_for('index'))
            
            print(f"Received sequence text of length: {len(fasta_content)}")
            
            # Check if input is in FASTA format
            if not fasta_content.startswith('>'):
                # Try to convert plain sequence to FASTA
                sequence = re.sub(r'[^ATGCN]', '', fasta_content.upper())
                if sequence:
                    fasta_content = f">Sequence\n{sequence}"
                    print("Converted plain sequence to FASTA format")
                else:
                    flash('Invalid sequence format!', 'error')
                    return redirect(url_for('index'))
            
            sequences = parse_fasta(fasta_content)
            
        elif sequence_input_method == 'upload':
            # Check if a file was uploaded
            if 'sequence_file' not in request.files:
                flash('No file uploaded!', 'error')
                return redirect(url_for('index'))
            
            file = request.files['sequence_file']
            
            # Check if the file is valid
            if file.filename == '':
                flash('No file selected!', 'error')
                return redirect(url_for('index'))
            
            if not allowed_file(file.filename):
                flash('Invalid file type! Please upload a FASTA file (.fasta, .fa, or .txt)', 'error')
                return redirect(url_for('index'))
            
            # Read and parse the file
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            
            with open(filepath, 'r') as f:
                fasta_content = f.read()
            
            # Clean up the uploaded file
            os.remove(filepath)
            
            sequences = parse_fasta(fasta_content)
        
        if not sequences:
            flash('No valid sequences found in the input!', 'error')
            return redirect(url_for('index'))
        
        print(f"Parsed {len(sequences)} sequences: {list(sequences.keys())}")
        
        # Get primer design parameters
        params = {
            'primer_opt_size': request.form.get('primer_opt_size', 20),
            'primer_min_size': request.form.get('primer_min_size', 18),
            'primer_max_size': request.form.get('primer_max_size', 27),
            'primer_opt_tm': request.form.get('primer_opt_tm', 60.0),
            'primer_min_tm': request.form.get('primer_min_tm', 57.0),
            'primer_max_tm': request.form.get('primer_max_tm', 63.0),
            'primer_min_gc': request.form.get('primer_min_gc', 20.0),
            'primer_max_gc': request.form.get('primer_max_gc', 80.0),
            'product_min_size': request.form.get('product_min_size', 100),
            'product_max_size': request.form.get('product_max_size', 1000),
            'num_primer_pairs': request.form.get('num_primer_pairs', 5),
            'include_probe': 'include_probe' in request.form,
        }
        
        print(f"Design parameters: {params}")
        
        # Check if parameters are valid
        for key in params:
            if key.startswith(('primer_', 'probe_', 'product_')) and key.endswith(('_size', '_tm', '_gc')):
                try:
                    if key.endswith('_size'):
                        params[key] = int(params[key])
                    else:
                        params[key] = float(params[key])
                except ValueError:
                    flash(f'Invalid value for parameter: {key}', 'error')
                    return redirect(url_for('index'))
        
        # Add probe parameters if needed
        if params['include_probe']:
            params.update({
                'probe_opt_size': request.form.get('probe_opt_size', 20),
                'probe_min_size': request.form.get('probe_min_size', 18),
                'probe_max_size': request.form.get('probe_max_size', 27),
                'probe_opt_tm': request.form.get('probe_opt_tm', 65.0),
                'probe_min_tm': request.form.get('probe_min_tm', 62.0),
                'probe_max_tm': request.form.get('probe_max_tm', 68.0),
                'probe_min_gc': request.form.get('probe_min_gc', 20.0),
                'probe_max_gc': request.form.get('probe_max_gc', 80.0),
            })
            
            # Convert probe parameters to proper types
            for key in params:
                if key.startswith('probe_') and key.endswith(('_size', '_tm', '_gc')):
                    try:
                        if key.endswith('_size'):
                            params[key] = int(params[key])
                        else:
                            params[key] = float(params[key])
                    except ValueError:
                        flash(f'Invalid value for parameter: {key}', 'error')
                        return redirect(url_for('index'))
        
        # Process each sequence
        all_results = {}
        for seq_id, sequence in sequences.items():
            print(f"Designing primers for sequence: {seq_id} (length: {len(sequence)})")
            primer_results = design_primers(sequence, params)
            if primer_results:
                print(f"Successfully designed primers for {seq_id}")
                processed_results = process_primer_results(primer_results, params['include_probe'])
                all_results[seq_id] = {
                    'sequence': sequence,
                    'results': processed_results
                }
                print(f"Found {len(processed_results)} primer pairs for {seq_id}")
            else:
                print(f"Failed to design primers for {seq_id}")
                flash(f'Failed to design primers for sequence: {seq_id}', 'warning')
        
        if not all_results:
            flash('No primers could be designed with the given parameters!', 'error')
            return redirect(url_for('index'))
        
        # Store results in session for download later
        # We'll save CSV content to temp files and keep track of them
        temp_files = {}
        for seq_id, result in all_results.items():
            csv_content = generate_csv(result['results'], seq_id, params['include_probe'])
            temp_file = NamedTemporaryFile(delete=False, suffix='.csv')
            with open(temp_file.name, 'w') as f:
                f.write(csv_content)
            temp_files[seq_id] = temp_file.name
        
        print("Rendering results template")
        return render_template('results.html', 
                               all_results=all_results, 
                               temp_files=temp_files,
                               include_probe=params['include_probe'])
        
    except Exception as e:
        print(f"Critical error in design route: {str(e)}")
        import traceback
        traceback.print_exc()
        flash(f'An error occurred: {str(e)}', 'error')
        return redirect(url_for('index'))

@app.route('/api/design', methods=['POST'])
def api_design():
    """API endpoint for primer design (for AJAX requests)"""
    try:
        print("=== Starting API primer design ===")
        data = request.get_json()
        
        # Extract sequence
        fasta_content = data.get('sequence', '')
        if not fasta_content:
            return {"error": "No sequence provided"}, 400
        
        print(f"Received sequence of length {len(fasta_content)}")
        
        # Check if input is in FASTA format
        if not fasta_content.startswith('>'):
            # Try to convert plain sequence to FASTA
            sequence = re.sub(r'[^ATGCN]', '', fasta_content.upper())
            if sequence:
                fasta_content = f">Sequence\n{sequence}"
                print("Converted plain sequence to FASTA format")
            else:
                return {"error": "Invalid sequence format"}, 400
        
        # Parse the sequence
        sequences = parse_fasta(fasta_content)
        if not sequences:
            return {"error": "No valid sequences found in the input"}, 400
        
        print(f"Parsed {len(sequences)} sequences: {list(sequences.keys())}")
        
        # Get parameters from request data
        params = {
            'primer_opt_size': int(data.get('primer_opt_size', 20)),
            'primer_min_size': int(data.get('primer_min_size', 18)),
            'primer_max_size': int(data.get('primer_max_size', 27)),
            'primer_opt_tm': float(data.get('primer_opt_tm', 60.0)),
            'primer_min_tm': float(data.get('primer_min_tm', 57.0)),
            'primer_max_tm': float(data.get('primer_max_tm', 63.0)),
            'primer_min_gc': float(data.get('primer_min_gc', 20.0)),
            'primer_max_gc': float(data.get('primer_max_gc', 80.0)),
            'product_min_size': int(data.get('product_min_size', 100)),
            'product_max_size': int(data.get('product_max_size', 1000)),
            'num_primer_pairs': int(data.get('num_primer_pairs', 5)),
            'include_probe': data.get('include_probe', False),
        }
        
        print(f"Design parameters: {params}")
        
        # Add probe parameters if needed
        if params['include_probe']:
            params.update({
                'probe_opt_size': int(data.get('probe_opt_size', 20)),
                'probe_min_size': int(data.get('probe_min_size', 18)),
                'probe_max_size': int(data.get('probe_max_size', 27)),
                'probe_opt_tm': float(data.get('probe_opt_tm', 65.0)),
                'probe_min_tm': float(data.get('probe_min_tm', 62.0)),
                'probe_max_tm': float(data.get('probe_max_tm', 68.0)),
                'probe_min_gc': float(data.get('probe_min_gc', 20.0)),
                'probe_max_gc': float(data.get('probe_max_gc', 80.0)),
            })
        
        # Process each sequence
        all_results = {}
        for seq_id, sequence in sequences.items():
            print(f"Designing primers for sequence: {seq_id} (length: {len(sequence)})")
            primer_results = design_primers(sequence, params)
            if primer_results:
                print(f"Successfully designed primers for {seq_id}")
                processed_results = process_primer_results(primer_results, params['include_probe'])
                all_results[seq_id] = {
                    'sequence': sequence[:50] + "..." if len(sequence) > 50 else sequence,  # Truncate for response
                    'sequence_length': len(sequence),
                    'results': processed_results
                }
                print(f"Found {len(processed_results)} primer pairs for {seq_id}")
        
        if not all_results:
            return {"error": "No primers could be designed with the given parameters"}, 400
        
        return {"results": all_results}, 200
        
    except Exception as e:
        print(f"API error: {str(e)}")
        import traceback
        traceback.print_exc()
        return {"error": str(e)}, 500

@app.route('/download/<sequence_id>')
def download(sequence_id):
    try:
        print(f"Attempting to download CSV for sequence: {sequence_id}")
        temp_files_str = request.args.get('temp_files', '')
        
        if not temp_files_str:
            flash('No results to download!', 'error')
            return redirect(url_for('index'))
        
        # Safely evaluate the string representation of the dictionary
        temp_files = eval(temp_files_str)
        
        # Get the temp file path
        temp_file_path = temp_files.get(sequence_id)
        if not temp_file_path or not os.path.exists(temp_file_path):
            flash('Results not found or expired!', 'error')
            return redirect(url_for('index'))
        
        print(f"Sending file: {temp_file_path}")
        
        # Send the file
        return send_file(temp_file_path, 
                        as_attachment=True, 
                        download_name=f'primers_{secure_filename(sequence_id)}.csv',
                        mimetype='text/csv')
    except Exception as e:
        print(f"Error in download: {str(e)}")
        import traceback
        traceback.print_exc()
        flash(f'Error downloading file: {str(e)}', 'error')
        return redirect(url_for('index'))

@app.teardown_appcontext
def cleanup_temp_files(exception=None):
    """Clean up temporary files at the end of the request"""
    # This would need to be improved in a production app to avoid file leaks
    pass

@app.route('/api/blast', methods=['POST'])
def api_blast():
    """API endpoint for BLAST analysis"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        
        if not sequence:
            return {"error": "No sequence provided"}, 400
        
        # Perform BLAST search
        blast_results = perform_blast_search(sequence)
        
        if blast_results is None:
            return {"error": "BLAST search failed"}, 500
        
        return {"results": blast_results}, 200
        
    except Exception as e:
        print(f"API error: {str(e)}")
        return {"error": str(e)}, 500

if __name__ == '__main__':
    app.run(debug=True)