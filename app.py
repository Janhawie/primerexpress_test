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

def calculate_primer_score(primer_pair, include_probe=False):
    """Calculate a comprehensive score for a primer pair with detailed quality checks
    
    Args:
        primer_pair (dict): Dictionary containing primer pair information
        include_probe (bool): Whether to include probe in scoring
        
    Returns:
        dict: Dictionary containing total score and detailed breakdown of scores
    """
    # Initialize scores dictionary for detailed breakdown
    scores = {
        'length_score': 0,
        'gc_score': 0,
        'tm_score': 0,
        'self_comp_score': 0,
        'product_size_score': 0,
        'specificity_score': 0,
        'efficiency_score': 0,
        'total_score': 0
    }
    
    # Score forward primer
    forward = primer_pair['forward_primer']
    forward_seq = forward['sequence']
    
    # Length scoring (20 points)
    length_score = 20 * (1 - abs(forward['length'] - 20) / 10)
    scores['length_score'] += length_score
    
    # GC content scoring (20 points)
    gc_content = forward['gc_percent']
    gc_score = 20 * (1 - abs(gc_content - 50) / 50)
    scores['gc_score'] += gc_score
    
    # Melting temperature scoring (20 points)
    tm = forward['tm']
    tm_score = 20 * (1 - abs(tm - 60) / 10)
    scores['tm_score'] += tm_score
    
    # Score reverse primer
    reverse = primer_pair['reverse_primer']
    reverse_seq = reverse['sequence']
    
    # Add reverse primer scores
    scores['length_score'] += 20 * (1 - abs(reverse['length'] - 20) / 10)
    scores['gc_score'] += 20 * (1 - abs(reverse['gc_percent'] - 50) / 50)
    scores['tm_score'] += 20 * (1 - abs(reverse['tm'] - 60) / 10)
    
    # Self-complementarity and secondary structure checks
    def check_secondary_structures(seq):
        score = 20  # Start with perfect score
        issues = []
        
        # Check for simple repeats
        for i in range(len(seq)-3):
            if seq[i:i+4] in seq[i+4:]:
                score -= 5
                issues.append(f"Repeat found: {seq[i:i+4]}")
        
        # Check for GC-rich regions
        gc_count = seq.count('G') + seq.count('C')
        if gc_count > len(seq) * 0.7:
            score -= 5
            issues.append("GC-rich region detected")
        
        # Check for potential hairpins
        for i in range(len(seq)-4):
            for j in range(i+4, len(seq)-3):
                if seq[i:i+4] == ''.join(reversed(seq[j:j+4])):
                    score -= 5
                    issues.append(f"Potential hairpin: {seq[i:i+4]}")
        
        return max(0, score), issues
    
    # Check both primers for secondary structures
    forward_comp_score, forward_issues = check_secondary_structures(forward_seq)
    reverse_comp_score, reverse_issues = check_secondary_structures(reverse_seq)
    scores['self_comp_score'] = forward_comp_score + reverse_comp_score
    
    # Product size scoring
    product_size = primer_pair['product_size']
    if 100 <= product_size <= 1000:
        product_score = 20 * (1 - abs(product_size - 500) / 500)
    else:
        product_score = 0
    scores['product_size_score'] = product_score
    
    # Specificity scoring based on primer characteristics
    specificity_score = 0
    specificity_issues = []
    
    # Check for balanced base composition
    for seq in [forward_seq, reverse_seq]:
        base_counts = {'A': seq.count('A'), 'T': seq.count('T'), 
                      'G': seq.count('G'), 'C': seq.count('C')}
        if max(base_counts.values()) > len(seq) * 0.4:
            specificity_score -= 5
            specificity_issues.append("Unbalanced base composition")
    
    # Check for 3' end stability
    for seq in [forward_seq, reverse_seq]:
        if seq[-2:] in ['GG', 'CC', 'GC', 'CG']:
            specificity_score += 5
        elif seq[-2:] in ['AA', 'TT', 'AT', 'TA']:
            specificity_score -= 5
            specificity_issues.append("Weak 3' end")
    
    scores['specificity_score'] = specificity_score + 20  # Base score of 20
    
    # Efficiency scoring
    efficiency_score = 0
    efficiency_issues = []
    
    # Check for primer-dimer formation
    def check_primer_dimer(seq1, seq2):
        for i in range(len(seq1)-3):
            for j in range(len(seq2)-3):
                if seq1[i:i+4] == ''.join(reversed(seq2[j:j+4])):
                    return True
        return False
    
    if check_primer_dimer(forward_seq, reverse_seq):
        efficiency_score -= 10
        efficiency_issues.append("Potential primer-dimer formation")
    
    # Check for internal stability
    for seq in [forward_seq, reverse_seq]:
        if seq.count('G') + seq.count('C') < len(seq) * 0.3:
            efficiency_score -= 5
            efficiency_issues.append("Low internal stability")
    
    scores['efficiency_score'] = efficiency_score + 20  # Base score of 20
    
    # Calculate total score
    scores['total_score'] = sum([
        scores['length_score'],
        scores['gc_score'],
        scores['tm_score'],
        scores['self_comp_score'],
        scores['product_size_score'],
        scores['specificity_score'],
        scores['efficiency_score']
    ])
    
    # Add probe scoring if included
    if include_probe and 'probe' in primer_pair:
        probe = primer_pair['probe']
        probe_seq = probe['sequence']
        
        # Add probe scores
        scores['length_score'] += 20 * (1 - abs(probe['length'] - 20) / 10)
        scores['gc_score'] += 10 * (1 - abs(probe['gc_percent'] - 50) / 50)
        scores['tm_score'] += 10 * (1 - abs(probe['tm'] - 65) / 10)
        
        # Check probe for secondary structures
        probe_comp_score, probe_issues = check_secondary_structures(probe_seq)
        scores['self_comp_score'] += probe_comp_score
        
        # Update total score
        scores['total_score'] = sum([
            scores['length_score'],
            scores['gc_score'],
            scores['tm_score'],
            scores['self_comp_score'],
            scores['product_size_score'],
            scores['specificity_score'],
            scores['efficiency_score']
        ])
    
    # Round all scores to 2 decimal places
    for key in scores:
        scores[key] = round(scores[key], 2)
    
    # Add detailed issues to the scores dictionary
    scores['issues'] = {
        'forward_primer': forward_issues,
        'reverse_primer': reverse_issues,
        'specificity': specificity_issues,
        'efficiency': efficiency_issues
    }
    
    return scores

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
            
            # Calculate and add the primer scores
            pair['scores'] = calculate_primer_score(pair, include_probe)
            
            processed_results.append(pair)
        except Exception as e:
            print(f"Error processing primer pair {i}: {e}")
            continue
    
    # Sort results by pair_id in ascending order
    processed_results.sort(key=lambda x: x['pair_id'])
    
    print(f"Returning {len(processed_results)} processed primer pairs")
    return processed_results

def generate_csv(results, sequence_id, include_probe=False):
    """Generate CSV content from primer design results with detailed scoring"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Write headers
    headers = ['Pair ID', 'Sequence ID', 'Total Score',
               'Length Score', 'GC Score', 'Tm Score', 'Self-Comp Score',
               'Product Size Score', 'Specificity Score', 'Efficiency Score',
               'Forward Primer', 'Forward Position', 'Forward Length', 'Forward Tm', 'Forward GC%',
               'Reverse Primer', 'Reverse Position', 'Reverse Length', 'Reverse Tm', 'Reverse GC%',
               'Product Size', 'Pair Penalty', 'Issues']
    
    if include_probe:
        headers.extend(['Probe Sequence', 'Probe Position', 'Probe Length', 'Probe Tm', 'Probe GC%'])
    
    writer.writerow(headers)
    
    # Write data rows
    for pair in results:
        # Format issues as a string
        issues = []
        for category, category_issues in pair['scores']['issues'].items():
            if category_issues:
                issues.append(f"{category}: {', '.join(category_issues)}")
        issues_str = '; '.join(issues) if issues else 'No issues'
        
        row = [
            pair['pair_id'], sequence_id, f"{pair['scores']['total_score']:.2f}",
            f"{pair['scores']['length_score']:.2f}",
            f"{pair['scores']['gc_score']:.2f}",
            f"{pair['scores']['tm_score']:.2f}",
            f"{pair['scores']['self_comp_score']:.2f}",
            f"{pair['scores']['product_size_score']:.2f}",
            f"{pair['scores']['specificity_score']:.2f}",
            f"{pair['scores']['efficiency_score']:.2f}",
            pair['forward_primer']['sequence'], pair['forward_primer']['position'], pair['forward_primer']['length'],
            f"{pair['forward_primer']['tm']:.2f}", f"{pair['forward_primer']['gc_percent']:.2f}",
            pair['reverse_primer']['sequence'], pair['reverse_primer']['position'], pair['reverse_primer']['length'],
            f"{pair['reverse_primer']['tm']:.2f}", f"{pair['reverse_primer']['gc_percent']:.2f}",
            pair['product_size'], f"{pair['pair_penalty']:.2f}",
            issues_str
        ]
        
        if include_probe and 'probe' in pair:
            row.extend([
                pair['probe']['sequence'], pair['probe']['position'], pair['probe']['length'],
                f"{pair['probe']['tm']:.2f}", f"{pair['probe']['gc_percent']:.2f}"
            ])
        
        writer.writerow(row)
    
    return output.getvalue()

def perform_blast_search(sequence, program="blastn", database="nt", hitlist_size=15, expect_value=1e-3):
    try:
        # Check sequence quality first
        quality_check = check_sequence_quality(sequence)
        
        if quality_check['quality_issues']:
            return {
                "error": "Sequence quality issues found",
                "quality_report": quality_check
            }
        
        # Use the cleaned sequence for BLAST
        cleaned_sequence = quality_check['cleaned_sequence']
        
        # Clean and validate sequence while preserving header
        cleaned_sequence = cleaned_sequence.strip()
        
        # Split into header and sequence
        if cleaned_sequence.startswith('>'):
            header, seq = cleaned_sequence.split('\n', 1)
            seq = seq.replace('\n', '').strip().upper()
        else:
            header = ">Sequence"
            seq = cleaned_sequence.replace('\n', '').strip().upper()
        
        # Detailed sequence validation
        if not seq:
            return {"error": "Empty sequence provided"}
            
        if not all(c in 'ATCGN' for c in seq):
            invalid_chars = set(c for c in seq if c not in 'ATCGN')
            return {"error": f"Invalid sequence characters found: {invalid_chars}. Only A, T, C, G, and N are allowed."}
        
        if len(seq) < 20:
            return {"error": f"Sequence is too short ({len(seq)} nucleotides). Minimum length should be 20 nucleotides."}
            
        if len(seq) > 10000:
            return {"error": f"Sequence is too long ({len(seq)} nucleotides). Maximum length is 10000 nucleotides."}
            
        # Check for low complexity regions
        if seq.count('A')/len(seq) > 0.8 or seq.count('T')/len(seq) > 0.8 or \
           seq.count('G')/len(seq) > 0.8 or seq.count('C')/len(seq) > 0.8:
            return {"error": "Sequence appears to be low complexity (too many repeated bases)"}
        
        print(f"Performing BLAST search with sequence length: {len(seq)} and E-value threshold: {expect_value}")
        
        # Perform BLAST search with stricter parameters
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=seq,
            hitlist_size=hitlist_size,
            expect=expect_value,  # Stricter E-value
            word_size=7,  # Decrease from 11 to 7 for more matches
            gapcosts="5 2",
            filter="L",
            perc_ident=90  # Only show matches with 90% or higher identity
        )
        
        # Parse the results
        from Bio.Blast import NCBIXML
        blast_records = NCBIXML.parse(result_handle)
        
        # Process results with stricter filtering
        blast_results = []
        total_matches = 0
        filtered_matches = 0
        
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    total_matches += 1
                    # Calculate percent identity
                    percent_identity = (hsp.identities / hsp.align_length) * 100
                    
                    # Only include highly significant matches
                    if percent_identity >= 90 and hsp.expect <= 1e-5:
                        filtered_matches += 1
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
                            'alignment_length': hsp.align_length,
                            'percent_identity': percent_identity,
                            'query_header': header,
                            'is_reverse_strand': hsp.sbjct_start > hsp.sbjct_end,
                            'alignment_direction': 'Reverse' if hsp.sbjct_start > hsp.sbjct_end else 'Forward',
                            'match_quality': 'High' if hsp.expect < 1e-10 else 'Medium',
                            'significance': 'Very significant' if hsp.expect < 1e-10 else 'Significant'
                        })
        
        # Sort results by normalized score in descending order and take top 5
        blast_results.sort(key=lambda x: x['normalized_score'], reverse=True)
        top_results = blast_results[:5]
        
        # Add relevance ranking
        for i, result in enumerate(top_results, 1):
            result['relevance_rank'] = i
            result['relevance_description'] = {
                1: "Most relevant match",
                2: "Second most relevant match",
                3: "Third most relevant match",
                4: "Fourth most relevant match",
                5: "Fifth most relevant match"
            }[i]
        
        if not top_results:
            return {
                "message": "Excellent! No significant matches found in the database.",
                "details": {
                    "status": "success",
                    "message_type": "positive",
                    "implications": [
                        "Your primer is highly specific",
                        "Low risk of off-target amplification",
                        "Reduced chance of cross-reactivity",
                        "Primer appears to be unique to your target"
                    ],
                    "style": {
                        "color": "green",
                        "background": "#e8f5e9",
                        "border": "1px solid #81c784",
                        "padding": "15px",
                        "border-radius": "5px"
                    }
                }
            }, 200
        
        return top_results
    except Exception as e:
        print(f"BLAST search error: {str(e)}")
        return {"error": f"BLAST search failed: {str(e)}"}

def check_primer_specificity(sequence, program="blastn", database="nt", hitlist_size=15):
    """Check PCR primer specificity using BLAST
    
    Args:
        sequence (str): The primer sequence to check
        program (str): BLAST program to use (default: blastn)
        database (str): Database to search against (default: nt)
        hitlist_size (int): Number of results to return (default: 15)
    """
    try:
        # Clean and validate sequence while preserving header
        sequence = sequence.strip()
        
        # Split into header and sequence
        if sequence.startswith('>'):
            header, seq = sequence.split('\n', 1)
            seq = seq.replace('\n', '').strip().upper()
        else:
            header = ">Sequence"
            seq = sequence.replace('\n', '').strip().upper()
        
        # Detailed sequence validation
        if not seq:
            return {"error": "Empty sequence provided"}
            
        if not all(c in 'ATCGN' for c in seq):
            invalid_chars = set(c for c in seq if c not in 'ATCGN')
            return {"error": f"Invalid sequence characters found: {invalid_chars}. Only A, T, C, G, and N are allowed."}
        
        if len(seq) < 18:
            return {"error": f"Primer is too short ({len(seq)} nucleotides). Minimum length should be 18 nucleotides."}
            
        if len(seq) > 30:
            return {"error": f"Primer is too long ({len(seq)} nucleotides). Maximum length should be 30 nucleotides."}
        
        print(f"Checking primer specificity with sequence length: {len(seq)}")
        
        # Perform BLAST search with strict parameters for primer specificity
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=seq,
            hitlist_size=hitlist_size,
            expect=0.001,  # Strict E-value for primer specificity
            word_size=7,  # Smaller word size for more sensitive search
            gapcosts="5 2",  # Standard gap costs
            filter="L"  # Low complexity filter
        )
        
        # Parse the results
        from Bio.Blast import NCBIXML
        blast_records = NCBIXML.parse(result_handle)
        
        # Process results
        blast_results = []
        total_matches = 0
        filtered_matches = 0
        
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    total_matches += 1
                    # Calculate percent identity
                    percent_identity = (hsp.identities / hsp.align_length) * 100
                    
                    # Strict filtering for primer specificity
                    if percent_identity >= 90 and hsp.expect <= 0.001:
                        filtered_matches += 1
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
                            'alignment_length': hsp.align_length,
                            'percent_identity': percent_identity,
                            'query_header': header,
                            'specificity_warning': 'High specificity risk' if hsp.expect > 0.0001 else 'Acceptable specificity'
                        })
        
        if not blast_results:
            return {
                "message": "No significant matches found. Primer appears to be specific.",
                "sequence_info": {
                    "header": header,
                    "length": len(seq),
                    "gc_content": (seq.count('G') + seq.count('C')) / len(seq) * 100,
                    "n_content": seq.count('N') / len(seq) * 100,
                    "total_matches_found": total_matches,
                    "matches_after_filtering": filtered_matches
                }
            }
        
        return {
            "results": blast_results,
            "specificity_assessment": "Primer may not be specific enough for PCR. Consider redesigning." if filtered_matches > 1 else "Primer appears to be specific."
        }
    except Exception as e:
        print(f"BLAST search error: {str(e)}")
        return {"error": f"BLAST search failed: {str(e)}"}

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
            'primer_opt_size': 20,
            'primer_min_size': 18,
            'primer_max_size': 27,
            'primer_opt_tm': 60.0,
            'primer_min_tm': 57.0,
            'primer_max_tm': 63.0,
            'primer_min_gc': 20.0,
            'primer_max_gc': 80.0,
            'product_min_size': 100,
            'product_max_size': 1000,
            'num_primer_pairs': 5,
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
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        
        if not sequence:
            return {"error": "No sequence provided"}, 400
        
        # Check sequence quality
        quality_check = check_sequence_quality(sequence)
        
        if quality_check['quality_issues']:
            return {
                "message": "Sequence quality issues found",
                "quality_report": {
                    "issues": quality_check['quality_issues'],
                    "suggestions": quality_check['improvements'],
                    "base_composition": quality_check['base_composition'],
                    "sequence_length": quality_check['sequence_length']
                }
            }, 200
        
        # Perform BLAST search with cleaned sequence
        blast_results = perform_blast_search(quality_check['cleaned_sequence'])
        
        if isinstance(blast_results, dict) and "error" in blast_results:
            return {"error": blast_results["error"]}, 400
        
        if not blast_results:
            return {
                "message": "Excellent! No significant matches found in the database.",
                "details": {
                    "status": "success",
                    "message_type": "positive",
                    "implications": [
                        "Your primer is highly specific",
                        "Low risk of off-target amplification",
                        "Reduced chance of cross-reactivity",
                        "Primer appears to be unique to your target"
                    ],
                    "style": {
                        "color": "green",
                        "background": "#e8f5e9",
                        "border": "1px solid #81c784",
                        "padding": "15px",
                        "border-radius": "5px"
                    }
                }
            }, 200
            
        return {"results": blast_results}, 200
        
    except Exception as e:
        return {"error": str(e)}, 500

@app.route('/api/check_primer', methods=['POST'])
def check_primer():
    """API endpoint for primer specificity checking"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        
        if not sequence:
            return {"error": "No sequence provided"}, 400
        
        # Check primer specificity
        results = check_primer_specificity(sequence)
        
        # Check if the result is an error message
        if isinstance(results, dict) and "error" in results:
            return {"error": results["error"]}, 400
        
        return results, 200
        
    except Exception as e:
        print(f"API error: {str(e)}")
        return {"error": str(e)}, 500

@app.route('/api/improve_sequence', methods=['POST'])
def improve_sequence():
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        
        if not sequence:
            return {"error": "No sequence provided"}, 400
        
        quality_check = check_sequence_quality(sequence)
        
        return {
            "original_sequence": sequence,
            "cleaned_sequence": quality_check['cleaned_sequence'],
            "quality_report": {
                "issues": quality_check['quality_issues'],
                "suggestions": quality_check['improvements'],
                "base_composition": quality_check['base_composition'],
                "sequence_length": quality_check['sequence_length']
            }
        }, 200
        
    except Exception as e:
        return {"error": str(e)}, 500

def test_primer_scoring():
    """Test function to demonstrate primer scoring with a real example"""
    # Example primer pair
    test_pair = {
        'forward_primer': {
            'sequence': 'ATCGATCGATCGATCGATCG',  # 20bp, 50% GC
            'position': 1,
            'length': 20,
            'tm': 60.0,
            'gc_percent': 50.0
        },
        'reverse_primer': {
            'sequence': 'GCTAGCTAGCTAGCTAGCTA',  # 20bp, 50% GC
            'position': 100,
            'length': 20,
            'tm': 60.0,
            'gc_percent': 50.0
        },
        'product_size': 500
    }
    
    # Calculate scores
    scores = calculate_primer_score(test_pair)
    
    # Print detailed scoring breakdown
    print("\n=== Example Primer Pair Scoring ===")
    print("\nForward Primer:", test_pair['forward_primer']['sequence'])
    print("Reverse Primer:", test_pair['reverse_primer']['sequence'])
    print("\nScoring Breakdown:")
    print(f"1. Length Score: {scores['length_score']:.2f}/20")
    print(f"   - Forward: 20bp (optimal)")
    print(f"   - Reverse: 20bp (optimal)")
    
    print(f"\n2. GC Content Score: {scores['gc_score']:.2f}/20")
    print(f"   - Forward: 50% GC (optimal)")
    print(f"   - Reverse: 50% GC (optimal)")
    
    print(f"\n3. Melting Temperature Score: {scores['tm_score']:.2f}/20")
    print(f"   - Forward: 60°C (optimal)")
    print(f"   - Reverse: 60°C (optimal)")
    
    print(f"\n4. Self-Complementarity Score: {scores['self_comp_score']:.2f}/20")
    print("   - Checks for:")
    print("     * Repeated sequences")
    print("     * GC-rich regions")
    print("     * Potential hairpins")
    
    print(f"\n5. Product Size Score: {scores['product_size_score']:.2f}/20")
    print(f"   - Product size: {test_pair['product_size']}bp (optimal)")
    
    print(f"\n6. Specificity Score: {scores['specificity_score']:.2f}/20")
    print("   - Checks for:")
    print("     * Balanced base composition")
    print("     * Strong 3' end")
    
    print(f"\n7. Efficiency Score: {scores['efficiency_score']:.2f}/20")
    print("   - Checks for:")
    print("     * Primer-dimer formation")
    print("     * Internal stability")
    
    print(f"\nTotal Score: {scores['total_score']:.2f}/140")
    
    if scores['issues']:
        print("\nIssues Found:")
        for category, issues in scores['issues'].items():
            if issues:
                print(f"- {category}: {', '.join(issues)}")
    else:
        print("\nNo issues found - This is an excellent primer pair!")
    
    return scores

def test_problematic_primer_scoring():
    """Test function to demonstrate scoring of a problematic primer pair"""
    # Example of a problematic primer pair
    problematic_pair = {
        'forward_primer': {
            'sequence': 'ATATATATATATATATATAT',  # 20bp, 0% GC, repetitive
            'position': 1,
            'length': 20,
            'tm': 45.0,  # Too low
            'gc_percent': 0.0  # Too low
        },
        'reverse_primer': {
            'sequence': 'GCGCGCGCGCGCGCGCGCGC',  # 20bp, 100% GC, repetitive
            'position': 100,
            'length': 20,
            'tm': 75.0,  # Too high
            'gc_percent': 100.0  # Too high
        },
        'product_size': 1500  # Too large
    }
    
    # Calculate scores
    scores = calculate_primer_score(problematic_pair)
    
    # Print detailed scoring breakdown
    print("\n=== Problematic Primer Pair Scoring ===")
    print("\nForward Primer:", problematic_pair['forward_primer']['sequence'])
    print("Reverse Primer:", problematic_pair['reverse_primer']['sequence'])
    print("\nScoring Breakdown:")
    print(f"1. Length Score: {scores['length_score']:.2f}/20")
    print(f"   - Forward: 20bp (optimal length but poor sequence)")
    print(f"   - Reverse: 20bp (optimal length but poor sequence)")
    
    print(f"\n2. GC Content Score: {scores['gc_score']:.2f}/20")
    print(f"   - Forward: 0% GC (too low)")
    print(f"   - Reverse: 100% GC (too high)")
    
    print(f"\n3. Melting Temperature Score: {scores['tm_score']:.2f}/20")
    print(f"   - Forward: 45°C (too low)")
    print(f"   - Reverse: 75°C (too high)")
    
    print(f"\n4. Self-Complementarity Score: {scores['self_comp_score']:.2f}/20")
    print("   - Forward: Highly repetitive (AT repeats)")
    print("   - Reverse: Highly repetitive (GC repeats)")
    
    print(f"\n5. Product Size Score: {scores['product_size_score']:.2f}/20")
    print(f"   - Product size: {problematic_pair['product_size']}bp (too large)")
    
    print(f"\n6. Specificity Score: {scores['specificity_score']:.2f}/20")
    print("   - Forward: Weak 3' end (AT)")
    print("   - Reverse: Strong 3' end but unbalanced composition")
    
    print(f"\n7. Efficiency Score: {scores['efficiency_score']:.2f}/20")
    print("   - Forward: Low internal stability")
    print("   - Reverse: High internal stability but potential for secondary structures")
    
    print(f"\nTotal Score: {scores['total_score']:.2f}/140")
    
    print("\nIssues Found:")
    for category, issues in scores['issues'].items():
        if issues:
            print(f"- {category}: {', '.join(issues)}")
    
    return scores

def is_biologically_significant(hsp, alignment):
    # Check E-value
    if hsp.expect > 1e-5:
        return False
        
    # Check percent identity
    percent_identity = (hsp.identities / hsp.align_length) * 100
    if percent_identity < 90:
        return False
        
    # Check alignment length
    if hsp.align_length < 20:  # Minimum alignment length
        return False
        
    # Check for low complexity regions
    if 'low complexity' in alignment.title.lower():
        return False
        
    return True

def check_sequence_quality(sequence):
    """Check and improve sequence quality"""
    quality_issues = []
    improvements = []
    
    # Convert to uppercase
    sequence = sequence.upper()
    
    # Remove whitespace and newlines
    sequence = ''.join(sequence.split())
    
    # Check for invalid characters
    valid_bases = set('ATCGN')
    invalid_chars = set(c for c in sequence if c not in valid_bases)
    if invalid_chars:
        quality_issues.append(f"Invalid characters found: {invalid_chars}")
        # Remove invalid characters
        sequence = ''.join(c for c in sequence if c in valid_bases)
        improvements.append("Removed invalid characters")
    
    # Check for ambiguous bases (N)
    n_count = sequence.count('N')
    if n_count > 0:
        quality_issues.append(f"Found {n_count} ambiguous bases (N)")
        improvements.append("Consider replacing N with specific bases if known")
    
    # Check for low complexity regions
    for base in 'ATCG':
        if sequence.count(base)/len(sequence) > 0.8:
            quality_issues.append(f"High frequency of {base} bases")
            improvements.append("Sequence may have low complexity regions")
    
    # Check for minimum length
    if len(sequence) < 20:
        quality_issues.append("Sequence is too short")
        improvements.append("Sequence should be at least 20 bases long")
    
    # Check for maximum length
    if len(sequence) > 10000:
        quality_issues.append("Sequence is too long")
        improvements.append("Sequence should be less than 10000 bases")
    
    # Check for balanced base composition
    base_counts = {base: sequence.count(base) for base in 'ATCG'}
    total_bases = sum(base_counts.values())
    for base, count in base_counts.items():
        if count/total_bases < 0.1:
            quality_issues.append(f"Low frequency of {base} bases")
            improvements.append(f"Consider adding more {base} bases if possible")
    
    return {
        'cleaned_sequence': sequence,
        'quality_issues': quality_issues,
        'improvements': improvements,
        'base_composition': {base: f"{(count/total_bases)*100:.1f}%" for base, count in base_counts.items()},
        'sequence_length': len(sequence)
    }

if __name__ == '__main__':
    # Run both tests before starting the app
    print("\n=== Testing Good Primer Pair ===")
    test_primer_scoring()
    print("\n=== Testing Problematic Primer Pair ===")
    test_problematic_primer_scoring()
    app.run(debug=True)