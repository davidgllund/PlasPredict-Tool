#!/usr/bin/env python3
"""
Flask web application for Plasmid Host Range Prediction
"""
import os
import sys
import tempfile
import logging
import json
from pathlib import Path
from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename
import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter
from itertools import product
import glob
import subprocess
import re
from pathlib import Path

from config import get_config

app = Flask(__name__)

env = os.environ.get('FLASK_ENV', 'development')
config = get_config(env)
app.config.from_object(config)

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

logger.info(f"Loading model from {app.config['MODEL_PATH']}")
if not Path(app.config['MODEL_PATH']).exists():
    raise FileNotFoundError(f"Model not found at {app.config['MODEL_PATH']}")
model = joblib.load(str(app.config['MODEL_PATH']))
logger.info("Model loaded successfully")

def generate_possible_kmers(k):
    """Generate all possible k-mers of length k"""
    return [''.join(p) for p in product('ACGT', repeat=k)]


def get_kmers(seq, k):
    """Extract all k-mers from a sequence"""
    return [str(seq[i:i+k]) for i in range(len(seq) - k + 1)]


def get_kmer_distribution(kmers, possible_kmers):
    """Calculate k-mer frequency distribution"""
    if len(kmers) == 0:
        return np.zeros(len(possible_kmers))
    
    kmers_counted = Counter(kmers)
    kmer_frequencies = np.array([kmers_counted[kmer] / len(kmers) for kmer in possible_kmers])
    return kmer_frequencies


def get_canonical_kmers(df):
    comp = str.maketrans("ACGT", "TGCA")

    def revcomp(seq):
        return seq.translate(comp)[::-1]

    canonical_map = {}
    for kmer in df.columns:
        rc = revcomp(kmer)
        canonical = min(kmer, rc)
        canonical_map.setdefault(canonical, []).append(kmer)

    canonical_df = pd.DataFrame(index=df.index)
    for canon, kmers in canonical_map.items():
        canonical_df[canon] = df[kmers].sum(axis=1)

    return canonical_df


def calc_kmer_distributions(genome, k, possible_kmers):
    """Calculate k-mer distribution for entire genome"""
    genome_kmers = []
    for seq_entry in genome.values():
        genome_kmers.extend(get_kmers(seq_entry.seq, k))
    
    distribution = get_kmer_distribution(genome_kmers, possible_kmers)
    return distribution


def annotate_conjugation_system(filepath):
    """Identify conjugation systems in plasmid using MacSyFinder CONJScan"""
    conjugation_system = [] 
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        proteins = tmp / "proteins.faa"
        
        try:
            logger.info(f"Running prodigal on {filepath}")
            result = subprocess.run(
                ["prodigal", "-i", filepath, "-a", str(proteins), "-p", "meta"],
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                logger.error(f"Prodigal failed with return code {result.returncode}")
                logger.error(f"Prodigal stderr: {result.stderr}")
                logger.error(f"Prodigal stdout: {result.stdout}")
                return conjugation_system
            
            logger.info(f"Prodigal succeeded, proteins written to {proteins}")
            logger.info(f"Proteins file exists: {proteins.exists()}")

            outdir = tmp / "conjscan_results"
            tsv_path = outdir / "best_solution.tsv"

            cmd = [
                "macsyfinder",
                "--db-type", "ordered_replicon",
                "--sequence-db", str(proteins),
                "--models", "CONJScan/Plasmids",
                "all",
                "--out-dir", str(outdir)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error("MacSyFinder failed")
                logger.error(result.stderr)
                logger.error(result.stdout)
                return []
            
            if not tsv_path.exists():
                logger.error(f"MacSyFinder output missing: {tsv_path}")
                return []
            
            conjscan_df = pd.read_csv(tsv_path, sep="\t", comment="#")
            conjugation_system = list(set([x.split('_')[1] for x in conjscan_df['gene_name']]))
        
        except Exception as e:
            logger.error(f"Exception in conjugation annotation: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    logger.info(f"Conjugation annotation complete. Found: {conjugation_system}")
    return conjugation_system


def run_plasmidfinder(fasta_path, output_dir, db_path):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, "-m", "plasmidfinder",
        "-i", str(fasta_path),
        "-o", str(output_dir),
        "-p", str(db_path),
        "-x"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.warning(f"Plasmidfinder failed (returncode {result.returncode})")
        logger.warning(f"Plasmidfinder stderr: {result.stderr}")
        return None
 
    stdout = result.stdout
    json_start = stdout.find('{')
    if json_start == -1:
        logger.warning("No JSON found in plasmidfinder stdout")
        return None

    try:
        data = json.loads(stdout[json_start:])
        json_path = output_dir / "data.json"
        json_path.write_text(json.dumps(data))
        logger.info("Plasmidfinder completed successfully")
        return output_dir
    except json.JSONDecodeError as e:
        logger.warning(f"Could not parse plasmidfinder JSON output: {e}")
        return None


def parse_plasmidfinder_json(json_path):
    """Parse plasmidfinder v3 JSON output to extract incompatibility types"""
    try:
        with open(json_path) as f:
            data = json.load(f)

        rows = []
        for key, region in data.get("seq_regions", {}).items():
            if not isinstance(region, dict):
                continue
            rows.append({
                "plasmid":    region.get("name"),
                "identity":   region.get("identity"),
                "coverage":   region.get("coverage"),
                "contig":     region.get("query_id"),
                "reference_accession": region.get("ref_acc"),
                "note":       region.get("note"),
            })

        return pd.DataFrame(rows)

    except Exception as e:
        logger.warning(f"Could not parse plasmidfinder output: {e}")
        return pd.DataFrame()
    


def run_rgi(fasta_path, output_prefix, executable="rgi"):
    """Run RGI to identify antibiotic resistance genes"""
    fasta_path = Path(fasta_path).resolve()
    output_prefix = Path(output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running RGI on {fasta_path}")
    logger.info(f"RGI output prefix: {output_prefix}")

    cmd = [
        executable,
        "main",
        "--input_sequence", str(fasta_path),
        "--output_file", str(output_prefix),
        "--input_type", "contig",
        "--local"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    
    logger.info(f"RGI stdout: {result.stdout}")
    logger.info(f"RGI stderr: {result.stderr}")
    logger.info(f"RGI return code: {result.returncode}")

    if result.returncode != 0:
        logger.error(f"RGI failed with return code {result.returncode}")
        return None

    logger.info(f"RGI completed successfully")
    return output_prefix


def parse_rgi_output(file_path: str) -> pd.DataFrame:
    """Parse and filter RGI output to remove overlapping hits"""
    try:
        df = pd.read_csv(file_path, sep="\t")

        if df.empty:
            return df

        df["start_pos"] = df[["Start", "Stop"]].min(axis=1)
        df["end_pos"] = df[["Start", "Stop"]].max(axis=1)

        tier_mapping = {"Perfect": 3, "Strict": 2, "Loose": 1}
        df["tier_score"] = df["Cut_Off"].map(tier_mapping).fillna(0)

        df = df.sort_values(
            by=["Contig", "start_pos", "tier_score", "Best_Hit_Bitscore", "Best_Identities"],
            ascending=[True, True, False, False, False],
        ).reset_index(drop=True)

        keep_indices = []

        for contig, group in df.groupby("Contig"):
            current_end = -1
            current_tier = -1
            current_bitscore = -1
            current_idx = None

            for idx, row in group.iterrows():
                if row["start_pos"] <= current_end:
                    is_better_tier = row["tier_score"] > current_tier
                    is_equal_tier_better_score = (
                        row["tier_score"] == current_tier
                    ) and (row["Best_Hit_Bitscore"] > current_bitscore)

                    if is_better_tier or is_equal_tier_better_score:
                        if current_idx in keep_indices:
                            keep_indices.remove(current_idx)

                        current_idx = idx
                        current_end = max(current_end, row["end_pos"])
                        current_tier = row["tier_score"]
                        current_bitscore = row["Best_Hit_Bitscore"]
                        keep_indices.append(idx)
                    else:
                        continue
                else:
                    current_idx = idx
                    current_end = row["end_pos"]
                    current_tier = row["tier_score"]
                    current_bitscore = row["Best_Hit_Bitscore"]
                    keep_indices.append(idx)
 
        filtered_df = df.loc[keep_indices].copy()
 
        filtered_df = filtered_df.drop(
            columns=["start_pos", "end_pos", "tier_score"]
        )

        return filtered_df.sort_values(by=["Contig", "Start"]).reset_index(
            drop=True
        )
    except Exception as e:
        logger.warning(f"Could not parse RGI output: {e}")
        return pd.DataFrame()


def predict_host_range(fasta_content, isolation_sources):
    """
    Make host range prediction for a plasmid sequence
    
    Args:
        fasta_content: FASTA sequence as string or file path
        isolation_sources: List of isolation sources selected (optional, 0 or 1)
    
    Returns:
        Dictionary with prediction results
    """
    try:
        if os.path.isfile(fasta_content):
            plasmid_sequence = SeqIO.to_dict(
                SeqIO.parse(fasta_content, 'fasta')
            )
        else:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as f:
                f.write(fasta_content)
                temp_fasta = f.name
            plasmid_sequence = SeqIO.to_dict(
                SeqIO.parse(temp_fasta, 'fasta')
            )
            os.unlink(temp_fasta)
        
        if not plasmid_sequence:
            return {'error': 'No valid sequences found in FASTA file', 'success': False}
        
        df = pd.DataFrame(
            0,
            index=plasmid_sequence.keys(),
            columns=model.get_booster().feature_names
        )
        
        all_kmers = generate_possible_kmers(3)
        kmer_distributions = calc_kmer_distributions(plasmid_sequence, 3, all_kmers)
        kmer_distributions = get_canonical_kmers(kmer_distributions)
        kmer_df = pd.DataFrame(
            [kmer_distributions],
            index=plasmid_sequence.keys(),
            columns=all_kmers
        )
        
        kmer_cols = list(set(kmer_df.columns) & set(df.columns))
        df[kmer_cols] = kmer_df[kmer_cols]
        
        first_seq_id = list(plasmid_sequence.keys())[0]
        df.loc[first_seq_id, 'size'] = len(plasmid_sequence[first_seq_id].seq)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as f:
            for seq_id, seq_record in plasmid_sequence.items():
                f.write(f">{seq_id}\n{str(seq_record.seq)}\n")
            temp_fasta = f.name
        
        conjugation_system = annotate_conjugation_system(temp_fasta)
        conjugation_system_detected = [c if len(c) > 1 else f"MPF{c}" for c in conjugation_system]
        conjugation_system = list(set(conjugation_system) & set(df.columns))
        logger.info(f"Detected conjugation systems: {conjugation_system}")
        
        inc_types_detected = []
        inc_types = []
        try:
            with tempfile.TemporaryDirectory() as plsmd_tmpdir:
                plsmd_output = run_plasmidfinder(
                    temp_fasta,
                    plsmd_tmpdir,
                    app.config['PLASMIDFINDER_DB_PATH']
                )
                if plsmd_output:
                    json_path = Path(plsmd_output) / "data.json"
                    if json_path.exists():
                        plasmidfinder_output = parse_plasmidfinder_json(str(json_path))
                        if not plasmidfinder_output.empty:
                            inc_types_detected = list(set(plasmidfinder_output['plasmid']))
                            inc_type_cleaned_temp = [re.sub(r"\([^)]*\)", "", item) for item in inc_types_detected]
                            inc_types = list(set(inc_type_cleaned_temp) & set(df.columns))
                            logger.info(f"Detected Inc types (raw): {inc_types_detected}")
                            logger.info(f"Inc types for model: {inc_types}")
        except Exception as e:
            logger.warning(f"Could not run plasmidfinder: {e}")
        
        drug_classes = []
        resistance_mechanisms = []
        args_detected = []
        try:
            with tempfile.TemporaryDirectory() as rgi_tmpdir:
                rgi_prefix = Path(rgi_tmpdir) / "rgi_out"
                logger.info(f"RGI temporary directory: {rgi_tmpdir}")
                rgi_output_prefix = run_rgi(temp_fasta, str(rgi_prefix))
                if rgi_output_prefix:
                    rgi_txt = str(rgi_output_prefix) + ".txt"
                    logger.info(f"Looking for RGI output at: {rgi_txt}")
                    logger.info(f"RGI output file exists: {os.path.exists(rgi_txt)}")
                    
                    if os.path.exists(rgi_txt):
                        rgi_output = parse_rgi_output(rgi_txt)
                        logger.info(f"RGI output shape: {rgi_output.shape}")
                        
                        if not rgi_output.empty:
                            arg_list = [
                                arg.strip()
                                for row in rgi_output["Best_Hit_ARO"].dropna().unique()
                                for arg in row.split(";")
                            ] 
                            args_detected = list(set(arg_list))
                            logger.info(f"Detected antibiotic resistance genes: {args_detected}")

                            drug_class_list = [
                                drug.strip()
                                for row in rgi_output["Drug Class"].dropna().unique()
                                for drug in row.split(";")
                            ] 

                            drug_classes = list(set(drug_class_list) & set(df.columns))
                            logger.info(f"Drug classes for model: {drug_classes}")
                            
                            resistance_mech_list = [
                                mechanism.strip()
                                for row in rgi_output["Resistance Mechanism"].dropna().unique()
                                for mechanism in row.split(";")
                            ] 

                            resistance_mechanisms = list(set(resistance_mech_list) & set(df.columns))
                            logger.info(f"Resistance mechanisms for model: {resistance_mechanisms}")
                        else:
                            logger.info("RGI output was empty")
                    else:
                        logger.warning(f"RGI output file not found at {rgi_txt}")
                        logger.debug(f"Contents of {Path(rgi_txt).parent}: {list(Path(rgi_txt).parent.iterdir())}")
                else:
                    logger.warning("RGI did not return an output prefix")
        except Exception as e:
            logger.error(f"Could not run RGI: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        os.unlink(temp_fasta)
        
        inc_type_cleaned = [re.sub(r"\([^)]*\)", "", item) for item in inc_types]
        
        all_features_for_model = (
            drug_classes +
            resistance_mechanisms +
            inc_type_cleaned +
            isolation_sources +
            conjugation_system
        )
        
        valid_features = [f for f in all_features_for_model if f in df.columns]
        if valid_features:
            df.loc[first_seq_id, valid_features] = 1
        
        prediction_proba = model.predict_proba(df)[0]
        prediction_class = model.predict(df)[0]
        
        labels = ['≤ Species', 'Genus/Family', 'Order/Class', '≥ Phylum']
        
        result = {
            'success': True,
            'sequence_id': first_seq_id,
            'sequence_length': len(plasmid_sequence[first_seq_id].seq),
            'detected_features': {
                'conjugation_systems': conjugation_system_detected,
                'inc_types': inc_types_detected,
                'antibiotic_resistance_genes': args_detected
            },
            'predictions': {
                label: float(prob)
                for label, prob in zip(labels, prediction_proba)
            },
            'predicted_class': labels[int(prediction_class)],
            'top_prediction': max(
                zip(labels, prediction_proba),
                key=lambda x: x[1]
            )[0],
            'top_probability': float(max(prediction_proba))
        }
        
        return result
    
    except Exception as e:
        logger.error(f"Error in predict_host_range: {str(e)}")
        return {'error': str(e), 'success': False}


@app.route('/')
def index():
    """Render the main page"""
    return render_template('index.html')


@app.route('/api/predict', methods=['POST'])
def api_predict():
    """API endpoint for predictions"""
    try:
        seq_input = request.form.get('sequence_input', '').strip()
        seq_source = request.form.get('sequence_source', 'text')
        
        isolation_sources = request.form.getlist('isolation_source')
        
        if seq_source == 'file':
            if 'sequence_file' not in request.files:
                return jsonify({'error': 'No file provided', 'success': False}), 400
            
            file = request.files['sequence_file']
            if file.filename == '':
                return jsonify({'error': 'No file selected', 'success': False}), 400
            
            if not file.filename.endswith(('.fna', '.fasta', '.fa', '.txt')):
                return jsonify({
                    'error': 'File must be FASTA format (.fna, .fasta, .fa)',
                    'success': False
                }), 400
            
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            fasta_content = filepath
        
        else:
            if not seq_input:
                return jsonify({
                    'error': 'Please provide a sequence',
                    'success': False
                }), 400
            fasta_content = seq_input
        
        if len(isolation_sources) > 1:
            return jsonify({
                'error': 'Please select only one isolation source',
                'success': False
            }), 400
         
        if not isolation_sources:
            isolation_sources = []
 
        result = predict_host_range(fasta_content, isolation_sources)
        
        # Clean up file if uploaded
        if seq_source == 'file' and os.path.exists(filepath):
            os.unlink(filepath)
        
        if result.get('success', False):
            return jsonify(result), 200
        else:
            return jsonify(result), 400
    
    except Exception as e:
        return jsonify({'error': f'Server error: {str(e)}', 'success': False}), 500


@app.route('/api/metadata', methods=['GET'])
def api_metadata():
    return jsonify({
        'isolation_sources': [
            'Aquatic animal',
            'Clinical',
            'Domestic animal',
            'Fungus',
            'Human',
            'Insect',
            'Plant',
            'Soil',
            'Wastewater',
            'Water',
            'Wild animal'
        ]
    })


# Error handlers
@app.errorhandler(404)
def not_found(error):
    """Handle 404 errors"""
    logger.warning(f"404 error: {request.path}")
    return jsonify({'error': 'Not found', 'success': False}), 404


@app.errorhandler(413)
def request_entity_too_large(error):
    """Handle file size errors"""
    logger.warning(f"413 error: File too large")
    return jsonify({
        'error': 'File is too large. Maximum size is 50 MB',
        'success': False
    }), 413


@app.errorhandler(500)
def internal_error(error):
    """Handle 500 errors"""
    logger.error(f"500 error: {str(error)}")
    return jsonify({
        'error': 'Internal server error',
        'success': False
    }), 500


# Application context
@app.before_request
def before_request():
    """Log incoming requests"""
    logger.debug(f"Request: {request.method} {request.path}")


@app.after_request
def after_request(response):
    """Log response"""
    logger.debug(f"Response: {response.status_code}")
    return response


if __name__ == '__main__':
    # Create logs directory
    os.makedirs('logs', exist_ok=True)
    
    logger.info(f"Starting Plasmid Host Range Predictor in {env} mode")
    logger.info(f"Debug mode: {app.debug}")
    
    # Run the application
    app.run(
        host='0.0.0.0',
        port=8000,
        debug=app.debug,
        use_reloader=app.debug
    )
