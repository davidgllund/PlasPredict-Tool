#!/usr/bin/env python3
"""
Flask web application for Plasmid Host Range Prediction
"""
import os
import sys
import tempfile
import logging
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

# Import configuration
from config import get_config

# Create Flask app
app = Flask(__name__)

# Load configuration based on environment
env = os.environ.get('FLASK_ENV', 'development')
config = get_config(env)
app.config.from_object(config)

# Setup logging
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent

# Create uploads directory if it doesn't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Load model
logger.info(f"Loading model from {app.config['MODEL_PATH']}")
if not Path(app.config['MODEL_PATH']).exists():
    raise FileNotFoundError(f"Model not found at {app.config['MODEL_PATH']}")
model = joblib.load(str(app.config['MODEL_PATH']))
logger.info("Model loaded successfully")

# HMM models for conjugation
HMM_PATH = app.config['HMM_PATH']
logger.info(f"HMM models path: {HMM_PATH}")


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


def calc_kmer_distributions(genome, k, possible_kmers):
    """Calculate k-mer distribution for entire genome"""
    genome_kmers = []
    for seq_entry in genome.values():
        genome_kmers.extend(get_kmers(seq_entry.seq, k))
    
    distribution = get_kmer_distribution(genome_kmers, possible_kmers)
    return distribution


def annotate_conjugation_system(filepath, models):
    """Identify conjugation systems in plasmid using HMM models"""
    conjugation_system = []
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        proteins = tmp / "proteins.faa"
        
        try:
            # Gene prediction with prodigal
            subprocess.run(
                ["prodigal", "-i", filepath, "-a", str(proteins), "-p", "meta"],
                check=True,
                capture_output=True
            )
            
            # Search for conjugation markers
            for hmm in models:
                hits = subprocess.run(
                    f'hmmsearch -E 0.0000000001 "{hmm}" "{proteins}"',
                    shell=True,
                    capture_output=True,
                    text=True
                )
                
                if "[No hits detected that satisfy reporting thresholds]" not in hits.stdout:
                    element = hmm.split('/')[-1].split('.')[0].split('_')[0]
                    if element not in conjugation_system:
                        conjugation_system.append(element)
        
        except Exception as e:
            print(f"Warning: Could not run conjugation annotation: {e}")
    
    return conjugation_system


def predict_host_range(fasta_content, inc_types, drug_classes, 
                       resistance_mechanisms, isolation_sources):
    """
    Make host range prediction for a plasmid sequence
    
    Args:
        fasta_content: FASTA sequence as string or file path
        inc_types: List of incompatibility types selected
        drug_classes: List of drug classes selected
        resistance_mechanisms: List of resistance mechanisms selected
        isolation_sources: List of isolation sources selected
    
    Returns:
        Dictionary with prediction results
    """
    try:
        # Parse FASTA content
        if os.path.isfile(fasta_content):
            plasmid_sequence = SeqIO.to_dict(
                SeqIO.parse(fasta_content, 'fasta')
            )
        else:
            # Create temporary FASTA file from string
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as f:
                f.write(fasta_content)
                temp_fasta = f.name
            plasmid_sequence = SeqIO.to_dict(
                SeqIO.parse(temp_fasta, 'fasta')
            )
            os.unlink(temp_fasta)
        
        if not plasmid_sequence:
            return {'error': 'No valid sequences found in FASTA file'}
        
        # Create feature dataframe
        df = pd.DataFrame(
            0,
            index=plasmid_sequence.keys(),
            columns=model.get_booster().feature_names
        )
        
        # Calculate k-mer distributions
        all_kmers = generate_possible_kmers(3)
        kmer_distributions = calc_kmer_distributions(plasmid_sequence, 3, all_kmers)
        kmer_df = pd.DataFrame(
            [kmer_distributions],
            index=plasmid_sequence.keys(),
            columns=all_kmers
        )
        
        # Add k-mer features
        kmer_cols = list(set(kmer_df.columns) & set(df.columns))
        df[kmer_cols] = kmer_df[kmer_cols]
        
        # Add plasmid size
        first_seq_id = list(plasmid_sequence.keys())[0]
        df.loc[first_seq_id, 'size'] = len(plasmid_sequence[first_seq_id].seq)
        
        # Annotate conjugation systems
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as f:
            for seq_id, seq_record in plasmid_sequence.items():
                f.write(f">{seq_id}\n{str(seq_record.seq)}\n")
            temp_fasta = f.name
        
        hmm_files = glob.glob(str(HMM_PATH / '*'))
        conjugation_system = annotate_conjugation_system(temp_fasta, hmm_files)
        os.unlink(temp_fasta)
        
        # Clean inc_type names (remove parentheses)
        inc_type_cleaned = [re.sub(r"\([^)]*\)", "", item) for item in inc_types]
        
        # Set feature values for selected categories
        all_features = (
            drug_classes +
            resistance_mechanisms +
            inc_type_cleaned +
            isolation_sources +
            conjugation_system
        )
        
        # Only set features that exist in the model
        valid_features = [f for f in all_features if f in df.columns]
        if valid_features:
            df.loc[first_seq_id, valid_features] = 1
        
        # Make prediction
        prediction_proba = model.predict_proba(df)[0]
        prediction_class = model.predict(df)[0]
        
        labels = ['≤ Species', 'Genus/Family', 'Order/Class', '≥ Phylum']
        
        result = {
            'success': True,
            'sequence_id': first_seq_id,
            'sequence_length': len(plasmid_sequence[first_seq_id].seq),
            'conjugation_systems': conjugation_system,
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
        return {'error': str(e), 'success': False}


@app.route('/')
def index():
    """Render the main page"""
    return render_template('index.html')


@app.route('/api/predict', methods=['POST'])
def api_predict():
    """API endpoint for predictions"""
    try:
        # Parse request
        seq_input = request.form.get('sequence_input', '').strip()
        seq_source = request.form.get('sequence_source', 'text')
        
        inc_types = request.form.getlist('inc_type')
        drug_classes = request.form.getlist('drug_class')
        resistance_mechanisms = request.form.getlist('resistance_mechanism')
        isolation_sources = request.form.getlist('isolation_source')
        
        # Get sequence
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
        
        else:  # text input
            if not seq_input:
                return jsonify({
                    'error': 'Please provide a sequence',
                    'success': False
                }), 400
            fasta_content = seq_input
        
        # Validate selections
        # Inc types, drug classes, and resistance mechanisms can be empty
        # Only isolation sources is required (exactly one)
        if not isolation_sources or len(isolation_sources) != 1:
            return jsonify({
                'error': 'Please select exactly one isolation source',
                'success': False
            }), 400
        
        # Run prediction
        result = predict_host_range(
            fasta_content,
            inc_types,
            drug_classes,
            resistance_mechanisms,
            isolation_sources
        )
        
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
    """Get available options for dropdowns"""
    # Load inc types from file
    inc_types = []
    try:
        with open('inc_types.txt', 'r') as f:
            inc_types = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        # Fallback to hardcoded list if file not found
        inc_types = [
            'IncQ1',
            'IncHI1B',
            'IncHI1A',
            'IncC',
            'IncF',
            'IncX',
            'IncN',
            'IncP',
            'IncW',
            'IncM',
            'IncU',
            'IncT',
            'IncA',
            'IncB',
            'IncD',
            'IncE',
            'IncG',
            'IncH',
            'IncJ',
            'IncK',
            'IncL',
            'IncO',
            'IncS',
            'IncV',
            'IncZ'
        ]
    
    return jsonify({
        'inc_types': inc_types,
        'drug_classes': [
            'aminocoumarin antibiotic',
            'aminoglycoside antibiotic',
            'antibacterial free fatty acids',
            'bicyclomycin-like antibiotic',
            'carbapenem',
            'cephamycin',
            'diaminopyrimidine antibiotic',
            'disinfecting agents and antiseptics',
            'elfamycin antibiotic',
            'fluoroquinolone antibiotic',
            'glycopeptide antibiotic',
            'glycylcycline',
            'isoniazid-like antibiotic',
            'lincosamide antibiotic',
            'macrolide antibiotic',
            'monobactam',
            'mupirocin-like antibiotic',
            'nitrofuran antibiotic',
            'orthosomycin antibiotic',
            'oxazolidinone antibiotic',
            'penam',
            'penem',
            'peptide antibiotic',
            'phenicol antibiotic',
            'phosphonic acid antibiotic',
            'pleuromutilin antibiotic',
            'polyamine antibiotic',
            'pyrazine antibiotic',
            'rifamycin antibiotic',
            'salicylic acid antibiotic',
            'streptogramin antibiotic',
            'streptogramin A antibiotic',
            'streptogramin B antibiotic',
            'sulfonamide antibiotic',
            'tetracycline antibiotic',
            'thioamide antibiotic'
        ],
        'resistance_mechanisms': [
            'antibiotic efflux',
            'antibiotic inactivation',
            'antibiotic target replacement',
            'antibiotic target protection',
            'antibiotic target alteration',
            'reduced permeability to antibiotic'
        ],
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
        port=5000,
        debug=app.debug,
        use_reloader=app.debug
    )
