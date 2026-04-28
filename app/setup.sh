#!/bin/bash
# Quick start script for Plasmid Host Range Predictor Web Application

set -e

echo "========================================"
echo "Plasmid Host Range Predictor - Setup"
echo "========================================"
echo ""

# Check Python version
echo "Checking Python version..."
python_version=$(python3 --version 2>&1)
echo "✓ Found: $python_version"
echo ""

# Check for required system tools
echo "Checking for required system tools..."

if ! command -v prodigal &> /dev/null; then
    echo "⚠ prodigal not found. Installing..."
    sudo apt-get update
    sudo apt-get install -y prodigal
else
    echo "✓ prodigal is installed"
fi

if ! command -v hmmsearch &> /dev/null; then
    echo "⚠ hmmer not found. Installing..."
    sudo apt-get update
    sudo apt-get install -y hmmer
else
    echo "✓ hmmer is installed"
fi

echo ""
echo "Installing Python dependencies..."
pip install -q -r requirements.txt
echo "✓ Dependencies installed"
echo ""

# Check for model file
if [ ! -f "plaspredict_model.pkl" ]; then
    echo "⚠ Warning: plaspredict_model.pkl not found in project root"
    echo "  The web application will fail without this model file"
else
    echo "✓ Model file found: plaspredict_model.pkl"
fi

# Check for HMM models
if [ ! -d "files/conjscan_models" ]; then
    echo "⚠ Warning: HMM models directory not found"
    echo "  Conjugation system detection will not work"
else
    num_models=$(ls files/conjscan_models/*.hmm 2>/dev/null | wc -l)
    echo "✓ Found $num_models HMM models"
fi

echo ""
echo "========================================"
echo "Setup complete!"
echo "========================================"
echo ""
echo "To start the web application, run:"
echo ""
echo "    python3 app.py"
echo ""
echo "Then open your browser to:"
echo "    http://localhost:5000"
echo ""
echo "For production deployment, use:"
echo "    gunicorn -w 4 -b 0.0.0.0:5000 app:app"
echo ""
