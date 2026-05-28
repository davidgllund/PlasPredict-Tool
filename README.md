# Plasmid Host Range Predictor (PlasPart)

A machine learning-powered web application for predicting the host range of plasmids using k-mer composition analysis, conjugation system detection, and automated feature identification.

## Overview

PlasPart analyzes plasmid DNA sequences to predict their potential bacterial host range based on genomic features:
- **K-mer composition analysis** - Uses k-mer frequency distribution (k=3) extracted from plasmid sequences
- **Automated feature detection** - Automatically identifies incompatibility types (via Plasmidfinder), drug classes and resistance mechanisms (via RGI)
- **Machine learning prediction** - XGBoost-based model trained on known plasmid-host relationships
- **Conjugation system detection** - Identifies conjugative elements using HMM models from ConjScan

## Features

- 🧬 **Sequence Upload** - Support for FASTA/FASTA-like formats (.fna, .fasta, .fa, .txt)
- 🤖 **ML Prediction** - XGBoost model for host range prediction
- 🔍 **Conjugation Detection** - HMM-based identification of conjugative systems
- 📊 **Result Visualization** - Interactive display of predictions and detected systems
- 🐳 **Docker Deployment** - Container-ready with production-grade Gunicorn server
- ☁️ **Cloud-Ready** - GitHub Actions CI/CD for automated DockerHub publishing

## Quick Start

### Local Development

#### Prerequisites
- Python 3.11+
- System dependencies: `prodigal`, `hmmer`

#### Installation

```bash
# Clone the repository
git clone https://github.com/davidgllund/plaspart.git
cd plaspart

# Install system dependencies (Ubuntu/Debian)
sudo apt-get install -y prodigal hmmer

# Install Python dependencies
pip install -r app/requirements.txt
```

#### Running the Application

```bash
cd app
python app.py
```

The web application will be available at `http://localhost:5000`

### Docker Deployment

#### Build Locally

```bash
docker build -f app/Dockerfile -t plaspart:latest .
docker run -p 8000:8000 plaspart:latest
```

Access at `http://localhost:8000`

#### Using DockerHub Image

```bash
docker pull davidgllund/plaspart:latest
docker run -p 8000:8000 davidgllund/plaspart:latest
```

#### Docker Compose

```bash
cd app
docker-compose up
```

## Project Structure

```
plaspart/
├── app/
│   ├── app.py                    # Main Flask application
│   ├── config.py                 # Configuration management
│   ├── requirements.txt           # Python dependencies
│   ├── Dockerfile                 # Multi-stage Docker build
│   ├── docker-compose.yml        # Docker Compose configuration
│   ├── start-script.sh           # Container entrypoint (Gunicorn)
│   ├── models/
│   │   ├── plaspredict_model.pkl # XGBoost trained model
│   │   └── conjscan_models/      # HMM models for conjugation detection
│   ├── static/
│   │   ├── style.css            # Application styling
│   │   ├── script.js            # Frontend JavaScript
│   │   └── plaspredict_logo.png # Logo
│   ├── templates/
│   │   └── index.html           # Web interface HTML
│   ├── uploads/                 # User-uploaded sequence files
│   ├── logs/                    # Application logs
│   └── tmp/                     # Temporary analysis files
├── .github/
│   └── workflows/
│       └── docker-build-publish.yml  # GitHub Actions CI/CD
├── DOCKER_PUBLISH_SETUP.md      # Docker/DockerHub setup guide
└── README.md                    # This file
```

## Configuration

### Environment Variables

```bash
# Flask environment (development/production)
export FLASK_ENV=production

# Secret key (change in production!)
export SECRET_KEY=your-secret-key-here
```

### Application Settings

Edit `app/config.py` to modify:
- Maximum file upload size (default: 50MB)
- Maximum sequence length (default: 1MB)
- Prediction timeout (default: 120 seconds)
- Logging level and paths
- Model and HMM paths

## Usage

### Web Interface

1. **Upload a sequence file** - Select a FASTA/nucleotide file from your computer
2. **Submit for analysis** - Click the "Predict" button
3. **View results**:
   - Host range predictions with confidence scores
   - Detected conjugative elements and systems
   - Sequence statistics and analysis details

### API Endpoints

#### POST `/predict`
Submit a sequence for analysis.

**Request:**
```bash
curl -X POST -F "file=@plasmid.fasta" http://localhost:8000/predict
```

**Response:**
```json
{
  "success": true,
  "predictions": {
    "Escherichia": 0.92,
    "Bacillus": 0.45,
    "Pseudomonas": 0.78
  },
  "conjugation_systems": [
    {
      "system": "B_traE",
      "type": "Type B",
      "genes": ["traE", "traF"]
    }
  ],
  "sequence_length": 5280,
  "gc_content": 0.52
}
```

#### GET `/get-inc-types`
Retrieve all supported incompatibility (Inc) types.

**Response:**
```json
{
  "inc_types": ["IncA", "IncB", "IncC", ...]
}
```

## Model Details

### XGBoost Predictor

- **Training data**: Curated plasmid sequences with known host ranges
- **Features**: K-mer frequency distributions (k=1 to k=6)
- **Output**: Probability scores for bacterial host genera

### Conjugation System Detection

Uses HMM models from ConjScan database:
- **Type B** conjugation systems (tra genes)
- **Type C/G/F/FA/FATA/I** systems
- Gene prediction via Prodigal
- HMM matching via HMMER

## Docker & CI/CD

### GitHub Actions Workflow

Automatically builds and publishes Docker images on:
- **Push to main** - Tags: `latest`, commit SHA
- **Version tags** (e.g., `v1.0.0`) - Semantic versioning tags
- **Pull requests** - Build only (no push)

See [DOCKER_PUBLISH_SETUP.md](DOCKER_PUBLISH_SETUP.md) for detailed setup instructions.

### Build Configuration

The Dockerfile uses a multi-stage build for optimal image size:
1. **Builder stage** - Compiles Python dependencies and system tools
2. **Runtime stage** - Minimal image with only necessary components

Key features:
- Non-root user (UID 1000) for security
- Gunicorn production server
- System dependencies: prodigal, hmmer
- Caches Python packages

## Troubleshooting

### 502 Bad Gateway Error
**Symptoms:** Application shows 502 error when accessed
**Solution:** 
- Ensure Gunicorn is running on port 8000
- Check logs: `docker logs <container-id>`
- Verify all dependencies are installed

### Conjugation Systems Not Detected
**Symptoms:** Prediction results show empty conjugation systems
**Possible causes:**
- HMM models directory not mounted/copied correctly
- HMMER or Prodigal not installed in container
- Sequence too short or too low quality

**Solutions:**
- Verify model files exist: `/home/appuser/app/models/conjscan_models/`
- Check container logs for errors
- Test with a longer plasmid sequence

### Model Loading Fails
**Error:** `FileNotFoundError: Model not found at...`
**Solution:**
- Ensure `plaspredict_model.pkl` exists in `app/models/`
- Verify file permissions are readable
- Check paths in `app/config.py`

### File Upload Fails
**Error:** Request entity too large or file not accepted
**Solutions:**
- Increase `MAX_CONTENT_LENGTH` in `config.py` (max 50MB by default)
- Ensure file extension is in `ALLOWED_EXTENSIONS` (.fna, .fasta, .fa, .txt)
- Verify file size is under limit: `ls -lh your_file.fasta`

## Development

### Setting Up Development Environment

```bash
# Create virtual environment
python3.11 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r app/requirements.txt

# Install development tools
pip install pytest pytest-cov flake8

# Run tests
pytest

# Run linting
flake8 app/
```

### Building Docker Image for Testing

```bash
# Build with tag
docker build -f app/Dockerfile -t plaspredict:dev .

# Run with volume mount for development
docker run -v $(pwd)/app:/home/appuser/app -p 8000:8000 plaspredict:dev
```

## Performance Considerations

- **Prediction timeout**: 120 seconds (configurable)
- **Max sequence length**: 1MB (1,000,000 bp)
- **Gunicorn workers**: 2 (configurable for your hardware)
- **File upload size**: 50MB maximum
- **Memory usage**: ~500MB-1GB per worker

For high-throughput deployments, consider:
- Increasing worker count based on CPU cores
- Using Kubernetes for horizontal scaling
- Adding Redis for result caching
- Load balancing with Nginx

## Dependencies

### System Requirements
- Python 3.11 or higher
- Prodigal (gene prediction)
- HMMER3 (HMM-based sequence matching)

### Python Packages
- **Flask** 2.3.3 - Web framework
- **Gunicorn** 21.2.0 - WSGI HTTP server
- **XGBoost** 2.0.0 - Machine learning model
- **BioPython** 1.81 - Biological sequence handling
- **NumPy** 1.24.3 - Numerical computing
- **Pandas** 2.0.3 - Data analysis
- **joblib** 1.3.2 - Model serialization
- **Werkzeug** 2.3.7 - WSGI utilities

## Security Notes

### Production Deployment

1. **Change SECRET_KEY** - Update in environment variables
2. **Enable HTTPS** - Use reverse proxy (Nginx, Caddy)
3. **File Uploads** - Secure temporary directory handling
4. **Model Protection** - Keep trained models in private storage
5. **Logging** - Monitor for suspicious activities
6. **Rate Limiting** - Consider adding request throttling
7. **CORS** - Configure if serving API to other domains

### Container Security

- Non-root user (appuser, UID 1000)
- Read-only filesystems where possible
- No hardcoded credentials
- Regular dependency updates via GitHub Actions

## License

[Add your license here]

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Support

For issues, questions, or suggestions:
- Open an [GitHub issue](https://github.com/davidgllund/plaspredict/issues)
- Check existing documentation in the repository
- Review troubleshooting section above

## Citation

If you use PlasPredict in your research, please cite:

```
[Add citation information here]
```

## Authors

- **David G. Lund** - Initial development and deployment

## Acknowledgments

- ConjScan database for HMM models
- BioPython community
- Flask and XGBoost teams

## Changelog

### Version 1.0.0 (2026-04-28)
- Initial release with Docker support
- GitHub Actions CI/CD pipeline
- Multi-stage Dockerfile optimization
- Gunicorn production server
- Comprehensive documentation

---

**Last Updated:** April 28, 2026
**Repository:** https://github.com/davidgllund/plaspredict
**DockerHub:** https://hub.docker.com/r/davidgllund/plaspredict
