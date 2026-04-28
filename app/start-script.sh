#!/bin/bash
set -e

# Set production environment if not already set
export FLASK_ENV=${FLASK_ENV:-production}

# Run application with Gunicorn
exec gunicorn \
    --bind=0.0.0.0:8000 \
    --workers=2 \
    --timeout=120 \
    --access-logfile=- \
    --error-logfile=- \
    'app:app'
