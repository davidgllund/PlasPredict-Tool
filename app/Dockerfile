# ── Stage 1: build ──────────────────────────────────────────────────────────
FROM python:3.11-slim AS builder

# Install system dependencies for bioinformatics tools
RUN apt-get update && apt-get install -y \
    hmmer \
    prodigal \
    blast2 \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# ── Stage 2: runtime ─────────────────────────────────────────────────────────
FROM python:3.11-slim

# ── Non-root user (UID 1000 required by SciLifeLab Serve) ────────────────────
ENV USER=appuser
ENV HOME=/home/$USER
RUN useradd -m -u 1000 appuser

# Install runtime dependencies (without build tools)
RUN apt-get update && apt-get install -y \
    hmmer \
    prodigal \
    blast2 \
    && rm -rf /var/lib/apt/lists/*

# Copy Python packages from builder
COPY --from=builder /usr/local /usr/local

# Copy application code from context root
COPY . $HOME/app/

# Create the data directory that Serve will mount the volume onto
RUN mkdir -p $HOME/data \
    && chmod +x $HOME/app/start-script.sh \
    && chown -R $USER:$USER $HOME

USER appuser

WORKDIR $HOME/app

# ── Gunicorn / Dash settings ──────────────────────────────────────────────────
ENV PYTHONUNBUFFERED=1
ENV GUNICORN_CMD_ARGS="--bind=0.0.0.0:5000 --workers=2 --timeout=120 --forwarded-allow-ips='*' --access-logfile=-"

EXPOSE 5000

ENTRYPOINT ["./start-script.sh"]