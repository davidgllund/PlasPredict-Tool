# ── Stage 1: build ──────────────────────────────────────────────────────────
FROM python:3.11-slim AS builder

WORKDIR /build

# Install system dependencies
RUN apt-get update && apt-get install -y \
    prodigal \
    hmmer \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY app/requirements.txt .
RUN pip install --prefix=/install --no-cache-dir -r requirements.txt
 
# ── Stage 2: runtime ─────────────────────────────────────────────────────────
FROM python:3.11-slim

# ── Non-root user (UID 1000 required by SciLifeLab Serve) ────────────────────
ENV USER=appuser
ENV HOME=/home/$USER
RUN useradd -m -u 1000 appuser

# Copy application files
COPY --from=builder /install /usr/local

# Copy application code
COPY app/ $HOME/app/

# Copy start script into WORKDIR — Serve calls it as ./start-script.sh
COPY app/start-script.sh $HOME/app/start-script.sh

# Create the data directory that Serve will mount the volume onto
# Must exist in the image before Serve can bind-mount over it
RUN mkdir -p $HOME/data \
    && chmod +x $HOME/app/start-script.sh \
    && chown -R $USER:$USER $HOME

USER appuser

WORKDIR $HOME/app

# ── Gunicorn / Dash settings ──────────────────────────────────────────────────
ENV PYTHONUNBUFFERED=1
# Serve only allows ports 3000-9999; 8000 is used in all Serve examples
ENV GUNICORN_CMD_ARGS="--bind=0.0.0.0:5000 --workers=2 --timeout=120 --forwarded-allow-ips='*' --access-logfile=-"

USER $USER
EXPOSE 5000

ENTRYPOINT ["./start-script.sh"]
