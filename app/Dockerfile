# ── Stage 1: build ──────────────────────────────────────────────────────────
FROM python:3.11-slim AS builder

WORKDIR /build

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    hmmer \
    ncbi-blast+ \
    prodigal \
    && rm -rf /var/lib/apt/lists/*

# Install RGI (Resistance Gene Identifier) from GitHub
RUN pip install git+https://github.com/arpcard/rgi.git

# Install PlasmidFinder
RUN pip install plasmidfinder

# Copy application code
COPY app/requirements.txt /build/

# Install any other Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Create a dedicated staging directory for the database in builder
RUN mkdir -p /build/plasmidfinder_db

# Clone or download the database into the staging directory. 
RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git /build/plasmidfinder_db

# ── Stage 2: runtime ─────────────────────────────────────────────────────────
FROM python:3.11-slim

# ── Non-root user ────────────────────────────────────────────────────────────
ENV USER=appuser
ENV HOME=/home/$USER
RUN useradd -m -u 1000 appuser

# Install runtime system dependencies (including git for plasmidfinder)
RUN apt-get update && apt-get install -y \
    git \
    prodigal \
    hmmer \
    ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# Copy Python packages from builder
COPY --from=builder /usr/local /usr/local

# Create application directory structure (including the required DB path)
RUN mkdir -p $HOME/app $HOME/data $HOME/logs $HOME/uploads $HOME/app/db/plasmidfinder/database

# FIX: Copy the downloaded database explicitly from the builder stage
COPY --from=builder --chown=$USER:$USER /build/plasmidfinder_db/ $HOME/app/db/plasmidfinder/database/

# Copy application code
COPY --chown=$USER:$USER app/ $HOME/app/

# Double check all path ownership is clean
RUN chown -R $USER:$USER $HOME

# Verify critical files exist
RUN ls -la $HOME/app/models/ || echo "WARNING: models directory missing" && \
    ls -la $HOME/app/localDB/ || echo "WARNING: localDB directory missing" && \
    ls -la $HOME/app/db/plasmidfinder/database/ || echo "WARNING: PlasmidFinder DB missing"

USER $USER
WORKDIR $HOME/app

ENV PYTHONUNBUFFERED=1
ENV FLASK_ENV=production

EXPOSE 5000
ENTRYPOINT ["./start-script.sh"]