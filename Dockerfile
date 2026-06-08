# Synthetic Biology Metabolic Pathway Pipeline - Docker Image
FROM python:3.9-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Set work directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY *.py ./

# Create directories for output and logs
RUN mkdir -p /app/pipeline_output /app/logs

# Set volume mounts for persistent data
VOLUME ["/app/pipeline_output", "/app/logs"]

# Default command (can be overridden)
CMD ["python", "main_pipeline_runner.py", "--help"]

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import main_pipeline_runner; print('OK')" || exit 1

# Metadata
LABEL maintainer="synbio-pipeline-team" \
      version="1.0.0" \
      description="Industrial-grade metabolic pathway design pipeline for synthetic biology" \
      organism.support="ecoli,scerevisiae,bubtilis,cglutamicum,pputida" \
      molecule.support="lycopene,vanillin,lysine,riboflavin,pha"
