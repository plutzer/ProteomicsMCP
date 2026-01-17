# CPTAC Proteomics Explorer - Dockerfile
#
# Self-contained build with pre-downloaded CPTAC data included.
# No external data mounts required.

FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY local_data_loader.py .
COPY cptac_proteomics.py .
COPY GUI/ ./GUI/

# Copy CPTAC data (pre-downloaded, ~500MB)
COPY data/ ./data/

# Copy datasets directory (for PSP data)
COPY Datasets/ ./Datasets/

ENV PYTHONUNBUFFERED=1

# Start the GUI
EXPOSE 3838
CMD ["python3", "-u", "/app/GUI/cptac_explorer_app.py"]
