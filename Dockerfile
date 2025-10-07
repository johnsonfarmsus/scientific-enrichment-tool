FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for pymatgen and RDKit
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy API service
COPY enrichment_api.py .

# Expose API port
EXPOSE 8099

# Run the API service
CMD ["python", "enrichment_api.py"]
