#!/bin/bash
# Main script to download all GEO datasets

# Set up project directory
PROJECT_DIR="raw_data"
mkdir -p $PROJECT_DIR

echo "Starting Phase 1 - Data Acquisition"
echo "===================================="

# Make all scripts executable
chmod +x download_*.sh

# Download all datasets
echo "Downloading GSE13355 (Psoriasis)..."
./download_GSE13355.sh

echo "Downloading GSE14905 (Psoriasis)..."
./download_GSE14905.sh

echo "Downloading GSE75214 (CD)..."
./download_GSE75214.sh

echo "Downloading GSE102133 (CD)..."
./download_GSE102133.sh

echo "All downloads completed!"
echo "Data saved in: $PROJECT_DIR"