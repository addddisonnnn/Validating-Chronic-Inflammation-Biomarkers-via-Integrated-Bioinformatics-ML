#!/bin/bash
# download_GSE14905.sh

PROJECT_DIR="raw_data"
GSE="GSE14905"
OUTPUT_DIR="$PROJECT_DIR/$GSE"

mkdir -p $OUTPUT_DIR

echo "Downloading $GSE..."

# Download using wget
wget -q -O $OUTPUT_DIR/${GSE}_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=${GSE}&format=file"

if [ -f $OUTPUT_DIR/${GSE}_RAW.tar ]; then
    echo "Extracting $GSE..."
    tar -xf $OUTPUT_DIR/${GSE}_RAW.tar -C $OUTPUT_DIR/
    rm $OUTPUT_DIR/${GSE}_RAW.tar
    echo "$GSE download and extraction completed"
else
    echo "Failed to download $GSE"
fi