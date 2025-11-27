#!/bin/bash
# run_analysis.sh

echo "========================================="
echo "PSORIASIS-CROHN'S DISEASE ANALYSIS"
echo "Reproducing Li et al. 2025"
echo "========================================="

# Create directories
mkdir -p data/processed results figures logs

echo "1. Installing required packages..."
Rscript install_packages.R

echo "2. Running complete analysis pipeline..."
Rscript main_analysis.R

echo "3. Generating summary report..."
# This would generate a comprehensive report
Rscript -e "rmarkdown::render('scripts/summary_report.Rmd', output_file = '../results/analysis_summary.html')"

echo "========================================="
echo "ANALYSIS COMPLETED!"
echo "Check the following folders:"
echo "- results/ : Analysis results and tables"
echo "- figures/ : All generated plots"
echo "- logs/    : Analysis logs"
echo "========================================="
