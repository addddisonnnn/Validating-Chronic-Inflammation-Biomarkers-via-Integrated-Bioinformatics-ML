# main_analysis.R
# ==============================================================================
# MAIN CONTROLLER SCRIPT FOR PSORIASIS-CROHN'S DISEASE ANALYSIS
# This script coordinates the entire analysis pipeline
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Enable parallel processing for faster computation
enableWGCNAThreads(nThreads = 4)

# Create output directories
create_directories <- function() {
  dirs <- c("data/processed", "results", "figures", "logs")
  for(dir in dirs) {
    if(!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Main analysis function
run_complete_analysis <- function() {
  cat("=========================================\n")
  cat("PSORIASIS-CROHN'S DISEASE ANALYSIS PIPELINE\n")
  cat("Reproducing Li et al. 2025\n")
  cat("=========================================\n\n")
  
  # Create directories
  create_directories()
  
  # Start analysis log
  sink("logs/analysis_log.txt", append = FALSE)
  cat("Analysis started at:", Sys.time(), "\n")
  sink()
  
  tryCatch({
    # ==========================================================================
    # STEP 1: Data Preprocessing
    # ==========================================================================
    cat("Step 1: Data Preprocessing...\n")
    source("scripts/01_data_preprocessing.R")
    processed_data <- preprocess_geo_data()
    
    # ==========================================================================
    # STEP 2: Differential Expression Analysis
    # ==========================================================================
    cat("Step 2: Differential Expression Analysis...\n")
    source("scripts/02_differential_expression.R")
    deg_results <- perform_differential_expression(processed_data)
    
    # ==========================================================================
    # STEP 3: Gene Set Enrichment Analysis
    # ==========================================================================
    cat("Step 3: Gene Set Enrichment Analysis...\n")
    source("scripts/03_gsea_analysis.R")
    gsea_results <- perform_gsea_analysis(deg_results)
    
    # ==========================================================================
    # STEP 4: WGCNA Network Analysis
    # ==========================================================================
    cat("Step 4: WGCNA Network Analysis...\n")
    source("scripts/04_wgcna_analysis.R")
    wgcna_results <- perform_wgcna_analysis(processed_data, deg_results)
    
    # ==========================================================================
    # STEP 5: Functional Enrichment Analysis
    # ==========================================================================
    cat("Step 5: Functional Enrichment Analysis...\n")
    source("scripts/05_functional_enrichment.R")
    enrichment_results <- perform_functional_enrichment(wgcna_results$shared_genes)
    
    # ==========================================================================
    # STEP 6: PPI Network Construction
    # ==========================================================================
    cat("Step 6: PPI Network Analysis...\n")
    source("scripts/06_ppi_network.R")
    ppi_results <- build_ppi_network(wgcna_results$shared_genes)
    
    # ==========================================================================
    # STEP 7: Machine Learning Biomarker Selection
    # ==========================================================================
    cat("Step 7: Machine Learning Biomarker Selection...\n")
    source("scripts/07_machine_learning.R")
    ml_results <- perform_machine_learning(processed_data, ppi_results$hub_genes)
    
    # ==========================================================================
    # STEP 8: Diagnostic Validation
    # ==========================================================================
    cat("Step 8: Diagnostic Validation...\n")
    source("scripts/08_diagnostic_validation.R")
    validation_results <- validate_biomarkers(processed_data, ml_results$final_biomarkers)
    
    # ==========================================================================
    # STEP 9: Immune Infiltration Analysis
    # ==========================================================================
    cat("Step 9: Immune Infiltration Analysis...\n")
    source("scripts/09_immune_infiltration.R")
    immune_results <- analyze_immune_infiltration(processed_data, ml_results$final_biomarkers)
    
    # ==========================================================================
    # STEP 10: Single-Cell Analysis (if data available)
    # ==========================================================================
    cat("Step 10: Single-Cell Analysis...\n")
    source("scripts/10_single_cell_analysis.R")
    sc_results <- analyze_single_cell_data()
    
    # ==========================================================================
    # STEP 11: Drug Repurposing
    # ==========================================================================
    cat("Step 11: Drug Repurposing...\n")
    source("scripts/11_drug_repurposing.R")
    drug_results <- perform_drug_repurposing(ml_results$final_biomarkers)
    
    # ==========================================================================
    # STEP 12: Molecular Docking (Simplified)
    # ==========================================================================
    cat("Step 12: Molecular Docking Analysis...\n")
    source("scripts/12_molecular_docking.R")
    docking_results <- perform_molecular_docking(drug_results$top_drugs)
    
    # ==========================================================================
    # Generate Summary Report
    # ==========================================================================
    generate_summary_report(list(
      deg_results = deg_results,
      wgcna_results = wgcna_results,
      ml_results = ml_results,
      validation_results = validation_results,
      drug_results = drug_results
    ))
    
    cat("\n=========================================\n")
    cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
    cat("=========================================\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in analysis:", e$message, "\n")
    sink("logs/error_log.txt", append = FALSE)
    cat("Error occurred at:", Sys.time(), "\n")
    cat("Error message:", e$message, "\n")
    traceback()
    sink()
    return(FALSE)
  })
}

# Run the analysis
if(!interactive()) {
  success <- run_complete_analysis()
  quit(save = "no", status = ifelse(success, 0, 1))
}
