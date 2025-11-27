# install_packages.R
# ==============================================================================
# INSTALL ALL REQUIRED PACKAGES FOR THE ANALYSIS
# ==============================================================================

# List of all required packages
required_packages <- c(
  # Bioconductor packages
  "BiocManager",
  "limma", "WGCNA", "clusterProfiler", "enrichplot", "DOSE",
  "org.Hs.eg.db", "GSVA", "GSEABase", "Seurat", "SingleCellExperiment",
  
  # CRAN packages for data manipulation
  "ggplot2", "pheatmap", "reshape2", "dplyr", "tidyr", "stringr",
  "readr", "tibble", "purrr",
  
  # CRAN packages for machine learning
  "ROCR", "caret", "xgboost", "randomForest", "e1071", "glmnet",
  "tidymodels", "SHAPforxgboost", "pROC", "viridis",
  
  # CRAN packages for visualization
  "patchwork", "RColorBrewer", "corrplot", "ggpubr", "ggrepel",
  
  # CRAN packages for networks
  "igraph", "network", "sna",
  
  # Utility packages
  "here", "fs", "openxlsx"
)

# Function to install and load packages
install_and_load <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      # Install from CRAN
      install.packages(pkg, dependencies = TRUE)
      
      # Load the package
      if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
        # Try Bioconductor if CRAN fails
        if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
          BiocManager::install(pkg)
          library(pkg, character.only = TRUE)
        }
      }
    }
  }
}

# Install Bioconductor first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install and load all packages
cat("Installing and loading required packages...\n")
install_and_load(required_packages)

# Verify installations
cat("\n=== PACKAGE INSTALLATION VERIFICATION ===\n")
for(pkg in required_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗", pkg, "failed to load\n")
  }
}

cat("\n=== SETUP COMPLETED ===\n")
cat("All packages installed and ready for analysis!\n")
