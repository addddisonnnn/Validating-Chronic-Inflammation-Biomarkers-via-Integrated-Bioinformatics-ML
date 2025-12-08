# scripts/01_data_preprocessing.R
# STEP 1: DATA PREPROCESSING
# Load, normalize, and quality control of GEO datasets

preprocess_geo_data <- function() {
  cat("Starting data preprocessing...\n")
  
  # Sample information based on Li et al. 2025
  sample_info <- list(
    psoriasis = list(
      GSE13355 = list(patients = 58, controls = 64),
      GSE14905 = list(patients = 33, controls = 21)
    ),
    crohns = list(
      GSE75214 = list(patients = 67, controls = 11),
      GSE102133 = list(patients = 65, controls = 12)
    )
  )
  
  # Function to simulate/load dataset (replace with actual data loading)
  load_dataset <- function(dataset_name, n_patients, n_controls) {
    cat("Processing", dataset_name, "...\n")
    
    # In practice, you would load actual data here:
    # expr_matrix <- read.csv(paste0("data/raw/", dataset_name, "/expression_matrix.csv"))
    
    # For demonstration, create simulated data
    set.seed(123)
    n_genes <- 20000
    n_samples <- n_patients + n_controls
    
    # Create simulated expression matrix
    expr_matrix <- matrix(
      rnorm(n_genes * n_samples, mean = 8, sd = 2),
      nrow = n_genes,
      ncol = n_samples
    )
    
    # Add some disease-specific signals
    disease_genes <- sample(1:n_genes, 1000)
    expr_matrix[disease_genes, 1:n_patients] <- 
      expr_matrix[disease_genes, 1:n_patients] + rnorm(1000 * n_patients, mean = 1, sd = 0.5)
    
    rownames(expr_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(expr_matrix) <- paste0("Sample_", 1:n_samples)
    
    return(expr_matrix)
  }
  
  # Load all datasets
  datasets <- list()
  
  # Psoriasis datasets
  datasets$psoriasis$GSE13355 <- load_dataset("GSE13355", 58, 64)
  datasets$psoriasis$GSE14905 <- load_dataset("GSE14905", 33, 21)
  
  # Crohn's disease datasets
  datasets$crohns$GSE75214 <- load_dataset("GSE75214", 67, 11)
  datasets$crohns$GSE102133 <- load_dataset("GSE102133", 65, 12)
  
  # Normalization function
  normalize_dataset <- function(expr_matrix) {
    cat("Normalizing dataset...\n")
    
    # Remove lowly expressed genes (more than 50% samples with 0 expression)
    keep_genes <- rowSums(expr_matrix > 0) >= 0.5 * ncol(expr_matrix)
    expr_matrix <- expr_matrix[keep_genes, ]
    cat("  Kept", sum(keep_genes), "genes after filtering\n")
    
    # Quantile normalization
    expr_normalized <- normalizeBetweenArrays(expr_matrix, method = "quantile")
    
    # Log2 transformation if needed
    if (max(expr_normalized) > 100) {
      expr_normalized <- log2(expr_normalized + 1)
      cat("  Applied log2 transformation\n")
    }
    
    return(expr_normalized)
  }
  
  # Apply normalization to all datasets
  cat("Applying normalization...\n")
  for(disease in names(datasets)) {
    for(dataset in names(datasets[[disease]])) {
      datasets[[disease]][[dataset]] <- normalize_dataset(datasets[[disease]][[dataset]])
    }
  }
  
  # Quality control - create box plots
  create_qc_plots <- function(datasets) {
    cat("Creating quality control plots...\n")
    
    for(disease in names(datasets)) {
      for(dataset in names(datasets[[disease]])) {
        expr_matrix <- datasets[[disease]][[dataset]]
        
        # Create box plot
        p <- ggplot(data = reshape2::melt(expr_matrix), 
                   aes(x = Var2, y = value)) +
          geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
          theme_minimal() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) +
          labs(title = paste("QC Plot:", dataset),
               x = "Samples",
               y = "Expression Level")
        
        ggsave(paste0("figures/", dataset, "_qc.png"), p, width = 10, height = 6)
      }
    }
  }
  
  # Create QC plots
  if(!dir.exists("figures")) dir.create("figures")
  create_qc_plots(datasets)
  
  # Save processed data
  saveRDS(datasets, "data/processed/normalized_datasets.rds")
  
  cat("Data preprocessing completed!\n")
  cat("Summary:\n")
  cat("- Processed 4 datasets (2 psoriasis, 2 Crohn's)\n")
  cat("- Applied quantile normalization\n")
  cat("- Performed quality control\n")
  cat("- Saved normalized data to data/processed/normalized_datasets.rds\n")
  
  return(list(
    datasets = datasets,
    sample_info = sample_info
  ))
}
