# scripts/02_differential_expression.R
# ==============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS
# Identify significantly different genes between disease and control groups
# ==============================================================================

perform_differential_expression <- function(processed_data) {
  cat("Performing differential expression analysis...\n")
  
  # Function to perform limma analysis
  perform_limma_analysis <- function(expr_matrix, n_patients, n_controls, dataset_name) {
    cat("  Analyzing", dataset_name, "...\n")
    
    # Create design matrix
    group <- factor(c(rep("Patient", n_patients), rep("Control", n_controls)))
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
    # Fit linear model
    fit <- lmFit(expr_matrix, design)
    
    # Create contrast (Patient vs Control)
    contrast.matrix <- makeContrasts(Patient - Control, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Get results
    results <- topTable(fit2, number = Inf, adjust.method = "BH")
    
    return(results)
  }
  
  # Analyze psoriasis datasets
  cat("Analyzing psoriasis datasets...\n")
  deg_psoriasis <- perform_limma_analysis(
    processed_data$datasets$psoriasis$GSE13355,
    58, 64, "GSE13355"
  )
  
  # Analyze Crohn's datasets
  cat("Analyzing Crohn's disease datasets...\n")
  deg_crohns <- perform_limma_analysis(
    processed_data$datasets$crohns$GSE75214,
    67, 11, "GSE75214"
  )
  
  # Filter significant DEGs (adj.p < 0.05, |logFC| > 0.585)
  filter_significant_genes <- function(deg_results, dataset_name) {
    significant <- deg_results[deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 0.585, ]
    cat("  ", dataset_name, ":", nrow(significant), "significant DEGs\n")
    return(significant)
  }
  
  sig_psoriasis <- filter_significant_genes(deg_psoriasis, "Psoriasis")
  sig_crohns <- filter_significant_genes(deg_crohns, "Crohn's")
  
  # Find shared DEGs
  shared_genes <- intersect(rownames(sig_psoriasis), rownames(sig_crohns))
  cat("Shared DEGs between psoriasis and Crohn's:", length(shared_genes), "\n")
  
  # Create visualization functions
  create_volcano_plot <- function(deg_results, title) {
    deg_results$significant <- ifelse(
      deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 0.585,
      "Significant", "Not significant"
    )
    
    p <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("gray", "red")) +
      theme_minimal() +
      theme(legend.position = "top") +
      labs(title = title,
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           color = "Significance")
    
    return(p)
  }
  
  create_heatmap <- function(expr_matrix, deg_genes, n_patients, n_controls, title) {
    # Subset to significant genes
    if(length(deg_genes) > 0) {
      deg_expr <- expr_matrix[deg_genes, ]
      
      # Create annotation
      annotation <- data.frame(
        Group = c(rep("Patient", n_patients), rep("Control", n_controls))
      )
      rownames(annotation) <- colnames(deg_expr)
      
      # Create heatmap
      pheatmap(deg_expr,
               scale = "row",
               annotation_col = annotation,
               show_rownames = FALSE,
               show_colnames = FALSE,
               main = title,
               color = colorRampPalette(c("blue", "white", "red"))(50))
    }
  }
  
  # Create visualizations
  cat("Creating visualizations...\n")
  
  # Volcano plots
  p1 <- create_volcano_plot(deg_psoriasis, "Psoriasis - GSE13355")
  p2 <- create_volcano_plot(deg_crohns, "Crohn's Disease - GSE75214")
  
  # Save plots
  ggsave("figures/volcano_psoriasis.png", p1, width = 8, height = 6)
  ggsave("figures/volcano_crohns.png", p2, width = 8, height = 6)
  
  # Save results
  results <- list(
    psoriasis_degs = sig_psoriasis,
    crohns_degs = sig_crohns,
    shared_genes = shared_genes,
    full_results = list(
      psoriasis = deg_psoriasis,
      crohns = deg_crohns
    )
  )
  
  saveRDS(results, "results/differential_expression_results.rds")
  
  # Write summary table
  summary_table <- data.frame(
    Dataset = c("Psoriasis (GSE13355)", "Crohn's (GSE75214)", "Shared"),
    DEGs_Upregulated = c(
      sum(sig_psoriasis$logFC > 0),
      sum(sig_crohns$logFC > 0),
      NA
    ),
    DEGs_Downregulated = c(
      sum(sig_psoriasis$logFC < 0),
      sum(sig_crohns$logFC < 0),
      NA
    ),
    Total_DEGs = c(
      nrow(sig_psoriasis),
      nrow(sig_crohns),
      length(shared_genes)
    )
  )
  
  write.csv(summary_table, "results/deg_summary.csv", row.names = FALSE)
  
  cat("Differential expression analysis completed!\n")
  cat("Summary:\n")
  cat("- Psoriasis:", nrow(sig_psoriasis), "DEGs\n")
  cat("- Crohn's:", nrow(sig_crohns), "DEGs\n")
  cat("- Shared:", length(shared_genes), "DEGs\n")
  
  return(results)
}
