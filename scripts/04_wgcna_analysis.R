# scripts/04_wgcna_analysis.R
# ==============================================================================
# STEP 4: WGCNA NETWORK ANALYSIS
# Identify co-expression modules associated with disease traits
# ==============================================================================

perform_wgcna_analysis <- function(processed_data, deg_results) {
  cat("Performing WGCNA analysis...\n")
  
  # Function to perform WGCNA on a single dataset
  perform_single_wgcna <- function(expr_matrix, trait_vector, dataset_name) {
    cat("  Running WGCNA for", dataset_name, "...\n")
    
    # Step 1: Select top 25% most variable genes
    gene_variance <- apply(expr_matrix, 1, var)
    variance_threshold <- quantile(gene_variance, 0.75)
    high_var_genes <- gene_variance > variance_threshold
    expr_filtered <- expr_matrix[high_var_genes, ]
    
    cat("    Selected", sum(high_var_genes), "highly variable genes\n")
    
    # Step 2: Choose soft threshold power
    powers <- c(1:20)
    sft <- pickSoftThreshold(t(expr_filtered), 
                            powerVector = powers, 
                            verbose = 5,
                            networkType = "signed")
    
    # Plot scale independence and mean connectivity
    png(paste0("figures/", dataset_name, "_soft_threshold.png"), 
        width = 10, height = 5)
    par(mfrow = c(1,2))
    
    # Scale independence plot
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = "Soft Threshold (power)", 
         ylab = "Scale Free Topology Model Fit, signed R^2",
         main = paste("Scale independence -", dataset_name))
    abline(h = 0.90, col = "red")
    
    # Mean connectivity plot
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold (power)", 
         ylab = "Mean Connectivity", 
         main = paste("Mean connectivity -", dataset_name))
    
    dev.off()
    
    # Use optimal power (or default to 12 if not found)
    optimal_power <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)
    cat("    Optimal soft threshold power:", optimal_power, "\n")
    
    # Step 3: Construct network
    net <- blockwiseModules(t(expr_filtered),
                           power = optimal_power,
                           TOMType = "signed",
                           minModuleSize = 30,
                           reassignThreshold = 0,
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = paste0("results/", dataset_name, "_TOM"),
                           verbose = 3)
    
    # Convert labels to colors
    module_colors <- labels2colors(net$colors)
    
    # Step 4: Calculate module-trait relationships
    MEs <- net$MEs
    moduleTraitCor <- cor(MEs, trait_vector, use = "p")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(expr_filtered))
    
    # Create module-trait relationship heatmap
    png(paste0("figures/", dataset_name, "_module_trait_heatmap.png"),
        width = 8, height = 6)
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    
    labeledHeatmap(Matrix = moduleTraitCor,
                  xLabels = "Disease",
                  yLabels = names(MEs),
                  ySymbols = names(MEs),
                  colorLabels = FALSE,
                  colors = blueWhiteRed(50),
                  textMatrix = textMatrix,
                  setStdMargins = FALSE,
                  cex.text = 0.7,
                  zlim = c(-1,1),
                  main = paste("Module-trait relationships -", dataset_name))
    dev.off()
    
    return(list(
      net = net,
      moduleTraitCor = moduleTraitCor,
      moduleTraitPvalue = moduleTraitPvalue,
      optimal_power = optimal_power,
      module_colors = module_colors
    ))
  }
  
  # Prepare trait vectors (1 = patient, 0 = control)
  psoriasis_trait <- as.numeric(c(rep(1, 58), rep(0, 64)))
  crohns_trait <- as.numeric(c(rep(1, 67), rep(0, 11)))
  
  # Run WGCNA for both diseases
  cat("Running WGCNA for psoriasis...\n")
  wgcna_psoriasis <- perform_single_wgcna(
    processed_data$datasets$psoriasis$GSE13355,
    psoriasis_trait,
    "psoriasis"
  )
  
  cat("Running WGCNA for Crohn's disease...\n")
  wgcna_crohns <- perform_single_wgcna(
    processed_data$datasets$crohns$GSE75214,
    crohns_trait,
    "crohns"
  )
  
  # Find key modules (highest absolute correlation with disease)
  find_key_module <- function(wgcna_result, dataset_name) {
    cor_values <- wgcna_result$moduleTraitCor[,1]
    key_module_index <- which.max(abs(cor_values))
    key_module <- names(cor_values)[key_module_index]
    correlation <- cor_values[key_module_index]
    
    cat("  ", dataset_name, "key module:", key_module, 
        "(correlation =", round(correlation, 3), ")\n")
    
    return(list(
      module_name = key_module,
      correlation = correlation,
      module_index = as.numeric(gsub("ME", "", key_module))
    ))
  }
  
  psoriasis_key <- find_key_module(wgcna_psoriasis, "Psoriasis")
  crohns_key <- find_key_module(wgcna_crohns, "Crohn's")
  
  # Extract genes from key modules
  extract_module_genes <- function(wgcna_result, module_index, gene_names) {
    module_genes <- gene_names[wgcna_result$net$colors == module_index]
    return(module_genes)
  }
  
  psoriasis_genes <- extract_module_genes(
    wgcna_psoriasis, 
    psoriasis_key$module_index,
    rownames(processed_data$datasets$psoriasis$GSE13355)
  )
  
  crohns_genes <- extract_module_genes(
    wgcna_crohns,
    crohns_key$module_index,
    rownames(processed_data$datasets$crohns$GSE75214)
  )
  
  cat("  Psoriasis module genes:", length(psoriasis_genes), "\n")
  cat("  Crohn's module genes:", length(crohns_genes), "\n")
  
  # Find shared genes between WGCNA modules and DEGs
  shared_genes <- intersect(psoriasis_genes, crohns_genes)
  shared_genes <- intersect(shared_genes, deg_results$shared_genes)
  
  cat("Shared genes (WGCNA modules ∩ DEGs):", length(shared_genes), "\n")
  
  # Save results
  results <- list(
    psoriasis = wgcna_psoriasis,
    crohns = wgcna_crohns,
    key_modules = list(
      psoriasis = psoriasis_key,
      crohns = crohns_key
    ),
    module_genes = list(
      psoriasis = psoriasis_genes,
      crohns = crohns_genes
    ),
    shared_genes = shared_genes
  )
  
  saveRDS(results, "results/wgcna_results.rds")
  
  # Create summary table
  summary_table <- data.frame(
    Dataset = c("Psoriasis", "Crohn's Disease"),
    Key_Module = c(psoriasis_key$module_name, crohns_key$module_name),
    Module_Correlation = c(psoriasis_key$correlation, crohns_key$correlation),
    Module_Genes = c(length(psoriasis_genes), length(crohns_genes)),
    Shared_Genes = c(length(shared_genes), length(shared_genes))
  )
  
  write.csv(summary_table, "results/wgcna_summary.csv", row.names = FALSE)
  
  cat("WGCNA analysis completed!\n")
  cat("Key findings:\n")
  cat("- Psoriasis key module:", psoriasis_key$module_name, 
      "(correlation:", round(psoriasis_key$correlation, 3), ")\n")
  cat("- Crohn's key module:", crohns_key$module_name,
      "(correlation:", round(crohns_key$correlation, 3), ")\n")
  cat("- Shared key genes:", length(shared_genes), "\n")
  
  return(results)
}
