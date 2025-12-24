# Save as: scripts/R/wgcna_analysis.R

# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("WGCNA", "impute", "preprocessCore"))
install.packages(c("reshape2", "ggplot2", "dplyr", "corrplot"))

library(WGCNA)
library(ggplot2)
library(dplyr)
library(reshape2)
library(corrplot)

# Enable multithreading (optional)
enableWGCNAThreads()

# Set paths
base_dir <- "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper"
output_dir <- file.path(base_dir, "results/wgcna")

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load normalized expression data
cat("Loading expression data...\n")
psoriasis_expr <- read.csv(file.path(base_dir, "data/processed/GSE13355_normalized_gene_expression.csv"), 
                           row.names = 1, check.names = FALSE)
crohns_expr <- read.csv(file.path(base_dir, "data/processed/GSE75214_normalized_gene_expression.csv"), 
                        row.names = 1, check.names = FALSE)

cat("Psoriasis expression matrix:", dim(psoriasis_expr), "\n")
cat("Crohn's expression matrix:", dim(crohns_expr), "\n")

# Load sample metadata
psoriasis_meta <- read.csv(file.path(base_dir, "metadata/GSE13355_sample_groups.csv"))
crohns_meta <- read.csv(file.path(base_dir, "metadata/GSE75214_sample_groups_fixed.csv"))

# Prepare trait data (binary: disease = 1, control = 0)
prepare_traits <- function(expr_matrix, metadata) {
  samples <- colnames(expr_matrix)
  traits <- data.frame(row.names = samples)
  
  for (sample in samples) {
    # Extract GSM ID from column name (might have sample.1 suffix)
    gsm_id <- sub("\\..*$", "", sample)
    
    # Find group in metadata
    group <- metadata$Group[metadata$GSM_ID == gsm_id]
    
    if (length(group) > 0) {
      traits[sample, "Disease"] <- ifelse(group == "Psoriasis" | group == "CD", 1, 0)
    } else {
      traits[sample, "Disease"] <- NA
    }
  }
  
  return(traits)
}

psoriasis_traits <- prepare_traits(psoriasis_expr, psoriasis_meta)
crohns_traits <- prepare_traits(crohns_expr, crohns_meta)

cat("\nTrait data prepared:\n")
cat("Psoriasis - Disease samples:", sum(psoriasis_traits$Disease == 1), 
    "Control samples:", sum(psoriasis_traits$Disease == 0), "\n")
cat("Crohn's - Disease samples:", sum(crohns_traits$Disease == 1), 
    "Control samples:", sum(crohns_traits$Disease == 0), "\n")

# Function to perform WGCNA
perform_wgcna <- function(expr_data, traits, dataset_name) {
  cat("Performing WGCNA for", dataset_name, "\n")

  # Transpose expression data (genes as columns, samples as rows)
  datExpr <- t(expr_data)
  
  # Check for missing values
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  
  if (!gsg$allOK) {
    # Remove genes and samples with too many missing values
    if (sum(!gsg$goodGenes) > 0)
      datExpr <- datExpr[, gsg$goodGenes]
    if (sum(!gsg$goodSamples) > 0)
      datExpr <- datExpr[gsg$goodSamples, ]
    cat("Removed", sum(!gsg$goodGenes), "genes and", sum(!gsg$goodSamples), "samples with missing data\n")
  }
  
  # Match samples between expression and trait data
  sample_names <- rownames(datExpr)
  trait_samples <- rownames(traits)
  common_samples <- intersect(sample_names, trait_samples)
  
  datExpr <- datExpr[common_samples, ]
  datTraits <- traits[common_samples, , drop = FALSE]
  
  cat("Final data dimensions:", dim(datExpr), "\n")
  cat("Matched samples:", length(common_samples), "\n")
  
  # 1. Choose soft threshold power
  cat("\n1. Choosing soft threshold power...\n")
  powers <- c(1:20)
  
  # Call network topology analysis function
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
                           networkType = "signed", corFnc = "bicor")
  
  # Plot scale independence and mean connectivity
  png(file.path(output_dir, paste0(dataset_name, "_soft_threshold.png")), 
      width = 10, height = 5, units = "in", res = 300)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index vs. soft threshold power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.80, col = "red")  # Red line at R² = 0.80
  
  # Mean connectivity vs. soft threshold power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  dev.off()
  
  # Choose power based on scale-free topology criterion
  power <- sft$powerEstimate
  if (is.na(power)) {
    power <- 6  # Default if automatic selection fails
    cat("Automatic power selection failed, using power =", power, "\n")
  } else {
    cat("Selected soft threshold power:", power, "\n")
  }
  
  # 2. Construct network and detect modules
  cat("\n2. Constructing network and detecting modules...\n")
  
  # Calculate adjacency
  adjacency <- adjacency(datExpr, power = power, type = "signed", corFnc = "bicor")
  
  # Transform adjacency into TOM
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM <- 1 - TOM
  
  # Cluster genes using TOM-based dissimilarity
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  # Module identification using dynamic tree cut
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
  
  # Convert numeric labels to colors
  dynamicColors <- labels2colors(dynamicMods)
  
  # Calculate module eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  
  # Merge similar modules
  MEDissThres <- 0.25  # Merge modules with correlation > 0.75
  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs
  
  # Plot dendrogram with module colors
  png(file.path(output_dir, paste0(dataset_name, "_dendrogram.png")), 
      width = 12, height = 8, units = "in", res = 300)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste(dataset_name, "gene dendrogram and module colors"))
  dev.off()
  
  # 3. Relate modules to traits
  cat("\n3. Relating modules to traits...\n")
  
  # Define numbers of genes and samples
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  # Recalculate MEs with color labels
  MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes
  MEs <- orderMEs(MEs0)
  
  # Calculate module-trait relationships
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create correlation heatmap
  png(file.path(output_dir, paste0(dataset_name, "_module_trait_heatmap.png")), 
      width = 8, height = 10, units = "in", res = 300)
  
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.7,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships in", dataset_name))
  dev.off()
  
  # 4. Identify significant modules
  cat("\n4. Identifying significant modules...\n")
  
  # Find modules with significant correlation to disease
  significant_modules <- which(moduleTraitPvalue[, "Disease"] < 0.05)
  
  if (length(significant_modules) > 0) {
    cat("Significant modules (p < 0.05):\n")
    for (i in significant_modules) {
      module_name <- rownames(moduleTraitCor)[i]
      correlation <- moduleTraitCor[i, "Disease"]
      p_value <- moduleTraitPvalue[i, "Disease"]
      cat(sprintf("  %s: r = %.3f, p = %.3e\n", module_name, correlation, p_value))
    }
  } else {
    cat("No modules significantly correlated with disease (p < 0.05)\n")
  }
  
  # 5. Extract genes from significant modules
  cat("\n5. Extracting genes from significant modules...\n")
  
  module_genes <- list()
  if (length(significant_modules) > 0) {
    for (i in significant_modules) {
      module_color <- gsub("ME", "", rownames(moduleTraitCor)[i])
      module_genes[[module_color]] <- colnames(datExpr)[mergedColors == module_color]
    }
  }
  
  # Save results
  results <- list(
    datExpr = datExpr,
    datTraits = datTraits,
    mergedColors = mergedColors,
    MEs = MEs,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    significant_modules = rownames(moduleTraitCor)[significant_modules],
    module_genes = module_genes,
    soft_power = power
  )
  
  # Save to RDS for later use
  saveRDS(results, file.path(output_dir, paste0(dataset_name, "_wgcna_results.rds")))
  
  # Save module assignments
  gene_module_df <- data.frame(
    Gene = colnames(datExpr),
    Module = mergedColors,
    stringsAsFactors = FALSE
  )
  write.csv(gene_module_df, 
            file.path(output_dir, paste0(dataset_name, "_gene_modules.csv")),
            row.names = FALSE)
  
  # Save module-trait correlations
  module_cor_df <- data.frame(
    Module = rownames(moduleTraitCor),
    Correlation = moduleTraitCor[, "Disease"],
    Pvalue = moduleTraitPvalue[, "Disease"]
  )
  write.csv(module_cor_df, 
            file.path(output_dir, paste0(dataset_name, "_module_correlations.csv")),
            row.names = FALSE)
  
  cat("\nWGCNA analysis for", dataset_name, "completed!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(results)
}

# Run WGCNA for both datasets
cat("STARTING WGCNA ANALYSIS FOR BOTH DATASETS\n")

psoriasis_results <- perform_wgcna(psoriasis_expr, psoriasis_traits, "psoriasis_GSE13355")
crohns_results <- perform_wgcna(crohns_expr, crohns_traits, "crohns_GSE75214")

# Compare results between datasets
cat("COMPARING RESULTS BETWEEN DATASETS\n")

# Load shared DEGs
shared_genes <- read.csv(file.path(base_dir, "results/shared_degs_training.csv"))$Gene

# Find which modules contain shared DEGs
analyze_shared_in_modules <- function(results, dataset_name, shared_genes) {
  gene_module_df <- data.frame(
    Gene = colnames(results$datExpr),
    Module = results$mergedColors
  )
  
  # Filter to shared genes
  shared_in_modules <- gene_module_df[gene_module_df$Gene %in% shared_genes, ]
  
  cat("\nShared DEGs in", dataset_name, "modules:\n")
  if (nrow(shared_in_modules) > 0) {
    module_counts <- table(shared_in_modules$Module)
    print(module_counts)
    
    # Save to file
    write.csv(shared_in_modules, 
              file.path(output_dir, paste0(dataset_name, "_shared_genes_in_modules.csv")),
              row.names = FALSE)
  } else {
    cat("No shared DEGs found in", dataset_name, "expression matrix\n")
  }
  
  return(shared_in_modules)
}

psoriasis_shared_modules <- analyze_shared_in_modules(psoriasis_results, "psoriasis", shared_genes)
crohns_shared_modules <- analyze_shared_in_modules(crohns_results, "crohns", shared_genes)

# Find overlapping modules between datasets
cat("FINDING OVERLAPPING MODULES\n")

# Get significant modules for each disease
psoriasis_sig_modules <- psoriasis_results$significant_modules
crohns_sig_modules <- crohns_results$significant_modules

cat("Psoriasis significant modules:", if(length(psoriasis_sig_modules) > 0) paste(psoriasis_sig_modules, collapse=", ") else "None", "\n")
cat("Crohn's significant modules:", if(length(crohns_sig_modules) > 0) paste(crohns_sig_modules, collapse=", ") else "None", "\n")

if (length(psoriasis_sig_modules) > 0 && length(crohns_sig_modules) > 0) {
  # Extract module colors (remove "ME" prefix)
  psoriasis_colors <- gsub("ME", "", psoriasis_sig_modules)
  crohns_colors <- gsub("ME", "", crohns_sig_modules)
  
  overlapping_colors <- intersect(psoriasis_colors, crohns_colors)
  
  if (length(overlapping_colors) > 0) {
    cat("\nOverlapping module colors:", paste(overlapping_colors, collapse=", "), "\n")
  } else {
    cat("\nNo overlapping module colors found\n")
  }
}

cat("WGCNA ANALYSIS COMPLETE!\n")
cat("Results saved in:", output_dir, "\n")
cat("\nNext steps:")
cat("\n1. Examine module-trait correlation heatmaps")
cat("\n2. Extract hub genes from significant modules")
cat("\n3. Perform functional enrichment on module genes")
cat("\n4. Compare with paper's findings (blue module for psoriasis, brown module for CD)\n")