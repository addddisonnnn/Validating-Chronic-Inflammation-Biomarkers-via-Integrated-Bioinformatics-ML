# complete_deg_analysis_GSE13355.R
# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install oligo and dependencies
BiocManager::install("oligo")

# Also install annotation package for your array platform
# For GSE13355 (hgu133plus2):
BiocManager::install("hgu133plus2.db")

BiocManager::install(c("ggplot2", "dplyr"))

# complete_analysis_GSE13355_final.R
library(oligo)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(dplyr)

# Set paths
cel_path <- "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/data/raw/GSE13355"
output_dir <- "/projectnb/bf528/students/addisony/project1/Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/results"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Load group assignments from our metadata analysis
group_info <- read.csv("GSE13355_sample_groups.csv")
analysis_samples <- group_info %>%
  filter(!is.na(Group)) %>%
  filter(Group %in% c("Psoriasis", "Control"))

cat("=== SAMPLE DISTRIBUTION ===\n")
print(table(analysis_samples$Group))
cat("Total samples for analysis:", nrow(analysis_samples), "\n")

# 2. Get CEL files
cat("\nFinding CEL files...\n")
cel_files <- list.files(cel_path, pattern = "\\.CEL(\\.gz)?$", 
                        full.names = TRUE, ignore.case = TRUE)

# Extract GSM IDs from filenames
file_gsm <- gsub("\\.CEL(\\.gz)?$", "", basename(cel_files))
names(cel_files) <- file_gsm

# Match with our analysis samples
matched_gsm <- intersect(analysis_samples$GSM_ID, file_gsm)
cat("Matched samples:", length(matched_gsm), "/", nrow(analysis_samples), "\n")

if (length(matched_gsm) < nrow(analysis_samples)) {
  missing <- setdiff(analysis_samples$GSM_ID, file_gsm)
  cat("Missing CEL files for:", length(missing), "samples\n")
  print(missing[1:5])
}

# 3. Subset to matched samples
selected_samples <- analysis_samples %>%
  filter(GSM_ID %in% matched_gsm)

selected_cel_files <- cel_files[selected_samples$GSM_ID]

# Ensure order matches
selected_samples <- selected_samples[match(names(selected_cel_files), selected_samples$GSM_ID), ]

# 4. Read and normalize
cat("\nReading and normalizing data...\n")
raw_data <- read.celfiles(selected_cel_files)
norm_data <- rma(raw_data)
expr_mat <- exprs(norm_data)

cat("Normalized expression matrix dimensions:", dim(expr_mat), "\n")

# 5. Map probes to gene symbols
cat("Mapping probes to gene symbols...\n")
probe_ids <- rownames(expr_mat)
gene_symbols <- mapIds(hgu133plus2.db, 
                       keys = probe_ids, 
                       column = "SYMBOL", 
                       keytype = "PROBEID")

# Remove probes without gene symbol
keep <- !is.na(gene_symbols)
expr_mat <- expr_mat[keep, ]
gene_symbols <- gene_symbols[keep]

cat("After filtering NA gene symbols:", dim(expr_mat), "\n")

# 6. Aggregate by gene symbol
cat("Aggregating by gene symbol...\n")
unique_genes <- unique(gene_symbols)
gene_mat <- matrix(NA, nrow = length(unique_genes), ncol = ncol(expr_mat))
rownames(gene_mat) <- unique_genes
colnames(gene_mat) <- colnames(expr_mat)

pb <- txtProgressBar(min = 0, max = length(unique_genes), style = 3)
for (i in 1:length(unique_genes)) {
  gene <- unique_genes[i]
  idx <- which(gene_symbols == gene)
  if (length(idx) > 1) {
    gene_mat[gene, ] <- colMeans(expr_mat[idx, ])
  } else {
    gene_mat[gene, ] <- expr_mat[idx, ]
  }
  if (i %% 1000 == 0) setTxtProgressBar(pb, i)
}
close(pb)

# 7. Save normalized data
output_file <- file.path(output_dir, "GSE13355_normalized_gene_expression.csv")
write.csv(gene_mat, output_file)
cat("\nSaved normalized expression matrix\n")

# 8. Differential expression
cat("\nPerforming differential expression analysis...\n")

# Create group vector (ensure correct order)
group_vector <- factor(selected_samples$Group)
cat("Final group distribution:\n")
print(table(group_vector))

design <- model.matrix(~0 + group_vector)
colnames(design) <- levels(group_vector)

fit <- lmFit(gene_mat, design)
contrast_matrix <- makeContrasts(Psoriasis - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 9. Extract results
deg_table <- topTable(fit2, number = Inf, adjust.method = "BH")
deg_table$Gene <- rownames(deg_table)

# Apply paper's filters (|logFC| > 0.585, adj.P.Val < 0.05)
deg_sig <- subset(deg_table, abs(logFC) > 0.585 & adj.P.Val < 0.05)
deg_sig <- deg_sig[order(deg_sig$adj.P.Val), ]

cat("\n=== DEG RESULTS ===\n")
cat("Total genes tested:", nrow(deg_table), "\n")
cat("Significant DEGs (|logFC| > 0.585, adj.P.Val < 0.05):", nrow(deg_sig), "\n")
cat("Upregulated in Psoriasis:", sum(deg_sig$logFC > 0), "\n")
cat("Downregulated in Psoriasis:", sum(deg_sig$logFC < 0), "\n")

# 10. Save DEG results
deg_file <- file.path(output_dir, "GSE13355_DEGs.csv")
write.csv(deg_sig, deg_file, row.names = FALSE)
cat("\nSaved DEG results to:", deg_file, "\n")

# 11. Create volcano plot
cat("\nCreating volcano plot...\n")
deg_table$Significant <- deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 0.585
deg_table$Direction <- ifelse(deg_table$logFC > 0, "Up", "Down")
deg_table$Direction[!deg_table$Significant] <- "Not Sig"

# Get top 10 DEGs for labeling
top_up <- deg_table %>% 
  filter(Significant, logFC > 0) %>% 
  arrange(adj.P.Val) %>% 
  head(5)

top_down <- deg_table %>% 
  filter(Significant, logFC < 0) %>% 
  arrange(adj.P.Val) %>% 
  head(5)

top_genes <- rbind(top_up, top_down)

volcano_plot <- ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Sig" = "gray", "Up" = "red", "Down" = "blue")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text(data = top_genes, aes(label = Gene), vjust = -0.5, size = 3, color = "black") +
  theme_minimal() +
  labs(title = "GSE13355: Psoriasis vs Control",
       subtitle = paste(nrow(deg_sig), "significant DEGs (|logFC| > 0.585, adj.P.Val < 0.05)"),
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)") +
  theme(legend.position = "bottom")

plot_file <- file.path(output_dir, "GSE13355_volcano_plot.png")
ggsave(plot_file, volcano_plot, width = 10, height = 8, dpi = 300)
cat("Saved volcano plot to:", plot_file, "\n")

# 12. Save full DEG table for reference
full_deg_file <- file.path(output_dir, "GSE13355_full_DEG_table.csv")
write.csv(deg_table, full_deg_file, row.names = FALSE)

# 13. Summary statistics
summary_stats <- data.frame(
  Metric = c("Total samples", "Psoriasis samples", "Control samples", 
             "Total genes", "Significant DEGs", "Upregulated", "Downregulated"),
  Count = c(nrow(selected_samples), 
            sum(group_vector == "Psoriasis"), 
            sum(group_vector == "Control"),
            nrow(deg_table), nrow(deg_sig),
            sum(deg_sig$logFC > 0), sum(deg_sig$logFC < 0))
)

summary_file <- file.path(output_dir, "GSE13355_analysis_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved in:", output_dir, "\n")
cat("\nFiles created:\n")
cat("1. Normalized expression matrix\n")
cat("2. Significant DEGs\n")
cat("3. Full DEG table\n")
cat("4. Volcano plot\n")
cat("5. Analysis summary\n")