# Figure 2: Differential Expression Analysis

## Summary

### Dataset Information
- **Psoriasis (GSE13355)**: 58 disease vs 64 control samples
- **Crohn's Disease (GSE75214)**: 67 disease vs 11 control samples

### Psoriasis DEGs: 1728
### Crohn's Disease DEGs: 1006
### Shared DEGs: 246
- Psoriasis unique: 1482
- Crohn's unique: 760

## Files in this Directory

1. **Figure2_AB_PCA.png** - Principal Component Analysis plots
   - A) Psoriasis (GSE13355)
   - B) Crohn's Disease (GSE75214)

2. **Figure2_CD_Volcano.png** - Volcano plots of differentially expressed genes
   - C) Psoriasis DEGs
   - D) Crohn's Disease DEGs

3. **Figure2E_Venn.png** - Venn diagram of shared DEGs
   - E) Overlap between psoriasis and Crohn's disease

4. **Figure2_FG_Heatmaps.png** - Expression heatmaps of shared DEGs
   - F) Psoriasis samples
   - G) Crohn's Disease samples

5. **Figure2_Combined_Summary.png** - Combined overview figure

6. **shared_deg_genes.csv** - List of shared differentially expressed genes

## Analysis Parameters

- **Differential expression criteria**: |log2FC| > 0.585 and adjusted p-value < 0.05
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Heatmap normalization**: Z-score by gene
- **PCA scaling**: Standard scaling (mean=0, variance=1)

## Comparison with Original Paper

The original paper (Li et al., 2025) reported:
- 223 shared DEGs between psoriasis and Crohn's disease
- Strong enrichment in cell cycle and immune response pathways
- Five key hub genes: KIF4A, DLGAP5, NCAPG, CCNB1, CEP55

Generated: 2025-12-08 20:37:16
