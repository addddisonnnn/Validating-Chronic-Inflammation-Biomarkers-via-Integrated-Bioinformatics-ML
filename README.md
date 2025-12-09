# Validation of Chronic Inflammation Biomarkers via Integrated ML & Bioinformatics

This project replicates the analysis from the paper "Identification and validation of shared biomarkers and drug repurposing in psoriasis and Crohn’s disease: integrating bioinformatics, machine learning, and experimental approaches" (Li et al., 2025). The study investigates shared molecular mechanisms between two chronic inflammatory diseases to identify common biomarkers and potential drug repurposing candidates.

## Problem Statement
Psoriasis and Crohn's Disease are chronic inflammatory conditions with significant similarities:
- Shared genetic risk factors
- Frequent co-occurrence in patients  
- Overlapping immune dysregulation
- Lack of shared diagnostic biomarkers and treatments

Instead of studying each disease separately, this research aims to find common molecular mechanisms that could improve diagnosis and identify existing drugs to treat both conditions.

## Key Findings (Paper)
Li et al. (2025) identified five shared overactive genes in both diseases: **KIF4A, DLGAP5, NCAPG, CCNB1, and CEP55**. These genes control cell division and are linked to immune system dysregulation in these conditions.

## Analysis Pipeline (Paper Methodology)

### Phase 1: Data Acquisition & Preprocessing **COMPLETE**
- **What**: Download and prepare gene expression datasets
- **Why**: Need sick vs. healthy comparison to find disease-related genes
- **Datasets**: 
  - **Primary Analysis**: GSE13355 (psoriasis), GSE75214 (Crohn's)
  - **Validation**: GSE14905 (psoriasis), GSE102133 (Crohn's)
  - **Single-cell**: GSE162183 (psoriasis skin), GSE214695 (Crohn's colon)
- **Preprocessing**: Probe-to-gene conversion, background correction, quantile normalization, log2 transformation

### Phase 2: Differential Expression Analysis **COMPLETE**
- **What**: Identify significantly different genes between disease and control groups
- **Why**: Find potential biomarkers
- **Method**: Linear models for microarray data (limma)
- **Visualization**: PCA plots, volcano plots, heatmaps

### Phase 3: Gene Set Enrichment Analysis
- **What**: Identify overrepresented biological pathways
- **Database**: MSigDB Hallmark gene sets
- **Why**: Understand biological context of DEGs and pathway similarities

### Phase 4: Weighted Gene Co-expression Network Analysis (WGCNA)
- **What**: Find functional gene modules
- **Why**: Discover gene networks and identify modules associated with disease traits
- **Found**: 79 key genes important in both diseases

### Phase 5: Functional Enrichment Analysis
- **What**: Understand biological functions of shared genes
- **Why**: Identify key biological processes
- **Method**: GO Analysis, KEGG Pathways with clusterProfiler
- **Found**: Genes mainly control cell division and immune responses

### Phase 6: Protein-Protein Interaction Networks
- **What**: Map protein interactions of shared genes
- **Why**: Identify central hub genes and prioritize biologically important genes
- **Method**: STRING Database, MCODE for hub detection
- **Found**: 18 core hub genes from the 79

### Phase 7: Machine Learning Biomarker Selection
- **What**: Use AI algorithms to select best diagnostic biomarkers
- **Why**: Objectively select most predictive genes and optimize diagnostic accuracy
- **Method**: RF, SVM, XGBoost, Decision Tree, LASSO, KNN
- **Found**: 5 final biomarkers from the 18 hub genes

### Phase 8: Diagnostic Validation
- **What**: Test if these 5 genes can accurately diagnose both diseases
- **Why**: Confirm clinical utility across independent datasets
- **Method**: ROC curves with AUC calculation
- **Result**: Excellent diagnostic accuracy (AUC > 0.9)

### Phase 9: Immune Infiltration Analysis
- **What**: Analyze immune cell composition and relationship to hub genes
- **Why**: Link genes to specific immune processes and explain inflammation
- **Method**: Single-sample GSEA, correlation analysis

### Phase 10: Single-Cell RNA Sequencing Analysis
- **What**: Identify cell types expressing these genes
- **Why**: Cellular resolution of gene expression and validate bulk findings
- **Method**: Seurat package, UMAP clustering
- **Found**: Mainly in skin cells (psoriasis) and gut lining cells (Crohn's)

### Phase 11: Drug Repurposing
- **What**: Search for existing drugs targeting the 5 hub genes
- **Why**: Leverage existing safety data for cost-effective drug discovery
- **Method**: DSigDB through Enrichr
- **Found**: Etoposide, Lucanthone, Piroxicam as top candidates

### Phase 12: Molecular Docking
- **What**: Simulate drug binding to target proteins
- **Why**: Predict binding affinity and support drug selection
- **Found**: Etoposide binds strongest

### Phase 13: Experimental Validation
- **What**: Test drugs on human cells in the lab
- **Why**: Confirm computational predictions and validate biological effects
- **Method**: CCK-8 viability, qRT-PCR gene expression
- **Found**: Etoposide reduced disease markers

## 📁 Project Structure
``` text
Validating-Chronic-Inflammation-Biomarkers-via-Integrated-Bioinformatics-ML-Paper/
├── 📂 data/ # All data files
│ ├── 📂 raw/ # Raw microarray data from GEO (4 datasets)
│ └── 📂 processed/ # Processed expression matrices & DEG results
│
├── 📂 scripts/ # Analysis scripts
│ ├── 📂 R/ # R analysis pipelines (4 complete scripts)
│ └── 📂 python/ # Python utilities for metadata creation
│
├── 📂 metadata/ # Sample annotations and group assignments
│ ├── GSE13355_sample_groups.csv # 58 Psoriasis | 64 Control
│ ├── GSE75214_sample_groups_fixed.csv # 67 CD | 11 Control (ileal only)
│ ├── GSE14905_sample_groups.csv # 33 Psoriasis | 21 Control
│ └── GSE102133_sample_groups.csv # 65 CD | 12 Control
│
├── 📂 results/ # Analysis outputs
│ └── 📂 plots/ # Visualization files (4 volcano plots)
│
├── Paper.pdf # Original research paper
└── README.md # This documentation
```


## Results Summary

### **All 4 Datasets Successfully Processed**

| Dataset | Type | Samples | Disease | Control | Genes Tested | Significant DEGs | Upregulated | Downregulated |
|---------|------|---------|---------|---------|--------------|------------------|-------------|---------------|
| **GSE13355** | Psoriasis Training | 122 | 58 | 64 | 21,358 | **1,728** | 929 | 799 |
| **GSE75214** | Crohn's Training | 78 | 67 | 11 | 19,998 | **1,006** | 562 | 444 |
| **GSE14905** | Psoriasis Validation | 54 | 33 | 21 | 21,358 | **2,431** | 1,290 | 1,141 |
| **GSE102133** | Crohn's Validation | 77 | 65 | 12 | 19,998 | **1,002** | 587 | 415 |
| **Total** | **All Datasets** | **331** | **223** | **108** | **-** | **6,167** | **3,368** | **2,799** |

### 🔍 **Key Observations:**
1. **Psoriasis shows more DEGs** than Crohn's disease (1,728 vs 1,006 in training sets)
2. **Validation sets show higher DEG counts**, possibly due to stricter sample selection
3. **Consistent up/down regulation patterns** across training and validation sets
4. **Good data quality**: All datasets show clear separation in PCA/volcano plots

### 📈 **Dataset Statistics:**
- **Total samples processed**: 331
- **Total disease samples**: 223 (67.4%)
- **Total control samples**: 108 (32.6%)
- **Total significant DEGs identified**: 6,167
- **Average DEGs per dataset**: 1,542

## Current Progress

### ✅ **Completed** (Phases 1-2)
- ✅ Downloaded and organized 4 GEO microarray datasets (331 total samples)
- ✅ Created comprehensive metadata for all samples  
- ✅ Implemented reproducible preprocessing pipeline
- ✅ Performed differential expression analysis for all datasets
- ✅ Generated normalized expression matrices and DEG results
- ✅ Created quality control visualizations (volcano plots)

### **Available Results Files:**
```text
data/processed/
├── GSE13355_normalized_gene_expression.csv # 21,358 × 122 matrix
├── GSE13355_DEGs.csv # 1,728 significant genes
├── GSE75214_normalized_gene_expression.csv # 19,998 × 78 matrix
├── GSE75214_DEGs.csv # 1,006 significant genes
├── GSE14905_normalized_gene_expression.csv # 21,358 × 54 matrix
├── GSE14905_DEGs.csv # 2,431 significant genes
├── GSE102133_normalized_gene_expression.csv # 19,998 × 77 matrix
└── GSE102133_DEGs.csv # 1,002 significant genes
```

## Next Steps
1. **Identify Shared DEGs** - Find overlapping genes between psoriasis and Crohn's
2. **Perform WGCNA** - Identify co-expression modules and hub genes  
3. **Functional Enrichment** - GO/KEGG pathway analysis of shared genes
4. **Machine Learning** - Apply ML algorithms to identify diagnostic biomarkers
5. **Validation** - Test biomarkers on independent datasets

## Tools & Technologies
- **R**: `limma`, `oligo`, `WGCNA`, `clusterProfiler`, `ggplot2`
- **Python**: `pandas` for metadata processing
- **Bioinformatics**: GEO database, Affymetrix microarrays
- **Analysis Platforms**: SCC HPC cluster

## Citation
Li X, Cao H, Niu M, Liu Q, Liang B, Hou J, Tu J, Gao J. Identification and validation of shared biomarkers and drug repurposing in psoriasis and Crohn's disease: integrating bioinformatics, machine learning, and experimental approaches. Front Immunol. 2025 May 8;16:1587705. doi: 10.3389/fimmu.2025.1587705. PMID: 40406126; PMCID: PMC12095375.

---

**Project Status**: Phase 1-2 Complete | **Next**: Shared DEG Analysis  
**Last Updated**: Dec 8, 2025 | **Total Samples Processed**: 331 | **Total DEGs**: 6,167