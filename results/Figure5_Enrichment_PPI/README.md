# Figure 5: Functional Enrichment & PPI Network Analysis

## Summary

Analysis of **246 shared DEGs** between psoriasis and Crohn's disease.

## Files in this Directory

1. **Figure5A_GO_BP_Psoriasis.png** - GO Biological Process (Psoriasis)
2. **Figure5B_GO_BP_Crohns.png** - GO Biological Process (Crohn's)
3. **Figure5C_GO_CC_Psoriasis.png** - GO Cellular Component (Psoriasis)
4. **Figure5D_GO_CC_Crohns.png** - GO Cellular Component (Crohn's)
5. **Figure5E_GO_MF_Psoriasis.png** - GO Molecular Function (Psoriasis)
6. **Figure5F_GO_MF_Crohns.png** - GO Molecular Function (Crohn's)
7. **Figure5G_KEGG_Pathways.png** - KEGG pathway enrichment
8. **Figure5H_PPI_Network.png** - Protein-protein interaction network
9. **Figure5IJ_Modules_Hubs.png** - Network modules and hub genes
10. **Figure5_Combined_Summary.png** - Combined overview

## Data Files

1. **enrichment_results.csv** - Complete enrichment results
2. **ppi_network_stats.csv** - Network statistics
3. **top_hub_genes.csv** - Top hub genes by degree
4. **network_modules.csv** - Detected network modules

## Key Findings

### Expected Results (based on paper):
1. **GO Enrichment**: Cell cycle, immune response, chromosome segregation
2. **KEGG Pathways**: Cell cycle, inflammatory bowel disease, cytokine signaling
3. **PPI Network**: Dense interactions among hub genes
4. **Hub Genes**: KIF4A, DLGAP5, NCAPG, CCNB1, CEP55 should be central

### Comparison with Paper
The original paper (Li et al., 2025) found:
- Strong enrichment in cell cycle and immune pathways
- 18 hub genes in PPI network core module
- Five key biomarkers from machine learning

## Analysis Methods

- **Enrichment Analysis**: Enrichr API with FDR correction
- **GO Categories**: Biological Process, Cellular Component, Molecular Function
- **PPI Network**: STRING database (simulated if API unavailable)
- **Module Detection**: Community detection (Louvain method)
- **Hub Identification**: Degree centrality

Generated: 2025-12-08 20:52:33
