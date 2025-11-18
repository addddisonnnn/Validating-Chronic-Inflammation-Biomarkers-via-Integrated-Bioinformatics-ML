# Validation-of-Chronic-Inflammation-Biomarkers-via-Integrated-ML-Bioinformatics-Paper
Documenting my progress towards validating and recreating a [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12095375/#s1) titled 'Identification and validation of shared biomarkers and drug repurposing in psoriasis and Crohn’s disease: integrating bioinformatics, machine learning, and experimental approaches'.

## Paper Summary

### Research Problem
Psoriases and Crohn's Disease are two chronic inflammatory conditions that have a myriad of similarities, including shared genetic risk factors, often occuring at the same time in paitents, have overlapping immune dysregulation, and lack shared diagnoistic biomarkers and treaments. Instead of studying each disease separately, scientists wanted to find common molecular mechanisms that could better diagnoses and find existing drugs to test and treat both conditions.

### Key Findings
Li et al. (2025) identified five shared overactive genes in both disease: KIF4A, DLGAP5, NCAPG, CCNB1, and CEP55. These genes are known for controlling cell division and have been linked to immune system problems in these condiitons.  

## Main Steps
1. Data procurement
    * What:
    * Why: 
2. Find Different Genes
What: Compared gene activity between patients and controls
Found: 223 genes that behave differently in BOTH diseases

3. Find Gene Teams
What: Used WGCNA to find groups of genes that work together
Found: 79 key genes that are important in both diseases

4. Understand What Genes Do
What: Analyzed the biological functions of these 79 genes
Found: They mainly control cell division and immune responses

5. Find Most Important Genes
What: Built interaction networks to find central "hub" genes
Found: 18 core hub genes from the 79

6. AI Selection
What: Used machine learning to pick the best diagnostic genes
Found: 5 final biomarkers from the 18 hub genes

7. Validation
What: Tested if these 5 genes can accurately diagnose both diseases
Result: Excellent diagnostic accuracy (AUC > 0.9)

8. Immune Analysis
What: Checked how these genes relate to immune cells
Found: Connected to T-cells and other immune cells in both diseases

9. Single-Cell Analysis
What: Looked at which specific cell types use these genes
Found: Mainly in skin cells (psoriasis) and gut lining cells (Crohn's)

10. Drug Discovery
What: Searched for existing drugs that target these 5 genes
Found: Etoposide, Lucanthone, Piroxicam

11. Computer Drug Testing
What: Simulated how well drugs bind to target proteins
Found: Etoposide binds strongest

12. Lab Experiments 
What: Tested drugs on human cells in the lab
Found: Etoposide reduced disease markers

### Citation
Li X, Cao H, Niu M, Liu Q, Liang B, Hou J, Tu J, Gao J. Identification and validation of shared biomarkers and drug repurposing in psoriasis and Crohn's disease: integrating bioinformatics, machine learning, and experimental approaches. Front Immunol. 2025 May 8;16:1587705. doi: 10.3389/fimmu.2025.1587705. PMID: 40406126; PMCID: PMC12095375.
