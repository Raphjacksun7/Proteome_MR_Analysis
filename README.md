# Proteome-wide MR Analysis & Colocalization Pipeline

This repository contains the R scripts used for the integrative proteogenomic analysis reported in our paper:

**Proteome-wide Mendelian randomization identifies circulating proteins causally associated with childhood obesity**

The analysis pipeline includes standard Mendelian Randomization (MR), colocalization, two‐step MR mediation, and reverse MR. Together, these scripts allow the investigation of the relationships between circulating protein levels, BMI, and childhood/adult obesity traits.

## Overview of Scripts

1. **MR_analysis.R**  
   Performs standard MR analysis by reading in exposure (protein) and outcome (e.g., BMI) GWAS data, harmonizing the datasets with the TwoSampleMR package, and compiling the MR estimates for each protein group.

2. **Colocalization_analysis.R**  
   Merges and processes exposure and outcome GWAS datasets to perform colocalization analysis using the `coloc` and `susieR` packages. This script tests whether the same causal variant influences both protein levels and obesity traits.

3. **TwoStep_MR_analysis.R**  
   Implements a two-step MR mediation analysis to test whether BMI mediates the effect of protein levels on obesity outcomes. It combines precomputed MR results (Protein → BMI and BMI → Outcome) and calculates the mediation effect using Sobel’s method.

4. **Reverse_MR_analysis.R**  
   Conducts reverse MR analysis, treating the outcome trait (BMI) as the exposure and protein levels as the outcome. This script selects genome-wide significant SNPs, formats the data for TwoSampleMR, harmonizes the exposure and outcome datasets, and performs the reverse MR analysis.

## Requirements

- **R (version 3.6 or higher)**
- **Key R packages:**
  - For MR Analysis: `MRInstruments`, `TwoSampleMR`, `MRPRESSO`, `LDlinkR`
  - For Colocalization: `coloc`, `tidyverse`, `GWAS.utils`, `reshape2`, `susieR`
  - For Two-Step MR & Mediation: `dplyr` (via tidyverse)
  - For Reverse MR: `TwoSampleMR`, `tidyverse`

You can install most packages from CRAN. For Bioconductor packages, try:

```r
install.packages(c("tidyverse", "reshape2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MRInstruments", "TwoSampleMR", "MRPRESSO", "LDlinkR", "coloc", "GWAS.utils", "susieR"))
```

## Data Requirements

- **Exposure GWAS Files (Protein Levels):**  
  These files should be tab-delimited with headers. Expected columns include identifiers for the protein, SNP, position, alleles, effect sizes, standard errors, p-values, etc.

- **Outcome GWAS Files (BMI):**  
  Similarly formatted files that include variant IDs, p-values, standard errors, allele frequencies, etc.

- **Additional Files:**  
  Adjust file paths in the scripts to match your local data organization.

## Repository Structure

```
.
├── MR_analysis.R              # Standard Mendelian Randomization analysis
├── Colocalization_analysis.R  # Colocalization analysis script
├── TwoStep_MR_analysis.R      # Two-step MR mediation analysis
├── Reverse_MR_analysis.R      # Reverse MR analysis
└── README.md                  # This README file
```

## How to Use

1. **MR_analysis.R:**  
   - Edit working directory and file paths to point to your exposure and outcome GWAS data.
   - Run the script to perform MR analysis for each protein exposure.
   - Output: A consolidated MR results file (e.g., saved as an RDS file).

2. **Colocalization_analysis.R:**  
   - Update parameters (e.g., SNP of interest, protein name, chromosome) as required.
   - Run the script to merge and process GWAS datasets and perform colocalization analysis.
   - Output: Colocalization summary statistics and, optionally, susie-based coloc results.

3. **TwoStep_MR_analysis.R:**  
   - Ensure MR results from Protein → BMI and BMI → Outcome analyses are available.
   - Run the script to calculate the mediation effect using Sobel’s method.
   - Output: Combined dataset with MR estimates and mediation analysis metrics.

4. **Reverse_MR_analysis.R:**  
   - Adjust file paths and parameters for the exposure (BMI) and protein GWAS data.
   - Run the script to perform reverse MR analysis.
   - Output: A list containing MR results, harmonized datasets, and the target protein identifier.

## Code Availability

The code implemented in these scripts is publicly available under an open-source license. In accordance with the code availability statement in our paper, the complete repository (including this README and all analysis scripts) is hosted on GitHub and archived in Zenodo for long-term preservation. This enables reproducibility and reuse of our methods in future research.

- **GitHub Repository:** [[https://github.com/Raphjacksun7/Proteome_MR_Analysis](https://github.com/Raphjacksun7/Proteome_MR_Analysis)]
- **Zenodo Archive DOI:** [10.5281/zenodo.14969499](https://doi.org/10.5281/zenodo.14969499)
Please refer to the “Code Availability” section in our paper for further details regarding the licensing and citation of the code.

## Citation

If you use this code in your research, please cite our paper:

> Avocegamou R, Jumentier B, Fagbemi K, et al. *Proteome-wide Mendelian randomization identifies circulating proteins causally associated with childhood obesity*. [Journal Name, Year, Volume, Pages]. DOI: [doi link]

## Contact

For any questions or feedback regarding the code or analyses, please contact:

**Despoina Manousaki, MD, PhD**  
Université de Montréal, CHU Sainte Justine  
Email: [despina.manousaki@umontreal.ca](mailto:despina.manousaki@umontreal.ca)
