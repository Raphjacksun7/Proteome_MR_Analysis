###########################################################################
# Generic Script for Reverse MR Analysis
#
# For proteins that have colocalized with outcomes (BMI) (H4 > 0.8),
# this script performs a reverse MR analysis where the exposure is the 
# outcome trait (BMI) and the outcome is the protein level.
#
# Requirements:
#   - TwoSampleMR and tidyverse packages
#   - Full GWAS data for the exposure trait (BMI)
#   - Full GWAS data for the protein level (as outcome)
###########################################################################

library(TwoSampleMR)
library(tidyverse)

#------------------ Reverse MR Analysis ------------------

# Load full GWAS data for the exposure trait (BMI)
# Replace "path_to_exposure_GWAS.txt" with your actual file path
expo <- read.table("exposure_GWAS", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter for genome-wide significant SNPs (p-value <= 5e-8)
expo <- expo %>% filter(P_value <= 5e-8)

# Format the exposure data for TwoSampleMR
expo <- format_data(expo)

# Clump the exposure dataset to obtain independent SNPs
expo <- clump_data(expo)

# Set the protein name (outcome variable)
protI <- "protein_name"  # Replace with the actual protein name

# Load full GWAS data for the outcome (protein levels)
# Replace "path_to_outcome_GWAS.txt" with your actual file path
outc <- read.table("path_to_outcome_GWAS.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Format the outcome data for TwoSampleMR (type = "outcome")
outc <- format_data(outc, type = "outcome")

# Harmonise the exposure and outcome datasets
dat1 <- harmonise_data(exposure_dat = expo, outcome_dat = outc, action = 1)

# Perform MR analysis on the harmonised data
res1 <- mr(dat1)

# Save the reverse MR results along with the harmonised data and protein identifier
result_list <- list(
  mr        = res1,
  harmonise = dat1,
  prot      = protI
)

# Save the result to an RDS file (adjust the file path as needed)
saveRDS(result_list, "data/reverse_MR_results.rds")
