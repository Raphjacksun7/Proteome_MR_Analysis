###########################################################################
# Generic Script for Colocalization Analysis
#
# This script performs colocalization analysis for a significant protein
# from our MR study by merging exposure and outcome GWAS datasets.
# It uses the coloc and susieR packages to estimate the probability of 
# shared causal variants between protein levels and the outcome trait.
#
# Requirements:
#   - coloc, tidyverse, GWAS.utils, reshape2, susieR packages
#   - Exposure GWAS file for protein levels (with headers)
#   - Outcome GWAS file (pre-cut by chromosome)
###########################################################################

# Load required libraries
library(coloc)
library(tidyverse)
library(GWAS.utils)
library(reshape2)
library(susieR)

# Initialize an empty data frame to store colocalization results
tab_res <- data.frame()

# ====================== TO MODIFY ==========================
# Set file paths and parameters for the analysis

# Exposure GWAS file (protein levels)
exposure_file <- "PROTEIN_GWAS"  # Replace with your actual file path

# Define SNP (rsid) of interest, protein name, chromosome, and study name
rsidI  <- "NAME_OF_RSID_OF_INTEREST"   # e.g., "rs123456"
protI  <- "NAME_OF_PROTEIN"             # e.g., "ProteinX"
chrI   <- "XXX"                        # e.g., "1"
studyI <- "Emilsson"                   # Optional: Name of GWAS study
# ============================================================

# ---------------------- Automatic Steps ---------------------

# Load the exposure GWAS data (assumed to have header with columns such as "rsids" and "Pos")
expo <- read.table(exposure_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Identify the position of the SNP of interest from the exposure GWAS
posI <- expo$Pos[expo$rsids == rsidI]

# Define a 1 Mb window around the SNP of interest (±500 kb)
pos_min <- posI - 500000
pos_max <- posI + 500000

# Subset the exposure data within the specified window
expo_cut <- expo %>% filter(Pos >= pos_min, Pos <= pos_max)

# Load the outcome GWAS data for the corresponding chromosome
# Assumes outcome file is pre-cut by chromosome and named accordingly.
outcome_file <- paste0("GWAS_OUTCOME_chr", chrI, ".rds")  # Adjust as needed
outc <- readRDS(outcome_file)

# Merge the two GWAS datasets by SNP ID:
# Assumption: expo_cut uses column "rsids" for SNP IDs and the outcome data uses "rsid".
tmp <- merge(expo_cut, outc, by.x = "rsids", by.y = "rsid")

# Rearrange and rename columns to standard names.
# Assumed columns in the merged data:
#   - "rsids": SNP identifier (to be renamed to "rsid")
#   - "chr": Chromosome
#   - "Pos": Position
#   - "eff.x", "se.x", "pv.x": Exposure effect, standard error, and p-value
#   - "N": Exposure sample size
#   - "maf": Exposure minor allele frequency
#   - "eff.y", "se.y", "pv.y": Outcome effect, standard error, and p-value
tmp <- tmp %>%
  rename(
    rsid       = rsids,
    position   = Pos,
    expo_eff   = eff.x,
    expo_se    = se.x,
    expo_pv    = pv.x,
    expo_N     = N,
    expo_maf   = maf,
    outc_eff   = eff.y,
    outc_se    = se.y,
    outc_pv    = pv.y
  ) %>%
  select(rsid, chr, position, expo_eff, expo_se, expo_pv, expo_N, expo_maf,
         outc_eff, outc_se, outc_pv)

# Calculate the variance of beta estimates (square of the standard errors)
tmp <- tmp %>%
  mutate(
    expo_varbeta = expo_se^2,
    outc_varbeta = outc_se^2
  )

# Print the row corresponding to the SNP of interest and the row with the lowest outcome p-value
print(tmp %>% filter(rsid == rsidI))
print(tmp %>% filter(outc_pv == min(outc_pv, na.rm = TRUE)))

# ------------------ (Optional) GWAS Plotting ------------------
# Uncomment and adjust the following code to plot GWAS data if needed.
#
# dat <- tmp %>% select(rsid, position, expo_pv, outc_pv) %>% 
#   melt(id.vars = c("rsid", "position"))
#
# dat <- dat %>%
#   mutate(
#     lab = ifelse(rsid == rsidI, rsid, NA),
#     size = ifelse(rsid == rsidI, 1, 0.5),
#     variable = recode(variable, expo_pv = protI, outc_pv = "Outcome")
#   )
#
# ggplot(dat, aes(x = position, y = -log10(value), color = lab, size = size)) +
#   geom_vline(xintercept = posI) +
#   geom_point() +
#   facet_wrap(~ variable, nrow = 2) +
#   ggtitle(paste0("RSID of Interest = ", rsidI)) +
#   xlab("Chromosomal Position") +
#   theme_bw() +
#   theme(legend.position = "none") +
#   scale_size_manual(values = c(0.5, 2))
#
# --------------------------------------------------------------

# ---------------------- COLOCALIZATION ANALYSIS ---------------------

# Prepare dataset for the exposure (protein levels)
D1_expo <- list(
  type     = "quant",
  snp      = tmp$rsid,
  position = tmp$position,
  beta     = tmp$expo_eff,
  varbeta  = tmp$expo_varbeta,
  sdY      = 1  # Protein levels are scaled; adjust if needed
)
check_dataset(D1_expo)

# Prepare dataset for the outcome
D2_outc <- list(
  type     = "quant",
  snp      = tmp$rsid,
  position = tmp$position,
  beta     = tmp$outc_eff,
  varbeta  = tmp$outc_varbeta,
  sdY      = "XXX"  # Replace "XXX" with the correct standard deviation (e.g., 1.3 for AAM, 4 for ANM)
)
check_dataset(D2_outc)
# Optionally, if any varbeta equals 0, replace them with the mean:
# D2_outc$varbeta[D2_outc$varbeta == 0] <- mean(D2_outc$varbeta)

# Run colocalization analysis using the coloc.abf method
res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)

# Append the colocalization summary results to tab_res
tab_res <- rbind(tab_res, c(studyI, protI, rsidI, as.vector(res$summary)))
colnames(tab_res) <- c("Study", "Gene", "Rsid", "nSNP", "H0", "H1", "H2", "H3", "H4")

# Convert hypothesis columns to numeric
tab_res <- tab_res %>% mutate_at(vars(H0:H4), as.numeric)

# Save the colocalization results (adjust the file path as needed)
# saveRDS(tab_res, "data/coloc_results.rds")

# ---------------------- SUSIE PLUG-IN ---------------------
# Read or generate the LD matrix for the tested SNPs.
ld <- "read_ld_files"  # Replace with actual code to load your LD matrix

# Calculate Z-scores for exposure and outcome data
z_expo <- tmp$expo_eff / tmp$expo_se
z_outc <- tmp$outc_eff / tmp$outc_se

# Estimate lambda (concordance between the LD matrix and Z-scores)
lambda_expo <- estimate_s_rss(z_expo, ld, n = tmp$expo_N[1])  # Assumes constant N across SNPs; adjust if needed
lambda_outc <- estimate_s_rss(z_outc, ld, n = tmp$expo_N[1])  # Replace with outcome N if available

# Run susie_rss for exposure and outcome data
fit_expo <- susie_rss(z_expo, ld, n = tmp$expo_N[1])
fit_outc <- susie_rss(z_outc, ld, n = tmp$expo_N[1])  # Replace with appropriate outcome sample size

# Run coloc.susie analysis
res_susie <- coloc.susie(fit_expo, fit_outc)

# Save susie coloc results (adjust the file path as needed)
# saveRDS(res_susie, "data/coloc_susie_results.rds")
