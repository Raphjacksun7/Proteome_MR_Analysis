###########################################################################
# Generic Script for Mendelian Randomisation (MR) Analysis
#
# This script reads in outcome and exposure data, formats them for the
# TwoSampleMR package, and performs the MR analysis. The consolidated
# results are then saved for further evaluation.
#
# Requirements:
#   - MRInstruments, TwoSampleMR, MRPRESSO, and LDlinkR packages
#   - Outcome file (e.g., BMI GWAS data; see manuscript for download link)
#   - Exposure files (list of files containing protein exposure data with headers)
###########################################################################

# Load required packages
library(MRInstruments)
library(TwoSampleMR)
library(MRPRESSO)
library(LDlinkR)

# Set working directory to where your files are located
setwd("/path")

# Read the outcome file (e.g., BMI GWAS data)
outcome <- read.table("outcome", sep = "\t", header = TRUE)

# Define column names for the results data frame
result_cols <- c("SNP", "Outcome", "Exposure", "Study", "Method", 
                 "Number.SNPs", "Beta", "SE", "Pval")

# Initialize an empty data frame for the consolidated MR results
results <- data.frame(matrix(ncol = length(result_cols), nrow = 0))
colnames(results) <- result_cols

# List of exposure files to process (each file must have a header row)
exposure_files <- c("file1", "file2", "...")

# Loop through each exposure file
for (exposure_file in exposure_files) {
  
  # Read the exposure file as a data frame using its header
  exposure_df <- read.table(exposure_file, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE)
  
  # Optionally, rename columns if they don't match your preferred names.
  # For example, if the columns are not already named appropriately, you might use:
  # colnames(exposure_df) <- c("Protein", "SNP", "col3", "col4", "Pos", 
  #                            "Allele.1", "Allele.2", "EAF", "col9", "Beta", "SE", "Pval")
  #
  # Then, select only the columns needed for the analysis:
  exposure_df <- exposure_df[, c("Protein", "SNP", "Pos", "Allele.1", 
                                 "Allele.2", "EAF", "Beta", "SE", "Pval")]
  
  # Group the exposure data by Protein using the "Protein" column
  exposure_data_list <- split(exposure_df, exposure_df$Protein)
  
  # Process each protein group within the exposure file
  for (protein in names(exposure_data_list)) {
    
    # Each group is already a data frame for that protein
    exposure_data <- exposure_data_list[[protein]]
    
    # Rename columns to match the format required by TwoSampleMR
    exposure_data$beta          <- exposure_data$Beta
    exposure_data$se            <- exposure_data$SE
    exposure_data$effect_allele <- exposure_data$Allele.1
    exposure_data$other_allele  <- exposure_data$Allele.2
    exposure_data$eaf           <- exposure_data$EAF
    exposure_data$id.exposure   <- exposure_data$SNP
    exposure_data$pval          <- exposure_data$Pval
    exposure_data$pval.exposure <- exposure_data$pval
    
    # Clean and convert the 'Pos' column into a numeric 'position' column
    exposure_data$position <- as.numeric(gsub("[\",]", "", exposure_data$Pos))
    
    # Format the exposure data for harmonisation
    exp_data <- format_data(exposure_data, snp_col = "SNP", type = "exposure")
    exp_data$id.exposure <- exposure_data$Protein[1]  # all rows share the same protein
    
    # Save the formatted exposure data to a file
    write.csv(exp_data, paste0(protein, ".tsv"), quote = FALSE, row.names = FALSE)
    
    # Clump the exposure data if there are at least 2 SNPs available
    new_exp_dat <- if (nrow(exp_data) >= 2) clump_data(exp_data) else exp_data
    
    # Prepare the outcome data by renaming columns for TwoSampleMR
    outcome$SNP  <- outcome$variant_id
    outcome$pval <- outcome$p_value
    outcome$se   <- outcome$standard_error
    outcome$eaf  <- outcome$effect_allele_frequency
    
    # Format the outcome data for harmonisation
    out_data <- format_data(outcome, snp_col = "SNP", type = "outcome")
    out_data$id.outcome <- "Trait"
    
    # Harmonise exposure and outcome data
    dat <- harmonise_data(exposure_dat = new_exp_dat, outcome_dat = out_data)
    
    # Perform MR analysis
    res <- mr(dat)
    
    # Optional: Uncomment and adjust the following lines to perform the Steiger Test
    # r_exp <- get_r_from_bsen(b = beta_expo, se = se_expo, n = N_expo)
    # r_out <- get_r_from_bsen(b = beta_out, se = se_out, n = N_out)
    # res <- mr_steiger(p_exp = p_val_expo, p_out = p_val_out, n_exp = N_expo,
    #                   n_out = N_out, r_exp = r_exp, r_out = r_out)
    
    # Skip to the next protein group if no results were produced
    if (length(res) == 0) next
    
    # Create a data frame for the MR results for the current protein
    new_res <- data.frame(
      SNP         = dat$SNP, 
      Outcome     = res$id.outcome, 
      Exposure    = res$id.exposure, 
      Study       = exposure_file, 
      Method      = res$method, 
      Number.SNPs = res$nsnp, 
      Beta        = res$b, 
      SE          = res$se, 
      Pval        = res$pval, 
      stringsAsFactors = FALSE
    )
    
    # Append the new results to the consolidated results data frame
    results <- rbind(results, new_res)
  }
}

# Save the consolidated MR results to file (adjust the path and filename as needed)
# saveRDS(results, "data/MR_results.rds")
