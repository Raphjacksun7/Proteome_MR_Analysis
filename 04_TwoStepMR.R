###########################################################################
# Generic Script for Two Step MR Analysis with Mediation Testing
#
# For proteins that have colocalized with outcomes (%BF) (H4 > 0.8),
# we test whether BMI mediates the relationship between protein levels 
# and the outcome.
#
# For BMI, adult %BF GWAS is used.
#
# Requirements:
#   - Previously computed MR results for:
#       (a) Protein --> BMI (mr_a)
#       (b) BMI --> Outcome: %BF (mr_b)
#
# The mediation function is adapted from the bda package.
# See our Manuscript for GWAS download links.
###########################################################################

library(dplyr)

#------------------ Mediation Function ------------------
# Calculates the indirect effect (a * b) along with its standard error,
# 95% confidence intervals, z-score, and p-value using Sobel's method.
mediation <- function(a, sa, b, sb) {
  # Estimate indirect effect
  ab <- a * b
  
  # Calculate variance of the indirect effect using Sobel's formula
  tmp1 <- b^2 * sa^2 + a^2 * sb^2
  
  # Calculate z-score and two-sided p-value
  zsob <- ab / sqrt(tmp1)
  psob <- 2 * pnorm(-abs(zsob))
  
  # Estimate standard error from beta and z-score
  se <- ab / zsob
  
  # Calculate 95% confidence intervals
  CI_sup <- ab + se * qnorm(0.975)
  CI_inf <- ab - se * qnorm(0.975)
  
  # Return results as a named vector
  res <- c(ab = ab, se = se, CI_inf = CI_inf, CI_sup = CI_sup, z = zsob, pv = psob)
  return(res)
}

#------------------ MR Results Preparation ------------------
# NOTE: 'mr_a' and 'mr_b' should be data frames containing your MR results.
# mr_a: Results from the MR analysis (Protein --> BMI)
# mr_b: Results from the MR analysis (BMI --> Outcome: %BF)

# Filter to include only the Inverse Variance Weighted (IVW) estimates
mr_a_ivw <- mr_a %>% filter(method == "Inverse variance weighted")
mr_b_ivw <- mr_b %>% filter(method == "Inverse variance weighted")

# Subset and rename columns for clarity.
# For mr_a, select columns 1:3 and then columns 8:12.
# Adjust these indices if your data frame structure differs.
mr_a_ivw <- mr_a_ivw[, c(1:3, 8:12)]
colnames(mr_a_ivw)[4:8] <- c("A_method", "A_nsnp", "A_beta", "A_se", "A_pval")

# For mr_b, select columns 5:9 and rename them.
mr_b_ivw <- mr_b_ivw[, 5:9]
colnames(mr_b_ivw) <- c("B_method", "B_nsnp", "B_beta", "B_se", "B_pval")

# Combine the MR results side-by-side
mr_combined <- cbind(mr_a_ivw, mr_b_ivw)

#------------------ Mediation Analysis ------------------
# Loop over each SNP/instrument in the combined MR results and 
# calculate the mediation effect.
mediation_results <- apply(mr_combined, 1, function(row) {
  mediation(a = as.numeric(row["A_beta"]), 
            sa = as.numeric(row["A_se"]), 
            b = as.numeric(row["B_beta"]), 
            sb = as.numeric(row["B_se"]))
})

# Convert mediation results to a data frame (transposing rows and columns)
mediation_results <- as.data.frame(t(mediation_results))

# Combine the MR estimates with the mediation results
RES <- cbind(mr_combined, mediation_results)

# Save the combined results to an RDS file
# Uncomment and adjust the file path as needed.
# saveRDS(RES, "data/mediation_results.rds")
