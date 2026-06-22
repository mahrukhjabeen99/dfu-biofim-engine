# ==============================================================================
# Script Name: Prepare_Outcome_Data.R
# Description: Script to download and format the Diabetic Foot Ulcer (DFU) 
#              GWAS summary statistics from the Pan-UK Biobank AWS repository.
# ==============================================================================

# Install data.table if not already installed for fast reading of large files
if (!require("data.table")) install.packages("data.table")
library(data.table)

# 1. Define the Pan-UKB AWS URL for Diabetic Foot Ulcer (Phecode 705.3)
# Note: Users can use wget, curl, or download.file to fetch this via HTTPS
url <- "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_metrics.txt.bgz"

# Download the file (uncomment the line below to download directly via R)
# download.file(url, destfile = dest_file, method = "curl")

# 2. Read the compressed TSV file
# Using data.table::fread which natively handles .bgz files in modern R environments
dat <- fread(dest_file)

# 3. Format the data to match the TwoSampleMR outcome requirements.
# Pan-UKB provides ancestry-specific columns; we extract the European (EUR) data.
# The alternate allele (alt) is typically used as the effect allele (A1).
gwas_txt <- data.frame(
  SNP = paste0("chr", dat$chr, ":", dat$pos), # Standardizing SNP ID to chr:pos
  A1 = dat$alt,                               # Effect allele
  A2 = dat$ref,                               # Reference allele
  freq = dat$af_EUR,                          # Allele frequency in EUR population
  b = dat$beta_EUR,                           # Effect size (beta) in EUR
  se = dat$se_EUR,                            # Standard error in EUR
  p = dat$pval_EUR,                           # P-value in EUR
  n = 420531                                  # Example Pan-UKB EUR total sample size (adjust if needed)
)

# 4. Clean the dataset (remove variants that were not tested in the EUR population, leading to NAs)
gwas_txt_clean <- na.omit(gwas_txt)

# 5. Save the formatted output as 'gwas.txt' for the downstream MR pipeline
write.table(gwas_txt_clean, file = "gwas.txt", sep = "\t", quote = FALSE, row.names = FALSE)