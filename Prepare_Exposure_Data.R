# ==============================================================================
# Script Name: Prepare_Exposure_Data.R
# Description: Script to parse, filter, and clump the 41 inflammatory cytokine 
#              GWAS files from the Bristol repository (Ahola-Olli et al., 2017)
#              to generate the 'merged_result.csv' and 'exposure.csv' input files.
#              The data could be download from 
#              https://doi.org/10.5523/bris.3g3i5smgghp0s2uvm1doflkx9x
# ==============================================================================

library(TwoSampleMR)

# NOTE: Users must download all 41 .txt files (e.g., meta_filt_B_NGF_noBMI.txt) 
# into the current working directory before running this script.

# Get a list of all 41 cytokine files in the directory
file_list <- list.files(pattern = "meta_filt_.*_noBMI.txt")

all_exposures_clumped <- data.frame()
trait_reference_list <- data.frame()

# Loop through each file to process the data
for (file in file_list) {
  
  # Extract the cytokine name from the file name
  cytokine_name <- gsub("meta_filt_(.*)_noBMI.txt", "\\1", file)
  
  print(paste("Processing:", cytokine_name))
  
  # 1. Read and format the raw data using exact column headers from the Bristol dataset
  dat <- read_exposure_data(
    filename = file,
    sep = "\t",                 # Files are tab-separated
    snp_col = "MarkerName",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    samplesize_col = "Total_N",
    phenotype_col = cytokine_name
  )
  
  # 2. Filter for genome-wide significant SNPs 
  # (Adjust the p-value threshold here based on the study protocol, e.g., 5e-8 or 5e-6)
  dat_sig <- subset(dat, pval.exposure < 5e-8)
  
  # Skip if no significant SNPs are found
  if(nrow(dat_sig) == 0) next
  
  # 3. Perform Linkage Disequilibrium (LD) clumping (r2 = 0.001, kb = 10000)
  dat_clumped <- clump_data(dat_sig)
  
  # 4. Append to the master data frames
  all_exposures_clumped <- rbind(all_exposures_clumped, dat_clumped)
  
  # Save trait metadata for the reference list
  trait_reference_list <- rbind(trait_reference_list, 
                                data.frame(id.exposure = dat_clumped$id.exposure[1], 
                                           exposure = cytokine_name))
}

# 5. Output the final required files for the MR pipeline
write.csv(all_exposures_clumped, "merged_result.csv", row.names = FALSE)
write.csv(unique(trait_reference_list), "exposure.csv", row.names = FALSE)