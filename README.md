# Data and code from: Exploring the causal relationship between inflammatory cytokines and diabetic foot ulcer: A bidirectional Mendelian randomization and clinical validation study

Dataset DOI: [10.5061/dryad.k98sf7mp4](https://doi.org/10.5061/dryad.k98sf7mp4)

## Description of the data and file structure

This dataset contains the primary data and code supporting the study "Exploring the causal relationship between inflammatory cytokines and diabetic foot ulcer: A bidirectional Mendelian randomization and clinical validation study." The study investigates the potential causal pathways between systemic inflammatory markers and the development of diabetic foot ulcers (DFU).

The experimental efforts and data collection encompass two main phases:

1. **Bidirectional Mendelian Randomization (MR):** We utilized publicly available Genome-Wide Association Study (GWAS) summary statistics to perform a two-sample bidirectional MR analysis. This phase aimed to infer the causal effects of various inflammatory cytokines (as exposures) on the risk of DFU (as the outcome), and vice versa, while rigorously screening for pleiotropy and heterogeneity.
2. **Clinical Validation:** To corroborate the genetic findings, we collected and analyzed clinical data from a cohort of patients. This phase involved measuring the serum levels of the causally implicated cytokines and assessing their correlation with clinical characteristics and wound healing outcomes in DFU patients.

The provided data files and R scripts allow for the full reproduction of the causal inference models and the clinical association analyses reported in the manuscript.

### Files and variables

#### File: Raw_data_generated_during_the_study.xlsx

**Description:** This Excel workbook contains four spreadsheets summarizing the final outputs of the Mendelian Randomization analyses (both forward and reverse), including effect estimates, heterogeneity, and horizontal pleiotropy test results.

##### Variables:

**Sheet 1:** **`table.heterogeneity`** Contains results from Cochran's Q test to assess heterogeneity among the instrumental variables (SNPs).

* **`id.exposure`** / **`id.outcome`**: Unique identifier for the exposure / outcome dataset.
* **`exposure`** / **`outcome`**: Name of the exposure / outcome trait.
* **`method`**: The MR method used for the heterogeneity test (e.g., MR Egger, Inverse variance weighted).
* **`Q`**: Cochran's Q statistic indicating the degree of heterogeneity.
* **`Q_df`**: Degrees of freedom associated with the Q statistic.
* **`Q_pval`**: P-value for Cochran's Q test (P < 0.05 indicates significant heterogeneity).
* **Sheet 2:** **`table.MRresult`** Contains the primary causal effect estimates of the exposure (inflammatory cytokines) on the outcome (diabetic foot ulcer) using various MR methods.
  * **`id.exposure`**, **`id.outcome`**, **`exposure`**, **`outcome`**: Identifiers and names for the analyzed traits.
  * **`method`**: The MR analytical method applied (e.g., Inverse variance weighted, Weighted median).
  * **`nsnp`**: Number of Single Nucleotide Polymorphisms (SNPs) used as instrumental variables.
  * **`b`**: Estimated causal effect size (beta coefficient).
  * **`se`**: Standard error of the causal effect estimate.
  * **`pval`**: P-value for the causal effect estimate.
  * **`lo_ci`** / **`up_ci`**: Lower and upper bounds of the 95% confidence interval for the beta coefficient.
  * **`or`**: Odds Ratio (OR) representing the causal effect.
  * **`or_lci95`** / **`or_uci95`**: Lower and upper bounds of the 95% confidence interval for the Odds Ratio.
  * **Sheet 3:** **`table.pleiotropy`** Contains results from the MR-Egger regression intercept test to assess directional horizontal pleiotropy.
    * **`id.exposure`**, **`id.outcome`**, **`exposure`**, **`outcome`**: Identifiers and names for the analyzed traits.
    * **`egger_intercept`**: The estimated intercept from the MR-Egger regression. A value significantly deviating from zero indicates the presence of directional pleiotropy.
    * **`se`**: Standard error of the MR-Egger intercept.
    * **`pval`**: P-value for the MR-Egger intercept test (P < 0.05 indicates significant directional pleiotropy).
    * **Sheet 4:** **`Reverse MR analysis`** Contains the causal effect estimates for the reverse MR analysis (where diabetic foot ulcer is the exposure and inflammatory cytokines are the outcomes).
      * **Variables** **`id.exposure`** **through** **`or_uci95`**: Identical definitions as listed in `table.MRresult` above, but reflecting the swapped exposure and outcome roles.
      * **`file_source`**: The origin or source identifier of the data file used for this specific reverse MR iteration.

#### File: MR_code.R

**Description:** R script containing the complete analytical pipeline for the forward Mendelian Randomization analysis. It performs data harmonization, calculates causal estimates (Inverse Variance Weighted, etc.), and generates sensitivity analyses (pleiotropy and heterogeneity tests) alongside visualization plots (scatter, forest, funnel, and leave-one-out plots).

##### Exposure Data Variables:

* **`SNP`**: Single Nucleotide Polymorphism identifier, typically an rsID (e.g., rs123456).
* **`beta.exposure`**: The effect size (beta coefficient) representing the direction and magnitude of the SNP's association with the exposure trait.
* **`se.exposure`**: Standard error of the beta coefficient.
* **`pval.exposure`**: The p-value for the genetic association test between the SNP and the exposure.
* **`effect_allele.exposure`**: The allele being tested or modeled for its effect (also known as the risk allele or coded allele).
* **`other_allele.exposure`**: The non-effect, alternative, or reference allele.
* **`eaf.exposure`**: Effect Allele Frequency; the frequency of the effect allele in the study population.
* **`exposure`**: The name of the exposure trait (e.g., specific inflammatory cytokines).
* **`id.exposure`**: A unique identifier for the exposure dataset.
* **`samplesize.exposure`**: The total number of participants (sample size) included in the exposure GWAS.

##### Outcome Data Variables

* **`SNP`**: Single Nucleotide Polymorphism identifier (rsID).
* **`beta.outcome`** **(`b`)**: The effect size (beta coefficient) of the SNP on the outcome trait (Diabetic foot ulcer).
* **`se.outcome`** **(`se`)**: Standard error of the outcome beta coefficient.
* **`effect_allele.outcome`** **(`A1`)**: The effect allele for the outcome. *Note: Must be carefully harmonized with the exposure effect allele.*
* **`other_allele.outcome`** **(`A2`)**: The reference allele for the outcome.
* **`pval.outcome`** **(`p`)**: The p-value for the association between the SNP and the outcome.
* **`eaf.outcome`** **(`freq`)**: Frequency of the effect allele in the outcome GWAS population.
* **`samplesize.outcome`** **(`n`)**: The total sample size of the outcome GWAS.
* **`outcome`**: The name of the outcome trait ("Diabetic foot ulcer").
* **`id.outcome`**: A unique identifier for the outcome dataset.
* **File Generation and Script Dependency:** Crucially, executing this script automatically generates an intermediate file named `outcome.csv` in the working directory. This file contains the subset of extracted and harmonized summary statistics from the outcome dataset. This dynamically generated `outcome.csv` serves as a mandatory input file for the subsequent reverse Mendelian Randomization script (`Reverse_MR_code.R`), establishing a direct sequential workflow dependency between the two scripts.

#### File: reverse_MR_code.R

**Description:** R script used for the reverse Mendelian Randomization analysis. It dynamically sets the diabetic foot ulcer data as the exposure and the inflammatory cytokines as the outcomes to test for reverse causality, filtering out traits with a significant reverse causal effect (P < 0.05).

##### Variables:

The core analytical variables and data structures generated in this script are identical to those described in the `MR code.R` section above. The only difference is that the roles of the traits are reversed: Diabetic foot ulcer serves as the exposure, and the inflammatory cytokines serve as the outcomes.

#### File: Prepare_Outcome_Data.R

**Description:** This script automates the formatting and data cleaning process for the outcome variable, Diabetic Foot Ulcer (DFU). It processes the raw, large-scale Genome-Wide Association Study (GWAS) summary statistics flat file from the Pan-UK Biobank, extracting and aligning the specific columns required for downstream Mendelian Randomization (MR) analysis.

**Generated Outputs:** Running this script will automatically generate the following file in the working directory:

* `gwas.txt`: The fully formatted and cleaned dataset containing the DFU summary statistics. This file serves as the mandatory outcome data input for the primary forward MR analysis script (`MR code.R`).

#### File: Prepare_Exposure_Data.R

* **Description:** This script automates the data extraction, formatting, and filtering pipeline for the exposure variables (inflammatory cytokines). It processes 41 raw Genome-Wide Association Study (GWAS) summary statistic files, applies genome-wide significance thresholds, and performs Linkage Disequilibrium (LD) clumping to synthesize the instrumental variables required for the Mendelian Randomization (MR) analysis.

  **Data Source & Prerequisites:** Due to file size constraints, the raw summary statistics are not hosted in this repository. Before executing this script, users must:
  * Download the 41 meta-analysis summary statistic files.
  * Place all downloaded `.txt` files directly into the same local working directory as this R script.
  * Ensure the `TwoSampleMR` R package is installed and loaded.
  * **Generated Outputs:** Running this script will automatically generate two files in the working directory:
    * `merged_result.csv`: The complete, clumped dataset containing the instrumental variables for all valid cytokines. This is the primary input file required for the forward MR analysis (`MR_code.R`).
    * `exposure.csv`: A reference list mapping the exposure IDs to their respective cytokine names. This file is required to iterate through traits during the reverse MR analysis (`reverse_MR_code.R`).

## Code/software

**1. Software to View Data:**

* **Tabular Data (`.csv`,** **`.txt`):** All tabular datasets and results can be viewed using any standard spreadsheet software, including Microsoft Excel, Apple Numbers, or open-source alternatives like LibreOffice Calc and Apache OpenOffice. `.txt` and `.csv` files can also be viewed in any basic text editor (e.g., Notepad, TextEdit).
* **Excel Data (`.xlsx`):** The summary result file (`Raw_data_generated_during_the_study.xlsx`) requires Microsoft Excel or an open-source equivalent (e.g., LibreOffice Calc) to view all sheets properly.
* **Code Files (`.R`):** The R scripts can be opened and viewed in any text editor, though an Integrated Development Environment (IDE) like **RStudio** is highly recommended for reading and execution.

**2. Software to Run Code:** To reproduce the Mendelian Randomization analyses, the following software and environment are required:

* **R Environment:** R version 4.3.1 (or higher).
* **IDE:** RStudio (recommended for workflow management).
* **Loaded R Packages:**
  * `TwoSampleMR`  – Core package for Mendelian Randomization analysis.
  * `VariantAnnotation` – Required for parsing and formatting genetic variants.
  * `gwasglue` – Used for connecting GWAS data with the TwoSampleMR pipeline.

**3. Included Scripts and Workflow:** These R scripts are included in this submission to fully reproduce the analytical pipeline:

* **Data Preparation Scripts:**
  * **`Prepare_Exposure_Data.R`**: Automates the extraction, formatting, and LD-clumping of the inflammatory cytokine GWAS data to generate the required exposure input files (`merged_result.csv` and `exposure.csv`).
  * **`Prepare_Outcome_Data.R`**: Formats the downloaded Pan-UK Biobank summary statistics into the standardized layout required for the outcome dataset (`gwas.txt`).
  * **Primary Analysis Scripts:**
  * **`MR code.R`**: This is the primary script. It loads the formatted exposure and outcome datasets, performs data harmonization, executes the forward MR analysis (estimating the causal effect of cytokines on diabetic foot ulcer), conducts sensitivity tests (heterogeneity and pleiotropy), and generates all visualization plots (scatter, forest, funnel, leave-one-out).
* **`Reverse_MR_code.R`**: This script performs the reverse causality analysis. It dynamically sets diabetic foot ulcer as the exposure and iteratively tests its effect against the inflammatory cytokines (as outcomes). It outputs the results and filters out traits that show significant reverse causality.

## Access information

Data was derived from the following sources: Genetic data on inflammatory cytokines can be accessed freely at [https://doi.org/10.5523/bris.3g3i5smgghp0s2uvm1doflkx9x](https://doi.org/10.5523/bris.3g3i5smgghp0s2uvm1doflkx9x). Genetic data on diabetic foot ulcer can be accessed freely at Pan-UK Biobank ([https://pan-ukb.org](https://www.google.com/search?q=https://pan-ukb.org)).

**License Information**

In accordance with Dryad's policies, all original clinical data and code provided in this repository are dedicated to the public domain under a **CC0 1.0 Universal (CC0 1.0) Public Domain Dedication** waiver. Users are free to copy, modify, distribute, and perform the work, even for commercial purposes, all without asking permission. The derived GWAS summary statistics were originally published under open-access terms compatible with the CC0 waiver.

## Human subjects data

All participants involved in the clinical validation phase of this study provided written informed consent, which included explicit consent for their de-identified clinical data to be published in the public domain for research purposes.

To protect participant privacy and comply with ethical guidelines, the dataset was fully anonymized prior to submission. All personally identifiable information (PII) and protected health information (PHI)—including names, birth dates, addresses, contact information, and specific hospital medical record numbers—were completely removed. Each participant is identified only by a randomly assigned, sequential study ID (e.g., Patient_1, Patient_2). 