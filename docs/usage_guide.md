# Detailed Usage Guide

## Table of Contents
1. [Data Preparation](#data-preparation)
2. [Running Analyses](#running-analyses)
3. [Interpreting Results](#interpreting-results)
4. [Expected output](#expected-output)

## Data Preparation

### 1. gRNA database format

Your gRNA database should contain:
- **ID**: Unique gRNA identifier
- **chr**: Chromosome name in UCSC or NCBI format
- **start/end**: Genomic coordinates
- **strand**: + or -
- **protospacer**: 20bp guide sequence
- **pam**: PAM sequence (3bp)
- **cut_site**: Cutting position
- **score_cfd, score_deephf, score_deepspcas9**: Scoring metrics

### 2. VCF files

- Use normalized VCF files (bcftools norm)
- Index with tabix: `tabix -p vcf variants.vcf.gz`
- Supported examples: DeepVariant and pbsv output VCFs

### 3. Homoeologous gene pair table

- **geneA**: First gene in pair
- **geneB**: Second gene in pair
- **perc_identity**: Identitiy pecentage from nucleotide sequence alignment

## Running Analyses

### Batch mode (Command line)

```bash
# Variant impact classification
Rscript scripts/grna_variant_impact_classifier.R \
  --grna data/grna_database.tsv \
  --snpvcf data/variants.vcf.gz \
  --svvcf data/structural_variants.vcf.gz \
  --output results/classification/ \
  --workers 8 \
  --batchsize 10000

# Statistics and vizualization
Rscript scripts/grna_statistics_visualization.R \
  -c results/classification/gRNA_classification_final.csv \
  -g data/grna_database.tsv \
  -p data/homoeolog_pairs_with_identity_percent.tsv \
  -o results/ \
  --cfd_threshold 0.05 \
  --plot_dpi 300
```

### Interactive plotting mode (RStudio)

```r
# Load script
source("scripts/grna_statistics_visualization.R")

# Example 1: Run complete analysis
results <- run_complete_analysis(
   classification_file = "results/classification/gRNA_classification_final.csv",
   grna_db = "data/grna_database.tsv",
   pairs_file = "data/homoeolog_pairs_with_identity_percent.tsv",
   output_dir = "results/",
   cfd_threshold = 0.05,
   save_plots = TRUE,
   save_csv = TRUE
 )

# Example 2: Load data once & generate specific plots
basic_data <- load_basic_data("results/classification/gRNA_classification_final.csv")
full_data <- load_full_data(
  classification_file = "results/classification/gRNA_classification_final.csv",
  grna_db = "data/grna_database.tsv",
  pairs_file = "data/homoeolog_pairs_with_identity_percent.tsv",
  cfd_threshold = 0.05
)

plot_global_category_combo(data = basic_data, save_plot = TRUE)
plot_pair_classification_combo(data = full_data, save_plot = TRUE)
plot_pam_status(data = basic_data, save_plot = TRUE)

threshold_analysis <- prepare_threshold_analysis(
  pairs_dt = full_data$pairs_dt,
  grna_all_dt = full_data$pairs_grna_all,
  grna_filtered_dt = full_data$eligible_grnas,
  thresholds = seq(0.1, 1.0, by = 0.05),  # Custom range
  score_columns = c("score_cfd", "score_deephf")
)
 plot_pair_threshold_sweep(threshold_analysis, save_plot = TRUE)

# etc.
```

## Interpreting Results

### Category definitions

- **Fully_Functional**: No variants affecting gRNA function
- **Reduced_Efficiency**: Variants in distal region or non-canonical PAM
- **Critical_Failure**: Variants in seed region or PAM disruption

### PAM status

- **Canonical**: NGG PAM sequences
- **Non-canonical**: NAG, NGA PAMs (reduced activity)
- **Disrupted**: Non-functional PAM

### Seed status

- **Intact**: No variants in seed region
- **Disrupted**: Critical variants in seed region

### Homoeologous pair classification

- **Both homoeologs**: Both genes have functional specific gRNAs
- **Only one gene in pair**: Only one gene of pair has specific gRNAs
- **No specific gRNA**: Neither gene has specific gRNA

## Expected output

### **Plots** (16 figures total)
- Individual analysis plots (13)
- Combined 2×2 panels (2)
- Threshold sweep analyses (2)

### **Tables** (15+ CSV files)
- Statistical summaries
- Filtered datasets
- Classification results

### **Logs**
- Master analysis log
- Per-process logs
- Error tracking
