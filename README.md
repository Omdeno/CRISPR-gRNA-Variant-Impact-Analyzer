# CRISPR gRNA Variant Impact Analyzer

## Overview

This repository contains R scripts for comprehensive analysis and visualization of CRISPR gRNA functionality with focus on variant impact assessment and homoeologous gene targeting.

Created for article "Optimization of functional genetics tools for a model tetraploid Capsella bursa-pastoris, with focus on homoeolog-aware gene editing"

## Features

### 1. Variant Impact Classification (`grna_variant_impact_classifier.R`)

- **Multi-Variant Support**: SNPs, MNPs, indels, and structural variants
- **Parallel Processing**: Efficient batch processing for large datasets
- **Comprehensive Logging**: Detailed process tracking and error reporting
- **Flexible Input**: VCF format support (DeepVariant, pbsv, etc.), gRNA library as table in TSV format (based on flattened output of crisprDesign (https://github.com/crisprVerse/crisprDesign))
- **Detailed Classification**: Functional, reduced efficiency, or critical failure

### 2. Statistical Analysis & Visualization (`grna_statistics_visualization.R`)

- **Global gRNA Classification**: Distribution analysis across functional categories based on Variant Impact Classification output
- **PAM Status Analysis**: Canonical, non-canonical, and disrupted SpCas9 PAM sequences
- **Seed Region Integrity**: Assessment of critical targeting regions
- **Homoeolog Analysis**: Gene pair targeting efficiency evaluation
- **Threshold Optimization**: Score-based filtering impact on gRNA retention (aggregated CFD, DeepHF and DeepSpCas9)
- **Publication-Ready Figures**: High-resolution plots using ggpubr

## Installation

### Prerequisites
R 4.4.2

Bioconductor 3.20

```r
# Install required CRAN packages
install.packages(c(
  "optparse", "data.table", "ggplot2", "scales",
  "forcats", "splitstackshape", "patchwork", "dplyr", 
  "tidyr", "purrr", "stringr", "ggpubr", "logger"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "VariantAnnotation", "StructuralVariantAnnotation", 
  "Biostrings", "GenomicRanges", "IRanges", 
  "BiocParallel", "GenomeInfoDb"
))
