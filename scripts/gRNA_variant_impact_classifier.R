#!/usr/bin/env Rscript
# =============================================================================
# gRNA Variant Impact Classifier
# =============================================================================
# Description: Comprehensive classification of CRISPR gRNA functionality
#              based on genomic variants (SNPs, indels, structural variants)
# Author: Denis Omelchenko
# Date: 2025
# Version: 1.0.0
# License: MIT
# =============================================================================

# ---- ENVIRONMENT SETUP ------------------------------------------------------

options(
  error = function() {
    traceback(3)
    quit(save = "no", status = 1)
  },
  scipen = 999,
  stringsAsFactors = FALSE
)

# Set locale for consistent behavior
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# ---- PACKAGE MANAGEMENT -----------------------------------------------------

#' Load Required Packages
#'
#' @description Installs and loads all required packages
#' @param packages Character vector of package names
#' @return NULL (packages loaded as side effect)
load_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
      # Bioconductor packages
      if (pkg %in% c(
        "VariantAnnotation", "StructuralVariantAnnotation", "Biostrings",
        "GenomicRanges", "IRanges", "BiocParallel", "GenomeInfoDb"
      )) {
        if (!suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        # CRAN packages
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

REQUIRED_PACKAGES <- c(
  "VariantAnnotation", "StructuralVariantAnnotation", "Biostrings",
  "GenomicRanges", "IRanges", "BiocParallel", "data.table",
  "logger", "optparse", "GenomeInfoDb"
)

load_required_packages(REQUIRED_PACKAGES)

# ---- COMMAND LINE ARGUMENTS -------------------------------------------------

option_list <- list(
  make_option(
    c("-w", "--workers"), 
    type = "integer",
    default = max(1L, parallel::detectCores(logical = TRUE) - 2L),
    help = "Number of worker processes [default: %default]"
  ),
  make_option(
    c("-g", "--grna"), 
    type = "character",
    default = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    help = "Path to gRNA database file [default: %default]"
  ),
  make_option(
    c("-s", "--snpvcf"), 
    type = "character",
    default = "variants.vcf.gz",
    help = "Path to SNP/indel VCF file [default: %default]"
  ),
  make_option(
    c("-v", "--svvcf"), 
    type = "character",
    default = "structural_variants.vcf.gz",
    help = "Path to structural variant VCF file [default: %default]"
  ),
  make_option(
    c("-o", "--output"), 
    type = "character",
    default = "gRNA_classification_results",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("-b", "--batchsize"), 
    type = "integer",
    default = 10000,
    help = "Batch size for progress logging [default: %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate worker count
total_cores <- parallel::detectCores(logical = TRUE)
opt$workers <- max(1L, min(opt$workers, total_cores))
progress_step <- max(1L, opt$batchsize)

# ---- LOGGING SETUP ----------------------------------------------------------

#' Prepare Logging Infrastructure
#'
#' @param output_dir Output directory for logs
#' @return NULL (logging configured as side effect)
prepare_logging <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "process_logs"), showWarnings = FALSE)
  dir.create(file.path(output_dir, "error_logs"), showWarnings = FALSE)
  
  log_file <- file.path(output_dir, "master_analysis.log")
  log_appender(appender_file(log_file))
  log_threshold(INFO)
}

prepare_logging(opt$output)

#' Safe Logging Function
#'
#' @param level Log level (DEBUG, INFO, WARN, ERROR)
#' @param message Log message
#' @param process_id Optional process ID for worker logs
#' @return NULL (logs written as side effect)
safe_log <- function(level, message, process_id = NULL) {
  tryCatch({
    if (!is.null(process_id)) {
      proc_log <- file.path(
        opt$output, "process_logs", sprintf("process_%s.log", process_id)
      )
      cat(sprintf("[%s] %s: %s\n", Sys.time(), level, message),
          file = proc_log, append = TRUE)
    } else {
      switch(
        level,
        "DEBUG" = log_debug(message),
        "INFO"  = log_info(message),
        "WARN"  = log_warn(message),
        "ERROR" = log_error(message),
        log_info(message)
      )
    }
  }, error = function(e) {
    cat(sprintf("[%s] %s: %s\n", Sys.time(), level, message), "\n")
  })
}

# ---- SEQUENCE MANIPULATION FUNCTIONS ----------------------------------------

# Cache for reverse complement sequences
REVERSE_COMPLEMENT_CACHE <- new.env(hash = TRUE, size = 10000L)

#' Get Reverse Complement of DNA Sequence
#'
#' @param sequence DNA sequence string
#' @return Reverse complement sequence
get_reverse_complement <- function(sequence) {
  if (is.null(sequence) || is.na(sequence) || sequence == "") {
    return("")
  }
  
  sequence <- toupper(sequence)
  
  # Check cache
  if (exists(sequence, envir = REVERSE_COMPLEMENT_CACHE, inherits = FALSE)) {
    return(get(sequence, envir = REVERSE_COMPLEMENT_CACHE, inherits = FALSE))
  }
  
  # Calculate reverse complement
  rc <- tryCatch(
    {
      if (grepl("[^ACGTN]", sequence)) {
        stop("Non-DNA character in sequence")
      }
      as.character(reverseComplement(DNAString(sequence)))
    },
    error = function(e) sequence
  )
  
  # Cache result
  assign(sequence, rc, envir = REVERSE_COMPLEMENT_CACHE)
  return(rc)
}

#' Convert Sequence to gRNA Orientation
#'
#' @param seq DNA sequence
#' @param strand Strand (+ or -)
#' @return Sequence in gRNA orientation
to_grna_orientation <- function(seq, strand) {
  if (is.null(seq) || is.na(seq) || seq == "") {
    return(seq)
  }
  
  seq <- toupper(seq)
  
  if (strand == "-") {
    return(get_reverse_complement(seq))
  }
  
  return(seq)
}

#' Apply Variant to Sequence
#'
#' @description Applies SNP/MNP variant to a sequence
#' @param sequence Original sequence
#' @param variant_pos Variant position (genomic)
#' @param variant_ref Reference allele
#' @param variant_alt Alternative allele
#' @param sequence_start Start position of sequence (genomic)
#' @return Modified sequence or original if variant cannot be applied
apply_variant_to_sequence <- function(
    sequence, 
    variant_pos, 
    variant_ref, 
    variant_alt, 
    sequence_start
) {
  
  if (is.null(sequence) || is.na(sequence) || sequence == "") {
    return(sequence)
  }
  
  if (any(is.na(c(variant_pos, variant_ref, variant_alt)))) {
    return(sequence)
  }
  
  # Normalize sequences
  sequence <- toupper(sequence)
  variant_ref <- toupper(variant_ref)
  variant_alt <- toupper(variant_alt)
  
  # Calculate relative position
  rel_pos <- variant_pos - sequence_start + 1L
  
  if (rel_pos < 1L) {
    return(sequence)
  }
  
  # Only handle SNPs and MNPs (same length)
  ref_len <- nchar(variant_ref)
  alt_len <- nchar(variant_alt)
  
  if (ref_len != alt_len) {
    return(sequence)
  }
  
  # Split sequence into characters
  seq_chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  end_pos <- rel_pos + ref_len - 1L
  
  if (end_pos > length(seq_chars)) {
    return(sequence)
  }
  
  # Verify reference allele matches
  target_segment <- paste(seq_chars[rel_pos:end_pos], collapse = "")
  
  if (target_segment != variant_ref) {
    return(sequence)
  }
  
  # Apply variant
  alt_chars <- strsplit(variant_alt, "", fixed = TRUE)[[1]]
  seq_chars[rel_pos:end_pos] <- alt_chars
  
  return(paste(seq_chars, collapse = ""))
}

# ---- CHROMOSOME STYLE CONVERSION --------------------------------------------

#' Convert Chromosome Style
#'
#' @param chr_vec Vector of chromosome names
#' @param target_style Target style ("UCSC" or "NCBI")
#' @return Converted chromosome names
convert_chr_style <- function(chr_vec, target_style) {
  if (is.null(target_style) || length(chr_vec) == 0L) {
    return(chr_vec)
  }
  
  chr_vec <- as.character(chr_vec)
  na_idx <- is.na(chr_vec) | chr_vec == ""
  result <- chr_vec
  
  if (target_style == "UCSC") {
    # Add "chr" prefix
    tmp <- sub("^chr", "", chr_vec, perl = TRUE)
    tmp[tmp %in% c("M", "MT", "chrM", "chrMT", "Mito")] <- "M"
    add_chr <- !na_idx & !grepl("^chr", chr_vec)
    tmp[na_idx] <- chr_vec[na_idx]
    result[add_chr] <- paste0("chr", tmp[add_chr])
    result[!add_chr & !na_idx] <- chr_vec[!add_chr & !na_idx]
    result[result == "chrMT"] <- "chrM"
    result[!grepl("^chr", result) & !na_idx] <- paste0("chr", result[!grepl("^chr", result) & !na_idx])
    
  } else if (target_style == "NCBI") {
    # Remove "chr" prefix
    result <- sub("^chr", "", chr_vec, perl = TRUE)
    result[result %in% c("M", "MT", "chrM", "chrMT", "Mito")] <- "MT"
  }
  
  result[na_idx] <- chr_vec[na_idx]
  return(result)
}

#' Get Sequence Style from Chromosome Names
#'
#' @param chromosomes Vector of chromosome names
#' @return "UCSC", "NCBI", or NULL if mixed/unknown
get_grna_seq_style <- function(chromosomes) {
  chr <- unique(as.character(chromosomes))
  chr <- chr[!is.na(chr) & chr != ""]
  
  if (!length(chr)) {
    return(NULL)
  }
  
  if (all(grepl("^chr", chr))) {
    return("UCSC")
  }
  
  if (all(!grepl("^chr", chr))) {
    return("NCBI")
  }
  
  return(NULL)
}

# ---- VARIANT CLASSIFICATION FUNCTIONS ---------------------------------------

#' Get Variant Type
#'
#' @param ref Reference allele
#' @param alt Alternative allele
#' @return Variant type (SNP, MNV, DEL, INS, COMPLEX, UNKNOWN)
get_variant_type <- function(ref, alt) {
  if (any(is.na(c(ref, alt))) || ref == "" || alt == "") {
    return("UNKNOWN")
  }
  
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  
  if (ref_len == 1L && alt_len == 1L) return("SNP")
  if (ref_len == alt_len && ref_len > 1L) return("MNV")
  if (ref_len > alt_len) return("DEL")
  if (ref_len < alt_len) return("INS")
  
  return("COMPLEX")
}

#' Classify PAM Variant Impact
#'
#' @param original_pam Original PAM sequence
#' @param alternative_pam Alternative PAM sequence after variant
#' @return Classification (Functional, Reduced_Efficiency, Non_Functional)
classify_pam_variant <- function(original_pam, alternative_pam) {
  alternative_pam <- toupper(alternative_pam)
  
  if (is.null(alternative_pam) || is.na(alternative_pam) || nchar(alternative_pam) != 3L) {
    return("Non_Functional")
  }
  
  # Canonical PAMs (NGG)
  canonical <- c("AGG", "CGG", "GGG", "TGG")
  
  # Non-canonical PAMs (NAG, NGA) - reduced efficiency
  non_canonical <- c("AGA", "CGA", "GGA", "TGA", "AAG", "CAG", "GAG", "TAG")
  
  if (alternative_pam %in% canonical) {
    return("Functional")
  }
  
  if (alternative_pam %in% non_canonical) {
    return("Reduced_Efficiency")
  }
  
  return("Non_Functional")
}

#' Get Variant Coordinates String
#'
#' @param variants GRanges object of variants
#' @return Formatted coordinate string
get_variant_coordinates <- function(variants) {
  if (length(variants) == 0L) {
    return("")
  }
  
  coords <- vapply(seq_along(variants), function(i) {
    variant <- variants[i]
    chr <- as.character(seqnames(variant))
    start_pos <- as.numeric(start(variant))
    end_pos <- as.numeric(end(variant))
    
    ref <- mcols(variant)[["REF"]]
    alt <- mcols(variant)[["ALT"]]
    svtype <- mcols(variant)[["SVTYPE"]]
    vartype <- mcols(variant)[["VARIANT_TYPE"]]
    
    ref <- if (!is.null(ref)) as.character(ref)[1] else NA_character_
    alt <- if (!is.null(alt)) as.character(alt)[1] else NA_character_
    svtype <- if (!is.null(svtype)) as.character(svtype)[1] else NA_character_
    vartype <- if (!is.null(vartype)) as.character(vartype)[1] else NA_character_
    
    # Format based on available information
    if (!is.na(svtype) && svtype != "") {
      if (!is.na(ref) && !is.na(alt) && ref != "" && alt != "") {
        return(sprintf("%s:%s_%s>%s_%s", chr, start_pos, ref, alt, svtype))
      }
      if (!is.na(end_pos) && end_pos != start_pos) {
        return(sprintf("%s:%s-%s_%s", chr, start_pos, end_pos, svtype))
      }
      return(sprintf("%s:%s_%s", chr, start_pos, svtype))
    }
    
    effective_type <- vartype
    if (is.na(effective_type) || effective_type == "") {
      effective_type <- get_variant_type(
        ifelse(is.na(ref), "", ref), 
        ifelse(is.na(alt), "", alt)
      )
    }
    
    if (!is.na(ref) && !is.na(alt) && ref != "" && alt != "") {
      if (!is.na(end_pos) && end_pos != start_pos) {
        return(sprintf("%s:%s-%s_%s>%s_%s", chr, start_pos, end_pos, ref, alt, effective_type))
      }
      return(sprintf("%s:%s_%s>%s_%s", chr, start_pos, ref, alt, effective_type))
    }
    
    if (!is.na(end_pos) && end_pos != start_pos) {
      return(sprintf("%s:%s-%s_%s", chr, start_pos, end_pos, effective_type))
    }
    
    sprintf("%s:%s_%s", chr, start_pos, effective_type)
  }, character(1))
  
  coords <- coords[coords != ""]
  
  if (!length(coords)) {
    return("")
  }
  
  return(paste(unique(coords), collapse = ";"))
}

# ---- gRNA REGION DEFINITION -------------------------------------------------

#' Get gRNA Functional Regions
#'
#' @description Defines PAM, seed, and distal regions for a gRNA
#' @param grna_row Single row from gRNA data.table
#' @return List of GRanges objects for each region
get_grna_regions <- function(grna_row) {
  tryCatch({
    # Validate required columns
    required_cols <- c("chr", "start", "end", "strand", "protospacer", "pam", "cut_site")
    missing_cols <- required_cols[!required_cols %in% names(grna_row)]
    
    if (length(missing_cols)) {
      stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Extract and validate data
    chr <- as.character(grna_row$chr)
    start_pos <- as.numeric(grna_row$start)
    end_pos <- as.numeric(grna_row$end)
    strand <- as.character(grna_row$strand)
    protospacer <- toupper(as.character(grna_row$protospacer))
    pam <- toupper(as.character(grna_row$pam))
    cut_site <- as.numeric(grna_row$cut_site)
    
    if (any(is.na(c(start_pos, end_pos, cut_site, chr))) ||
        is.na(strand) || is.na(protospacer) || is.na(pam)) {
      stop("NA values in gRNA definition")
    }
    
    # Ensure start < end
    if (start_pos > end_pos) {
      tmp <- start_pos
      start_pos <- end_pos
      end_pos <- tmp
    }
    
    # Define protospacer region
    protospacer_region <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = start_pos, end = end_pos),
      strand = strand
    )
    
    # Define regions based on strand
    if (strand == "+") {
      # Plus strand: PAM is downstream
      pam_start <- end_pos + 1L
      pam_end <- end_pos + 3L
      seed_start <- max(start_pos, end_pos - 11L)
      seed_end <- end_pos
      distal_start <- start_pos
      distal_end <- seed_start - 1L
      genomic_protospacer <- protospacer
      genomic_pam <- pam
      
    } else if (strand == "-") {
      # Minus strand: PAM is upstream
      pam_start <- max(1L, start_pos - 3L)
      pam_end <- start_pos - 1L
      seed_start <- start_pos
      seed_end <- min(end_pos, start_pos + 11L)
      distal_start <- seed_end + 1L
      distal_end <- end_pos
      genomic_protospacer <- get_reverse_complement(protospacer)
      genomic_pam <- get_reverse_complement(pam)
      
    } else {
      stop(sprintf("Invalid strand: %s", strand))
    }
    
    # Create region GRanges
    pam_region <- GRanges(
      seqnames = chr, 
      ranges = IRanges(start = pam_start, end = pam_end), 
      strand = strand
    )
    
    seed_region <- GRanges(
      seqnames = chr, 
      ranges = IRanges(start = seed_start, end = seed_end), 
      strand = strand
    )
    
    distal_region <- GRanges(
      seqnames = chr, 
      ranges = IRanges(start = distal_start, end = distal_end), 
      strand = strand
    )
    
    # Define search region (expanded for variant detection)
    search_start <- min(c(start_pos, start(pam_region)), na.rm = TRUE)
    search_end <- max(c(end_pos, end(pam_region)), na.rm = TRUE)
    search_region <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = max(1L, search_start - 2L), end = search_end + 2L),
      strand = "*"
    )
    
    return(list(
      pam_region = pam_region,
      seed_region = seed_region,
      distal_region = distal_region,
      protospacer_genomic = genomic_protospacer,
      pam_genomic = genomic_pam,
      cut_site = cut_site,
      full_protospacer_region = protospacer_region,
      search_region = search_region,
      chromosome = chr,
      strand = strand
    ))
    
  }, error = function(e) {
    stop(sprintf("Error in get_grna_regions: %s", e$message))
  })
}

# ---- VCF LOADING FUNCTIONS --------------------------------------------------

#' Extract Genotype from VCF
#'
#' @param vcf VCF object
#' @return Character vector of genotypes
extract_genotype <- function(vcf) {
  tryCatch({
    if (!"GT" %in% names(geno(vcf))) {
      return(rep("./.", length(vcf)))
    }
    
    gt <- geno(vcf)$GT
    
    if (is.null(dim(gt))) {
      return(as.vector(gt))
    }
    
    if (ncol(gt) == 0L) {
      return(rep("./.", nrow(gt)))
    }
    
    return(gt[, 1L])
    
  }, error = function(e) {
    return(rep("./.", length(vcf)))
  })
}

#' Load SNP/Indel VCF File
#'
#' @param vcf_file Path to VCF file
#' @return GRanges object with variants
load_snp_vcf <- function(vcf_file) {
  
  tryCatch({
    vcf <- readVcf(vcf_file)
    vcf_expanded <- VariantAnnotation::expand(vcf)
    gt <- extract_genotype(vcf_expanded)
    
    # Extract alternative alleles
    alt_sequences <- alt(vcf_expanded)
    if (is(alt_sequences, "DNAStringSetList")) {
      alt_chars <- as.character(unlist(alt_sequences))
    } else {
      alt_chars <- as.character(alt_sequences)
    }
    
    # Create data.table
    vcf_dt <- data.table(
      CHROM = as.character(seqnames(rowRanges(vcf_expanded))),
      POS = start(rowRanges(vcf_expanded)),
      REF = as.character(ref(vcf_expanded)),
      ALT = toupper(alt_chars),
      FILTER = as.character(filt(vcf_expanded)),
      GT = gt
    )
    
    vcf_dt[, VARIANT_TYPE := mapply(get_variant_type, REF, ALT)]
    
    # Convert to GRanges
    gr <- GRanges(
      seqnames = vcf_dt$CHROM,
      ranges = IRanges(
        start = vcf_dt$POS,
        end = vcf_dt$POS + pmax(nchar(vcf_dt$REF), 1L) - 1L
      ),
      REF = vcf_dt$REF,
      ALT = vcf_dt$ALT,
      VARIANT_TYPE = vcf_dt$VARIANT_TYPE,
      FILTER = vcf_dt$FILTER,
      GT = vcf_dt$GT
    )
    
    safe_log("INFO", sprintf("Loaded %s SNP/indel variants", length(gr)))
    return(gr)
    
  }, error = function(e) {
    safe_log("ERROR", sprintf("Error loading SNP VCF: %s", e$message))
    return(GRanges())
  })
}

#' Load Structural Variant VCF File
#'
#' @param vcf_file Path to VCF file
#' @return GRanges object with structural variants
load_sv_vcf <- function(vcf_file) {
  
  tryCatch({
    vcf <- readVcf(vcf_file)
    info_df <- info(vcf)
    
    chrom <- as.character(seqnames(rowRanges(vcf)))
    pos <- start(rowRanges(vcf))
    filter_vals <- as.character(filt(vcf))
    gt <- extract_genotype(vcf)
    
    # Extract SV-specific fields
    end_vals <- if ("END" %in% names(info_df)) {
      as.numeric(info_df$END)
    } else {
      rep(NA_real_, length(vcf))
    }
    
    svtype_vals <- if ("SVTYPE" %in% names(info_df)) {
      as.character(info_df$SVTYPE)
    } else {
      rep("UNKNOWN", length(vcf))
    }
    
    svlen_vals <- if ("SVLEN" %in% names(info_df)) {
      as.numeric(info_df$SVLEN)
    } else {
      rep(NA_real_, length(vcf))
    }
    
    # Create data.table
    sv_dt <- data.table(
      CHROM = chrom,
      POS = pos,
      END = end_vals,
      SVTYPE = svtype_vals,
      SVLEN = svlen_vals,
      FILTER = filter_vals,
      GT = gt
    )
    
    # Fix missing END values
    sv_dt[is.na(END) | END < POS, END := POS]
    sv_dt[SVTYPE == "BND" & (is.na(END) | END <= POS), END := POS + 1L]
    
    # Convert to GRanges
    gr <- GRanges(
      seqnames = sv_dt$CHROM,
      ranges = IRanges(start = sv_dt$POS, end = pmax(sv_dt$END, sv_dt$POS)),
      SVTYPE = sv_dt$SVTYPE,
      SVLEN = sv_dt$SVLEN,
      FILTER = sv_dt$FILTER,
      GT = sv_dt$GT,
      END = sv_dt$END
    )
    
    safe_log("INFO", sprintf("Loaded %s structural variants", length(gr)))
    return(gr)
    
  }, error = function(e) {
    safe_log("ERROR", sprintf("Error loading SV VCF: %s", e$message))
    return(GRanges())
  })
}

#' Filter Variants by gRNA Regions
#'
#' @description Pre-filters variants to those near gRNA regions
#' @param snp_gr GRanges of SNP/indel variants
#' @param sv_gr GRanges of structural variants
#' @param grna_data data.table of gRNA data
#' @return List with filtered snp_gr and sv_gr
filter_variants_by_grna_regions <- function(snp_gr, sv_gr, grna_data) {
  safe_log("INFO", "Pre-filtering variants by gRNA regions")
  
  if (nrow(grna_data) == 0L) {
    return(list(snp_gr = snp_gr, sv_gr = sv_gr))
  }
  
  # Create extended regions around gRNAs
  extended_start <- pmax(as.numeric(grna_data$start) - 5L, 1L)
  extended_end <- as.numeric(grna_data$end) + 5L
  
  all_regions <- GRanges(
    seqnames = as.character(grna_data$chr),
    ranges = IRanges(start = extended_start, end = extended_end),
    strand = "*"
  )
  
  all_regions <- reduce(all_regions, ignore.strand = TRUE)
  
  # Filter variants
  if (length(snp_gr) > 0L) {
    snp_gr <- subsetByOverlaps(snp_gr, all_regions, ignore.strand = TRUE)
    safe_log("INFO", sprintf("Retained %s SNP/indel variants", length(snp_gr)))
  }
  
  if (length(sv_gr) > 0L) {
    sv_gr <- subsetByOverlaps(sv_gr, all_regions, ignore.strand = TRUE)
    safe_log("INFO", sprintf("Retained %s SV variants", length(sv_gr)))
  }
  
  return(list(snp_gr = snp_gr, sv_gr = sv_gr))
}

#' Prepare Variant Index by Chromosome
#'
#' @param variant_gr GRanges of variants
#' @return List of GRanges split by chromosome
prepare_variant_index <- function(variant_gr) {
  if (length(variant_gr) == 0L) {
    return(list())
  }
  
  chr <- as.character(seqnames(variant_gr))
  valid <- !is.na(chr) & chr != ""
  
  if (!any(valid)) {
    return(list())
  }
  
  variant_gr <- variant_gr[valid]
  chr <- chr[valid]
  
  split_list <- split(variant_gr, chr, drop = TRUE)
  
  out <- vector("list", length(split_list))
  names(out) <- names(split_list)
  
  for (nm in names(split_list)) {
    out[[nm]] <- split_list[[nm]]
  }
  
  return(out)
}

# ---- STRUCTURAL VARIANT CLASSIFICATION --------------------------------------

#' Classify Structural Variant Impact on gRNA
#'
#' @param variant_record GRanges object (single variant)
#' @param grna_regions List of gRNA region GRanges
#' @param variant_type Type of structural variant
#' @return List with impact classification
classify_sv_impact <- function(variant_record, grna_regions, variant_type = NULL) {
  tryCatch({
    if (length(variant_record) > 1L) {
      variant_record <- variant_record[1L]
    }
    
    # Determine variant type
    if (is.null(variant_type) || is.na(variant_type)) {
      if ("SVTYPE" %in% names(mcols(variant_record))) {
        variant_type <- as.character(mcols(variant_record)[["SVTYPE"]])
      } else if ("VARIANT_TYPE" %in% names(mcols(variant_record))) {
        variant_type <- as.character(mcols(variant_record)[["VARIANT_TYPE"]])
      } else {
        variant_type <- "UNKNOWN"
      }
    }
    
    vt <- toupper(variant_type)
    
    # Check overlaps with functional regions
    has_seed <- length(grna_regions$seed_region) > 0 && width(grna_regions$seed_region) > 0
    has_pam <- length(grna_regions$pam_region) > 0 && width(grna_regions$pam_region) > 0
    has_distal <- length(grna_regions$distal_region) > 0 && width(grna_regions$distal_region) > 0
    
    overlaps_seed <- has_seed && overlapsAny(variant_record, grna_regions$seed_region, ignore.strand = TRUE)
    overlaps_pam <- has_pam && overlapsAny(variant_record, grna_regions$pam_region, ignore.strand = TRUE)
    overlaps_distal <- has_distal && overlapsAny(variant_record, grna_regions$distal_region, ignore.strand = TRUE)
    
    # Calculate region boundaries
    seed_pam_ranges_start <- numeric()
    seed_pam_ranges_end <- numeric()
    
    if (has_seed) {
      seed_pam_ranges_start <- c(seed_pam_ranges_start, as.numeric(start(grna_regions$seed_region)))
      seed_pam_ranges_end <- c(seed_pam_ranges_end, as.numeric(end(grna_regions$seed_region)))
    }
    if (has_pam) {
      seed_pam_ranges_start <- c(seed_pam_ranges_start, as.numeric(start(grna_regions$pam_region)))
      seed_pam_ranges_end <- c(seed_pam_ranges_end, as.numeric(end(grna_regions$pam_region)))
    }
    
    seed_pam_start <- if (length(seed_pam_ranges_start)) min(seed_pam_ranges_start) else NA_real_
    seed_pam_end <- if (length(seed_pam_ranges_end)) max(seed_pam_ranges_end) else NA_real_
    
    proto_pam_ranges_start <- c(as.numeric(start(grna_regions$full_protospacer_region)))
    proto_pam_ranges_end <- c(as.numeric(end(grna_regions$full_protospacer_region)))
    
    if (has_pam) {
      proto_pam_ranges_start <- c(proto_pam_ranges_start, as.numeric(start(grna_regions$pam_region)))
      proto_pam_ranges_end <- c(proto_pam_ranges_end, as.numeric(end(grna_regions$pam_region)))
    }
    
    proto_pam_start <- if (length(proto_pam_ranges_start)) min(proto_pam_ranges_start) else NA_real_
    proto_pam_end <- if (length(proto_pam_ranges_end)) max(proto_pam_ranges_end) else NA_real_
    
    # Get variant boundaries
    variant_start <- as.numeric(start(variant_record))
    variant_end <- as.numeric(end(variant_record))
    
    if (!is.na(variant_start) && !is.na(variant_end) && variant_start > variant_end) {
      tmp <- variant_start
      variant_start <- variant_end
      variant_end <- tmp
    }
    
    # Check if variant contains critical regions
    contains_seed_pam <- !is.na(seed_pam_start) && !is.na(seed_pam_end) &&
      !is.na(variant_start) && !is.na(variant_end) &&
      variant_start <= seed_pam_start && variant_end >= seed_pam_end
    
    contains_proto_pam <- !is.na(proto_pam_start) && !is.na(proto_pam_end) &&
      !is.na(variant_start) && !is.na(variant_end) &&
      variant_start <= proto_pam_start && variant_end >= proto_pam_end
    
    # Classify impact
    impact <- "None"
    reason <- "No_effect"
    affects_seed <- FALSE
    affects_pam <- FALSE
    affects_distal <- FALSE
    
    if (vt %in% c("DEL", "INS", "DUP", "INV", "TRA", "BND", "CNV", "UNK", "UNKNOWN")) {
      
      if (contains_proto_pam && vt %in% c("DUP", "INV")) {
        # Duplication/inversion containing entire gRNA - uncertain effect
        impact <- "None"
        reason <- sprintf("%s_contains_protospacer_to_PAM", vt)
        
      } else if (contains_seed_pam && vt %in% c("DUP", "INV")) {
        # Duplication/inversion of seed+PAM - likely reduced efficiency
        impact <- "Reduced"
        reason <- sprintf("%s_contains_seed_to_PAM", vt)
        affects_seed <- has_seed
        affects_pam <- has_pam
        
      } else if (contains_seed_pam && vt == "DEL") {
        # Deletion of seed+PAM - non-functional
        impact <- "Nonfunctional"
        reason <- "DEL_contains_seed_to_PAM"
        affects_seed <- has_seed
        affects_pam <- has_pam
        
      } else {
        # Check specific overlaps
        if (overlaps_seed || overlaps_pam) {
          impact <- "Nonfunctional"
          reason <- sprintf("%s_overlap_seed_or_PAM", vt)
          affects_seed <- overlaps_seed
          affects_pam <- overlaps_pam
          
        } else if (overlaps_distal) {
          impact <- "Reduced"
          reason <- sprintf("%s_overlap_distal", vt)
          affects_distal <- TRUE
          
        } else {
          impact <- "None"
          reason <- sprintf("%s_no_relevant_overlap", vt)
        }
      }
    }
    
    return(list(
      impact = impact,
      reason = reason,
      affects_seed = affects_seed,
      affects_pam = affects_pam,
      affects_distal = affects_distal,
      contains_seed_to_pam = contains_seed_pam,
      contains_protospacer_to_pam = contains_proto_pam
    ))
    
  }, error = function(e) {
    return(list(
      impact = "None",
      reason = sprintf("Error_in_classification:%s", e$message),
      affects_seed = FALSE,
      affects_pam = FALSE,
      affects_distal = FALSE,
      contains_seed_to_pam = FALSE,
      contains_protospacer_to_pam = FALSE
    ))
  })
}

# ---- gRNA CLASSIFICATION FUNCTION -------------------------------------------

#' Classify Single gRNA Based on Variants
#'
#' @description Main classification function for a single gRNA
#' @param grna_row Single row from gRNA data.table
#' @param snp_by_chr List of SNP/indel GRanges by chromosome
#' @param sv_by_chr List of SV GRanges by chromosome
#' @param process_id Process ID for logging
#' @return data.table with classification results
classify_grna_optimized <- function(grna_row, snp_by_chr, sv_by_chr, process_id = NULL) {
  tryCatch({
    # Convert to list for easier access
    grna_list <- as.list(grna_row)
    
    # Get gRNA regions
    regions <- get_grna_regions(grna_list)
    chr_key <- regions$chromosome
    
    # Get candidate variants for this chromosome
    snp_candidates <- if (!is.null(snp_by_chr[[chr_key]])) {
      subsetByOverlaps(snp_by_chr[[chr_key]], regions$search_region, ignore.strand = TRUE)
    } else {
      GRanges()
    }
    
    sv_candidates <- if (!is.null(sv_by_chr[[chr_key]])) {
      subsetByOverlaps(sv_by_chr[[chr_key]], regions$search_region, ignore.strand = TRUE)
    } else {
      GRanges()
    }
    
    # Initialize result
    original_pam_grna <- to_grna_orientation(regions$pam_genomic, regions$strand)
    original_protospacer_grna <- to_grna_orientation(regions$protospacer_genomic, regions$strand)
    
    result <- list(
      gRNA_id = ifelse(is.null(grna_list$ID) || is.na(grna_list$ID), "UNKNOWN", as.character(grna_list$ID)),
      chr = chr_key,
      category = "Fully_Functional",
      score = 100L,
      pam_status = "Canonical",
      seed_status = "Intact",
      sv_impact = "None",
      snp_indel_impact = "None",
      details = "",
      original_pam = original_pam_grna,
      altered_pam = original_pam_grna,
      original_protospacer = original_protospacer_grna,
      altered_protospacer = original_protospacer_grna,
      protospacer_length = nchar(original_protospacer_grna),
      seed_variant_coordinates = "",
      distal_variant_coordinates = "",
      pam_variant_coordinates = "",
      sv_variant_coordinates = "",
      strand = regions$strand,
      grna_coordinates = sprintf("%s:%s-%s", chr_key, grna_list$start, grna_list$end)
    )
    
    # Classification flags
    nonfunctional <- FALSE
    reduced <- FALSE
    snp_impacts <- character()
    sv_impacts <- character()
    
    # Current sequences (will be modified by variants)
    current_pam_genomic <- toupper(regions$pam_genomic)
    current_protospacer_genomic <- toupper(regions$protospacer_genomic)
    
    # Check original PAM
    original_pam_class <- classify_pam_variant(original_pam_grna, original_pam_grna)
    if (original_pam_class == "Functional") {
      result$pam_status <- "Canonical"
    } else if (original_pam_class == "Reduced_Efficiency") {
      result$pam_status <- "Non_canonical"
      reduced <- TRUE
    } else {
      result$pam_status <- "Disrupted"
      nonfunctional <- TRUE
    }
    
    # Process PAM variants
    if (length(regions$pam_region) > 0L && length(snp_candidates) > 0L) {
      pam_variants <- subsetByOverlaps(snp_candidates, regions$pam_region, ignore.strand = TRUE)
      
      if (length(pam_variants) > 0L) {
        pam_variants <- pam_variants[order(start(pam_variants))]
        result$pam_variant_coordinates <- get_variant_coordinates(pam_variants)
        pam_nonfunctional <- FALSE
        
        for (j in seq_along(pam_variants)) {
          variant <- pam_variants[j]
          ref <- as.character(mcols(variant)[["REF"]])
          alt <- as.character(mcols(variant)[["ALT"]])
          vartype <- mcols(variant)[["VARIANT_TYPE"]]
          
          if (is.null(vartype) || is.na(vartype)) {
            vartype <- get_variant_type(ref, alt)
          }
          vartype <- toupper(vartype)
          
          snp_impacts <- c(snp_impacts, sprintf("PAM_%s", vartype))
          
          if (is.null(current_pam_genomic) || is.na(current_pam_genomic) || current_pam_genomic == "") {
            pam_nonfunctional <- TRUE
            break
          }
          
          if (vartype %in% c("SNP", "SNV", "MNV", "MNP")) {
            new_pam <- apply_variant_to_sequence(
              current_pam_genomic,
              start(variant),
              ref,
              alt,
              start(regions$pam_region)
            )
            
            if (is.null(new_pam) || is.na(new_pam) || nchar(new_pam) != 3L) {
              pam_nonfunctional <- TRUE
              break
            }
            
            current_pam_genomic <- toupper(new_pam)
          } else {
            # Indel in PAM - non-functional
            pam_nonfunctional <- TRUE
            break
          }
        }
        
        if (pam_nonfunctional) {
          result$pam_status <- "Disrupted"
          current_pam_genomic <- NA_character_
          nonfunctional <- TRUE
        } else {
          pam_class <- classify_pam_variant(
            original_pam_grna,
            to_grna_orientation(current_pam_genomic, regions$strand)
          )
          
          if (pam_class == "Functional") {
            result$pam_status <- "Canonical"
          } else if (pam_class == "Reduced_Efficiency") {
            result$pam_status <- "Non_canonical"
            if (!nonfunctional) reduced <- TRUE
          } else {
            result$pam_status <- "Disrupted"
            nonfunctional <- TRUE
          }
        }
      }
    }
    
    # Process seed region variants
    if (length(regions$seed_region) > 0L && length(snp_candidates) > 0L) {
      seed_variants <- subsetByOverlaps(snp_candidates, regions$seed_region, ignore.strand = TRUE)
      
      if (length(seed_variants) > 0L) {
        result$seed_variant_coordinates <- get_variant_coordinates(seed_variants)
        result$seed_status <- "Disrupted"
        
        for (j in seq_along(seed_variants)) {
          variant <- seed_variants[j]
          ref <- as.character(mcols(variant)[["REF"]])
          alt <- as.character(mcols(variant)[["ALT"]])
          vartype <- mcols(variant)[["VARIANT_TYPE"]]
          
          if (is.null(vartype) || is.na(vartype)) {
            vartype <- get_variant_type(ref, alt)
          }
          vartype <- toupper(vartype)
          
          snp_impacts <- c(snp_impacts, sprintf("Seed_%s", vartype))
          
          # Apply SNP/MNP to protospacer
          if (vartype %in% c("SNP", "SNV", "MNV", "MNP")) {
            current_protospacer_genomic <- apply_variant_to_sequence(
              current_protospacer_genomic,
              start(variant),
              ref,
              alt,
              start(regions$full_protospacer_region)
            )
          }
        }
        
        nonfunctional <- TRUE
      }
    }
    
    # Process distal region variants
    if (length(regions$distal_region) > 0L && length(snp_candidates) > 0L) {
      distal_variants <- subsetByOverlaps(snp_candidates, regions$distal_region, ignore.strand = TRUE)
      
      if (length(distal_variants) > 0L) {
        result$distal_variant_coordinates <- get_variant_coordinates(distal_variants)
        
        for (j in seq_along(distal_variants)) {
          variant <- distal_variants[j]
          ref <- as.character(mcols(variant)[["REF"]])
          alt <- as.character(mcols(variant)[["ALT"]])
          vartype <- mcols(variant)[["VARIANT_TYPE"]]
          
          if (is.null(vartype) || is.na(vartype)) {
            vartype <- get_variant_type(ref, alt)
          }
          vartype <- toupper(vartype)
          
          snp_impacts <- c(snp_impacts, sprintf("Distal_%s", vartype))
          
          # Apply SNP/MNP to protospacer
          if (vartype %in% c("SNP", "SNV", "MNV", "MNP")) {
            current_protospacer_genomic <- apply_variant_to_sequence(
              current_protospacer_genomic,
              start(variant),
              ref,
              alt,
              start(regions$full_protospacer_region)
            )
          } else {
            # Indel in distal region - reduced efficiency
            if (!nonfunctional) reduced <- TRUE
          }
        }
        
        if (!nonfunctional) reduced <- TRUE
      }
    }
    
    # Process structural variants
    if (length(sv_candidates) > 0L) {
      overlapping_sv <- unique(sv_candidates)
      
      if (length(overlapping_sv) > 0L) {
        result$sv_variant_coordinates <- get_variant_coordinates(overlapping_sv)
        
        for (j in seq_along(overlapping_sv)) {
          variant <- overlapping_sv[j]
          svtype <- if ("SVTYPE" %in% names(mcols(variant))) {
            as.character(mcols(variant)[["SVTYPE"]])
          } else if ("VARIANT_TYPE" %in% names(mcols(variant))) {
            as.character(mcols(variant)[["VARIANT_TYPE"]])
          } else {
            "UNKNOWN"
          }
          
          sv_result <- classify_sv_impact(variant, regions, svtype)
          sv_impacts <- c(sv_impacts, sv_result$reason)
          
          if (sv_result$impact == "Nonfunctional") {
            nonfunctional <- TRUE
            result$sv_impact <- "Nonfunctional"
            
            if (sv_result$affects_seed || sv_result$contains_seed_to_pam) {
              result$seed_status <- "Disrupted"
            }
            if (sv_result$affects_pam || sv_result$contains_seed_to_pam) {
              result$pam_status <- "Disrupted"
            }
            
          } else if (sv_result$impact == "Reduced") {
            if (result$sv_impact != "Nonfunctional") {
              result$sv_impact <- "Reduced"
            }
            if (!nonfunctional) reduced <- TRUE
            
            if (sv_result$affects_seed) {
              result$seed_status <- "Disrupted"
            }
            if (sv_result$affects_pam && result$pam_status != "Disrupted") {
              result$pam_status <- "Non_canonical"
            }
          }
        }
      }
    }
    
    # Summarize SNP/indel impact
    if (!length(snp_impacts)) {
      result$snp_indel_impact <- "None"
    } else {
      result$snp_indel_impact <- paste(sort(unique(snp_impacts)), collapse = ";")
    }
    
    # Combine all details
    all_details <- unique(c(snp_impacts, sv_impacts))
    result$details <- if (length(all_details)) {
      paste(sort(unique(all_details)), collapse = ";")
    } else {
      ""
    }
    
    # Set altered sequences
    if (!is.null(current_pam_genomic) && !is.na(current_pam_genomic) && current_pam_genomic != "") {
      result$altered_pam <- to_grna_orientation(current_pam_genomic, regions$strand)
    } else {
      result$altered_pam <- NA_character_
    }
    
    if (!is.null(current_protospacer_genomic) && !is.na(current_protospacer_genomic) && 
        current_protospacer_genomic != "") {
      result$altered_protospacer <- to_grna_orientation(current_protospacer_genomic, regions$strand)
    } else {
      result$altered_protospacer <- NA_character_
    }
    
    # Final classification
    if (nonfunctional) {
      result$category <- "Critical_Failure"
      result$score <- 0L
    } else if (reduced) {
      result$category <- "Reduced_Efficiency"
      result$score <- 50L
    } else {
      result$category <- "Fully_Functional"
      result$score <- 100L
    }
    
    return(as.data.table(result))
    
  }, error = function(e) {
    error_msg <- sprintf(
      "Error in gRNA %s: %s",
      ifelse(is.null(grna_row$ID), "UNKNOWN", as.character(grna_row$ID)),
      e$message
    )
    safe_log("ERROR", error_msg, process_id)
    
    return(data.table(
      gRNA_id = ifelse(is.null(grna_row$ID) || is.na(grna_row$ID), "UNKNOWN", as.character(grna_row$ID)),
      chr = ifelse(is.null(grna_row$chr), "Error", as.character(grna_row$chr)),
      category = "Processing_Error",
      score = 0L,
      pam_status = "Error",
      seed_status = "Error",
      sv_impact = "Error",
      snp_indel_impact = "Error",
      details = sprintf("Error: %s", e$message),
      original_pam = "Error",
      altered_pam = "Error",
      original_protospacer = "Error",
      altered_protospacer = "Error",
      protospacer_length = -1L,
      seed_variant_coordinates = "Error",
      distal_variant_coordinates = "Error",
      pam_variant_coordinates = "Error",
      sv_variant_coordinates = "Error",
      strand = ifelse(is.null(grna_row$strand), "Error", as.character(grna_row$strand)),
      grna_coordinates = ifelse(
        is.null(grna_row$chr) || is.null(grna_row$start) || is.null(grna_row$end),
        "Error",
        sprintf("%s:%s-%s", grna_row$chr, grna_row$start, grna_row$end)
      )
    ))
  })
}

# ---- PARALLEL PROCESSING FUNCTIONS ------------------------------------------

#' Create Static Batch Indices
#'
#' @param n Total number of items
#' @param workers Number of workers
#' @return List of index vectors
static_batch_indices <- function(n, workers) {
  batch_size <- ceiling(n / workers)
  split(seq_len(n), ceiling(seq_along(seq_len(n)) / batch_size))
}

#' Process Batch of gRNAs
#'
#' @description Worker function for parallel processing
#' @param batch_indices Indices of gRNAs to process
#' @param grna_data Full gRNA data.table
#' @param snp_by_chr List of SNP GRanges by chromosome
#' @param sv_by_chr List of SV GRanges by chromosome
#' @param progress_step Progress logging interval
#' @return data.table with classification results
process_grna_batch_optimized <- function(
    batch_indices, 
    grna_data, 
    snp_by_chr,
    sv_by_chr, 
    progress_step
) {
  
  process_id <- Sys.getpid()
  safe_log("INFO", sprintf("Worker %s processing %s gRNAs", process_id, length(batch_indices)), process_id)
  
  # Load required packages in worker
  local_packages <- c("GenomicRanges", "IRanges", "data.table", "Biostrings")
  invisible(lapply(local_packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
  
  # Garbage collection for large batches
  if (length(batch_indices) > 500L) {
    gc(full = TRUE)
  }
  
  # Process each gRNA
  results <- vector("list", length(batch_indices))
  
  for (i in seq_along(batch_indices)) {
    idx <- batch_indices[[i]]
    grna_row <- grna_data[idx]
    results[[i]] <- classify_grna_optimized(grna_row, snp_by_chr, sv_by_chr, process_id)
    
    # Progress logging
    if (progress_step > 0L && (i %% progress_step == 0L)) {
      progress <- (i / length(batch_indices)) * 100
      safe_log("DEBUG", sprintf("Progress: %.2f%%", progress), process_id)
    }
  }
  
  safe_log("INFO", sprintf("Worker %s completed", process_id), process_id)
  
  return(rbindlist(results, use.names = TRUE, fill = TRUE))
}

# ---- MAIN ANALYSIS FUNCTION -------------------------------------------------

#' Main Analysis Function
#'
#' @description Executes complete gRNA variant impact analysis
#' @return NULL (results saved to files)
main <- function() {
  safe_log("INFO", strrep("=", 70))
  safe_log("INFO", "STARTING gRNA VARIANT IMPACT ANALYSIS")
  safe_log("INFO", strrep("=", 70))
  
  tryCatch({
    # Load gRNA data
    safe_log("INFO", sprintf("Loading gRNA data: %s", opt$grna))
    
    if (!file.exists(opt$grna)) {
      stop(sprintf("gRNA file not found: %s", opt$grna))
    }
    
    grna_data <- fread(opt$grna, showProgress = FALSE)
    safe_log("INFO", sprintf("Loaded %s gRNAs", nrow(grna_data)))
    
    # Validate required columns
    required_cols <- c("ID", "chr", "start", "end", "strand", "protospacer", "pam", "cut_site")
    missing_cols <- setdiff(required_cols, names(grna_data))
    
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Data type conversion
    grna_data[, ID := as.character(ID)]
    grna_data[, chr := as.character(chr)]
    grna_data[, strand := as.character(strand)]
    grna_data[, start := as.numeric(start)]
    grna_data[, end := as.numeric(end)]
    grna_data[, cut_site := as.numeric(cut_site)]
    grna_data[, protospacer := toupper(as.character(protospacer))]
    grna_data[, pam := toupper(as.character(pam))]
    
    # Data cleaning
    initial_total <- nrow(grna_data)
    
    grna_data <- grna_data[!is.na(chr) & chr != ""]
    removed_chr <- initial_total - nrow(grna_data)
    if (removed_chr > 0L) {
      safe_log("WARN", sprintf("Removed %s gRNAs with missing chromosome", removed_chr))
    }
    
    grna_data <- grna_data[!is.na(start) & !is.na(end) & start > 0 & end > 0]
    
    # Fix inverted coordinates
    invalid_range <- grna_data[start > end]
    if (nrow(invalid_range) > 0L) {
      safe_log("WARN", sprintf("Fixed %s gRNAs with start > end", nrow(invalid_range)))
      grna_data[start > end, c("start", "end") := .(end, start)]
    }
    
    grna_data <- grna_data[!is.na(cut_site)]
    grna_data <- grna_data[!is.na(protospacer) & protospacer != ""]
    grna_data <- grna_data[!is.na(pam) & pam != ""]
    
    invalid_strand <- grna_data[!strand %in% c("+", "-")]
    if (nrow(invalid_strand) > 0L) {
      safe_log("WARN", sprintf("Removed %s gRNAs with invalid strand", nrow(invalid_strand)))
      grna_data <- grna_data[strand %in% c("+", "-")]
    }
    
    if (nrow(grna_data) == 0L) {
      stop("No valid gRNA records after filtering")
    }
    
    safe_log("INFO", sprintf("After cleaning: %s valid gRNAs", nrow(grna_data)))
    
    # Determine chromosome style
    grna_style <- get_grna_seq_style(grna_data$chr)
    if (is.null(grna_style)) {
      safe_log("WARN", "Could not determine chromosome style; will attempt harmonization with VCF")
    } else {
      safe_log("INFO", sprintf("gRNA chromosome style: %s", grna_style))
    }
    
    # Load variant files
    safe_log("INFO", sprintf("Loading SNP/indel VCF: %s", opt$snpvcf))
    if (!file.exists(opt$snpvcf)) {
      stop(sprintf("SNP VCF not found: %s", opt$snpvcf))
    }
    snp_gr <- load_snp_vcf(opt$snpvcf)
    
    safe_log("INFO", sprintf("Loading SV VCF: %s", opt$svvcf))
    if (!file.exists(opt$svvcf)) {
      stop(sprintf("SV VCF not found: %s", opt$svvcf))
    }
    sv_gr <- load_sv_vcf(opt$svvcf)
    
    # Harmonize chromosome styles
    variant_styles <- unique(c(
      if (length(snp_gr) > 0L) seqlevelsStyle(snp_gr)[1] else NA_character_,
      if (length(sv_gr) > 0L) seqlevelsStyle(sv_gr)[1] else NA_character_
    ))
    variant_styles <- variant_styles[!is.na(variant_styles)]
    
    if (is.null(grna_style) && length(variant_styles) > 0L) {
      target_style <- variant_styles[1]
      safe_log("INFO", sprintf("Converting gRNA chromosomes to %s style", target_style))
      grna_data[, chr := convert_chr_style(chr, target_style)]
      grna_style <- get_grna_seq_style(grna_data$chr)
    } else {
      target_style <- grna_style
    }
    
    if (!is.null(target_style)) {
      if (length(snp_gr) > 0L) {
        tryCatch(
          { seqlevelsStyle(snp_gr) <- target_style },
          error = function(e) safe_log("WARN", sprintf("Could not set SNP VCF style: %s", e$message))
        )
      }
      if (length(sv_gr) > 0L) {
        tryCatch(
          { seqlevelsStyle(sv_gr) <- target_style },
          error = function(e) safe_log("WARN", sprintf("Could not set SV VCF style: %s", e$message))
        )
      }
    }
    
    # Filter variants by gRNA regions
    filtered_variants <- filter_variants_by_grna_regions(snp_gr, sv_gr, grna_data)
    snp_gr <- filtered_variants$snp_gr
    sv_gr <- filtered_variants$sv_gr
    
    # Prepare variant indices
    snp_by_chr <- prepare_variant_index(snp_gr)
    sv_by_chr <- prepare_variant_index(sv_gr)
    
    # Prepare batches for parallel processing
    batches <- static_batch_indices(nrow(grna_data), max(1L, opt$workers))
    num_batches <- length(batches)
    num_workers <- min(opt$workers, max(1L, num_batches))
    
    safe_log("INFO", sprintf("Starting parallel processing: %s batches, %s workers", num_batches, num_workers))
    
    # Configure parallel backend
    if (.Platform$OS.type == "unix") {
      bp_param <- MulticoreParam(
        workers = num_workers,
        progressbar = TRUE,
        stop.on.error = FALSE
      )
    } else {
      bp_param <- SnowParam(
        workers = num_workers,
        progressbar = TRUE,
        stop.on.error = FALSE
      )
    }
    
    # Run parallel processing
    start_time <- Sys.time()
    
    results_list <- bplapply(
      batches,
      function(batch_indices) {
        process_grna_batch_optimized(batch_indices, grna_data, snp_by_chr, sv_by_chr, progress_step)
      },
      BPPARAM = bp_param
    )
    
    # Combine results
    all_results <- if (length(results_list) > 0L) {
      rbindlist(results_list, use.names = TRUE, fill = TRUE)
    } else {
      data.table()
    }
    
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    safe_log("INFO", sprintf("Processing completed in %.2f seconds", processing_time))
    
    # Separate successful and error results
    if (!"category" %in% names(all_results)) {
      all_results[, category := "Processing_Error"]
    }
    
    success_results <- all_results[category != "Processing_Error"]
    error_results <- all_results[category == "Processing_Error"]
    
    # Save successful results
    if (nrow(success_results) > 0L) {
      output_file <- file.path(opt$output, "gRNA_classification_final.csv")
      fwrite(success_results, output_file)
      safe_log("INFO", sprintf("Saved %s successful results: %s", nrow(success_results), output_file))
      
      # Generate statistics
      stats_file <- file.path(opt$output, "analysis_statistics.txt")
      stats_text <- capture.output({
        cat(strrep("=", 70), "\n")
        cat("gRNA VARIANT IMPACT ANALYSIS - STATISTICS\n")
        cat(strrep("=", 70), "\n\n")
        cat("Analysis Date:", as.character(Sys.time()), "\n")
        cat("Processing Time:", sprintf("%.2f seconds", processing_time), "\n\n")
        cat("INPUT DATA:\n")
        cat("  Total gRNAs:", scales::comma(nrow(grna_data)), "\n")
        cat("  Successfully processed:", scales::comma(nrow(success_results)), "\n")
        cat("  Errors:", scales::comma(nrow(error_results)), "\n\n")
        
        if ("chr" %in% names(success_results)) {
          cat("DISTRIBUTION BY CHROMOSOME:\n")
          print(table(success_results$chr))
          cat("\n")
        }
        
        if ("category" %in% names(success_results)) {
          cat("FUNCTIONAL CATEGORIES:\n")
          cat_table <- table(success_results$category)
          for (cat_name in names(cat_table)) {
            pct <- cat_table[cat_name] / nrow(success_results) * 100
            cat(sprintf("  %-25s: %10s (%5.2f%%)\n", 
                        cat_name, 
                        scales::comma(cat_table[cat_name]), 
                        pct))
          }
          cat("\n")
        }
        
        if ("pam_status" %in% names(success_results)) {
          cat("PAM STATUS:\n")
          pam_table <- table(success_results$pam_status)
          for (pam_name in names(pam_table)) {
            pct <- pam_table[pam_name] / nrow(success_results) * 100
            cat(sprintf("  %-25s: %10s (%5.2f%%)\n", 
                        pam_name, 
                        scales::comma(pam_table[pam_name]), 
                        pct))
          }
          cat("\n")
        }
        
        if ("seed_status" %in% names(success_results)) {
          cat("SEED REGION STATUS:\n")
          seed_table <- table(success_results$seed_status)
          for (seed_name in names(seed_table)) {
            pct <- seed_table[seed_name] / nrow(success_results) * 100
            cat(sprintf("  %-25s: %10s (%5.2f%%)\n", 
                        seed_name, 
                        scales::comma(seed_table[seed_name]), 
                        pct))
          }
          cat("\n")
        }
        
        if ("sv_impact" %in% names(success_results)) {
          cat("STRUCTURAL VARIANT IMPACT:\n")
          sv_table <- table(success_results$sv_impact)
          for (sv_name in names(sv_table)) {
            pct <- sv_table[sv_name] / nrow(success_results) * 100
            cat(sprintf("  %-25s: %10s (%5.2f%%)\n", 
                        sv_name, 
                        scales::comma(sv_table[sv_name]), 
                        pct))
          }
          cat("\n")
        }
        
        cat(strrep("=", 70), "\n")
      })
      
      writeLines(stats_text, stats_file)
      safe_log("INFO", sprintf("Saved statistics: %s", stats_file))
      
    } else {
      safe_log("WARN", "No successful results to save")
    }
    
    # Save error results
    if (nrow(error_results) > 0L) {
      error_file <- file.path(opt$output, "gRNA_processing_errors.csv")
      fwrite(error_results, error_file)
      safe_log("WARN", sprintf("Saved %s errors: %s", nrow(error_results), error_file))
    } else {
      safe_log("INFO", "No processing errors")
    }
    
    safe_log("INFO", strrep("=", 70))
    safe_log("INFO", "ANALYSIS COMPLETED SUCCESSFULLY")
    safe_log("INFO", strrep("=", 70))
    
  }, error = function(e) {
    safe_log("ERROR", sprintf("CRITICAL ERROR: %s", e$message))
    traceback(1)
    quit(save = "no", status = 1)
  })
}

# ---- EXECUTE MAIN FUNCTION --------------------------------------------------

main()

# =============================================================================
# END OF SCRIPT
# =============================================================================
