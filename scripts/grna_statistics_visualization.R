#!/usr/bin/env Rscript
# =============================================================================
# gRNA Functionality Statistics & Visualization
# =============================================================================
# Description: Analysis and visualization of CRISPR gRNA functionality
#              based on classification by genomic variants 
# Author: Denis Omelchenko
# Date: 2025
# Version: 1.0.0
# License: MIT
# =============================================================================

# ---- PACKAGE MANAGEMENT -----------------------------------------------------

#' Load Required Packages
#'
#' @description Installs and loads all required packages
#' @param packages Character vector of package names
#' @return NULL (packages loaded as side effect)
load_required_packages <- function(packages) {
  missing_pkgs <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Missing required packages: %s\nPlease install them first:\ninstall.packages(c(%s))",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
    ))
  }
  
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

REQUIRED_PACKAGES <- c(
  "optparse", "data.table", "ggplot2", "scales",
  "forcats", "splitstackshape", "patchwork", "dplyr", "tidyr",
  "purrr", "stringr", "ggpubr"
)

load_required_packages(REQUIRED_PACKAGES)

# ---- GLOBAL SETTINGS AND PALETTES -------------------------------------------

# Color palettes (matching original)
PALETTE_CONSISTENT <- c(
  "#2E8B57",  # Green
  "#F4A460",  # Orange
  "#A50026"   # Red
)

PALETTE_CATEGORIES <- c(
  Fully_Functional   = "#2E8B57",
  Reduced_Efficiency = "#F4A460",
  Critical_Failure   = "#A50026",
  Processing_Error   = "#8E6C8A",
  Unclassified       = "#5D5A58"
)

PALETTE_PAM <- c(
  Canonical       = "#2E8B57",
  `Non-canonical` = "#F4A460",
  Disrupted       = "#A50026"
)

PALETTE_SEED <- c(
  Intact    = "#2E8B57",
  Disrupted = "#A50026"
)

PALETTE_PAIRS <- c(
  "Both homoeologs"       = "#2E8B57",
  "Only one gene in pair" = "#F4A460",
  "No specific gRNA"      = "#A50026"
)

PALETTE_SCORES <- c(
  score_cfd        = "#D62728",
  score_deephf     = "#2CA02C",
  score_deepspcas9 = "#1F77B4"
)

# ggpubr theme (matching original)
THEME_PUBLICATION <- theme_pubclean() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# ---- DATA LOADING FUNCTIONS -------------------------------------------------

#' Load Basic Classification Data
#'
#' @description Loads and validates gRNA classification data
#' @param classification_file Path to classification CSV file
#' @return data.table with classification data
#' @export
load_basic_data <- function(classification_file = "gRNA_classification_final.csv") {
  message("Loading classification file: ", classification_file)
  
  if (!file.exists(classification_file)) {
    stop("Classification file not found: ", classification_file)
  }
  
  classification_dt <- fread(classification_file, showProgress = FALSE)
  
  # Validate required columns
  if (!"gRNA_id" %in% names(classification_dt)) {
    stop("Column 'gRNA_id' not found in classification file")
  }
  if (!"category" %in% names(classification_dt)) {
    stop("Column 'category' not found in classification file")
  }
  
  # Remove duplicates
  classification_dt <- classification_dt[!duplicated(gRNA_id)]
  
  # Clean category names
  classification_dt[, category_clean := ifelse(
    is.na(category) | category == "",
    "Unclassified",
    category
  )]
  
  message("Loaded ", nrow(classification_dt), " unique gRNAs")
  
  return(classification_dt)
}

#' Load Full Data with Homoeologous Pairs
#'
#' @description Loads classification, gRNA database, and pairs data
#' @param classification_file Path to classification CSV
#' @param grna_db Path to gRNA database TSV
#' @param pairs_file Path to homoeologous pairs file
#' @param cfd_threshold CFD score threshold for filtering
#' @return List with multiple data.tables
#' @export
load_full_data <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05
) {
  
  # Load classification data
  classification_dt <- load_basic_data(classification_file)
  
  # Load gRNA database
  message("Loading gRNA database: ", grna_db)
  if (!file.exists(grna_db)) {
    stop("gRNA database file not found: ", grna_db)
  }
  grna_db_data <- fread(grna_db, showProgress = FALSE)
  
  # Validate database columns
  required_db_cols <- c("ID", "gene_id", "score_cfd", "score_cfd_offs", "cds_offs")
  missing_db_cols <- setdiff(required_db_cols, names(grna_db_data))
  if (length(missing_db_cols) > 0) {
    stop("Missing required columns in gRNA database: ", paste(missing_db_cols, collapse = ", "))
  }
  
  # Load pairs file
  message("Loading pairs file: ", pairs_file)
  if (!file.exists(pairs_file)) {
    stop("Pairs file not found: ", pairs_file)
  }
  pairs_raw <- fread(pairs_file, header = FALSE, fill = TRUE, showProgress = FALSE)
  
  if (ncol(pairs_raw) < 2) {
    stop("Pairs file must contain at least two columns (geneA, geneB)")
  }
  
  setnames(pairs_raw, 1:2, c("geneA", "geneB"))
  if (ncol(pairs_raw) >= 5) {
    setnames(pairs_raw, 5, "perc_identity")
  }
  
  keep_cols <- intersect(c("geneA", "geneB", "perc_identity"), names(pairs_raw))
  pairs_dt <- pairs_raw[, ..keep_cols]
  pairs_dt <- pairs_dt[!is.na(geneA) & !is.na(geneB)]
  pairs_dt[, pair_id := .I]
  
  if ("perc_identity" %in% names(pairs_dt)) {
    pairs_dt[, perc_identity := suppressWarnings(as.numeric(perc_identity))]
  }
  
  # Get genes from pairs
  pair_genes <- unique(c(pairs_dt$geneA, pairs_dt$geneB))
  pairs_grna_db <- grna_db_data[gene_id %in% pair_genes]
  pairs_grna_db <- pairs_grna_db[!duplicated(ID)]
  
  # Filter by off-target CFD threshold
  message("Filtering gRNAs by off-target CFD threshold: ", cfd_threshold)
  mask_dt <- pairs_grna_db[, .(ID, gene_id, cds_offs, score_cfd_offs, score_cfd)]
  long_dt <- cSplit(
    mask_dt,
    c("cds_offs", "score_cfd_offs"),
    sep = ";",
    direction = "long",
    type.convert = FALSE
  )
  long_dt[, score_cfd_offs := suppressWarnings(as.numeric(score_cfd_offs))]
  long_dt <- long_dt[!is.na(score_cfd_offs)]
  
  # Identify bad spacers (off-target hits)
  offtarget_dt <- long_dt[gene_id != cds_offs & score_cfd_offs >= cfd_threshold]
  bad_spacers <- unique(offtarget_dt$ID)
  eligible_grnas <- pairs_grna_db[!ID %in% bad_spacers]
  
  message("Filtered out ", length(bad_spacers), " gRNAs with off-target hits")
  message("Retained ", nrow(eligible_grnas), " eligible gRNAs")
  
  # Merge with classification
  filtered_classification_pairs <- merge(
    classification_dt,
    eligible_grnas,
    by.x = "gRNA_id",
    by.y = "ID",
    all.x = FALSE
  )
  
  # Classify pairs
  specific_genes <- filtered_classification_pairs[!is.na(gene_id), unique(gene_id)]
  pair_detail <- copy(pairs_dt)
  pair_detail[, `:=`(
    geneA_has_specific = geneA %in% specific_genes,
    geneB_has_specific = geneB %in% specific_genes
  )]
  pair_detail[, pair_category := ifelse(
    geneA_has_specific & geneB_has_specific, "Both homoeologs",
    ifelse(geneA_has_specific | geneB_has_specific, "Only one gene in pair", "No specific gRNA")
  )]
  
  # Calculate pair counts
  pair_counts <- pair_detail[, .(pair_count = .N), by = pair_category]
  total_pairs <- nrow(pair_detail)
  if (total_pairs > 0) {
    pair_counts[, percentage := round(pair_count / total_pairs * 100, 2)]
  }
  
  return(list(
    classification = classification_dt,
    pairs_grna_filtered = filtered_classification_pairs,
    pair_classification = pair_detail,
    pair_counts = pair_counts,
    pairs_dt = pairs_dt,
    eligible_grnas = eligible_grnas,
    bad_spacers = bad_spacers,
    pairs_grna_all = pairs_grna_db
  ))
}

# ---- HELPER FUNCTIONS -------------------------------------------------------

#' Check if Labels Need Rotation
#'
#' @param labels Character vector of labels
#' @param max_chars Maximum characters before rotation
#' @return Logical
needs_rotation <- function(labels, max_chars = 15) {
  any(nchar(as.character(labels)) > max_chars) || length(labels) > 5
}

#' Calculate Bar Width
#'
#' @param labels Character vector of labels
#' @param base_width Base width value
#' @param min_width Minimum width
#' @param max_width Maximum width
#' @return Numeric width value
calculate_bar_width <- function(labels, base_width = 0.4, min_width = 0.3, max_width = 0.8) {
  max_label_width <- max(nchar(as.character(labels)), na.rm = TRUE)
  max(min_width, min(max_width, base_width * (10 / max_label_width)))
}

#' Format Count and Percentage Label
#'
#' @param count Numeric count
#' @param percentage Numeric percentage
#' @param digits Number of decimal places
#' @param multiline Whether to use multiline format
#' @return Character label
format_count_pct_label <- function(count, percentage, digits = 2, multiline = TRUE) {
  pct_fmt <- ifelse(
    is.na(percentage),
    NA_character_,
    format(round(percentage, digits), nsmall = digits, trim = TRUE)
  )
  formatted_count <- scales::comma(count)
  
  if (multiline) {
    ifelse(
      is.na(pct_fmt),
      formatted_count,
      paste0(formatted_count, "\n(", pct_fmt, "%)")
    )
  } else {
    ifelse(
      is.na(pct_fmt),
      formatted_count,
      paste0(formatted_count, " (", pct_fmt, "%)")
    )
  }
}

#' Get Color Palette for Levels
#'
#' @param levels_vec Vector of levels
#' @param palette_lookup Named vector of colors
#' @param fallback_pal Fallback palette function
#' @return Character vector of colors
get_palette <- function(levels_vec, palette_lookup, fallback_pal = scales::hue_pal()) {
  if (is.factor(levels_vec)) {
    lv <- levels(levels_vec)
  } else {
    lv <- unique(as.character(levels_vec))
  }
  lv <- lv[!is.na(lv)]
  
  if (length(lv) == 0) {
    return(character())
  }
  
  colors <- palette_lookup[lv]
  missing <- is.na(colors)
  
  if (any(missing)) {
    colors[missing] <- fallback_pal(sum(missing))
  }
  
  return(unname(colors))
}

#' Make Threshold Scale
#'
#' @param values Numeric vector of threshold values
#' @return ggplot2 scale object
make_threshold_scale <- function(values) {
  if (is.null(values) || length(values) == 0) {
    return(scale_x_continuous())
  }
  
  values <- sort(unique(values))
  
  if (length(values) >= 2) {
    scale_x_continuous(
      breaks = values,
      limits = range(values),
      labels = scales::number_format(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.02))
    )
  } else {
    scale_x_continuous(
      breaks = values,
      labels = scales::number_format(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.02))
    )
  }
}

# ---- THRESHOLD ANALYSIS FUNCTIONS -------------------------------------------

#' Calculate Pair Category Counts by Threshold
#'
#' @description Calculates how many pairs fall into each category at different thresholds
#' @param grna_dt data.table of gRNA data
#' @param pairs_dt data.table of homoeologous pairs
#' @param score_columns Character vector of score column names
#' @param thresholds Numeric vector of threshold values
#' @param dataset_label Label for this dataset
#' @return data.table with pair counts by threshold
calculate_pair_category_counts_by_threshold <- function(
    grna_dt,
    pairs_dt,
    score_columns,
    thresholds,
    dataset_label
) {
  
  if (is.null(pairs_dt) || nrow(pairs_dt) == 0) {
    message("No homoeologous pairs available for threshold sweep analysis")
    return(data.table())
  }
  
  if (length(score_columns) == 0) {
    return(data.table())
  }
  
  pairs_dt_local <- as.data.table(copy(pairs_dt))
  if (!all(c("geneA", "geneB") %in% names(pairs_dt_local))) {
    stop("pairs_dt must contain 'geneA' and 'geneB' columns")
  }
  
  grna_dt_local <- as.data.table(copy(grna_dt))
  if (!"gene_id" %in% names(grna_dt_local)) {
    stop("gRNA data must contain 'gene_id' column")
  }
  
  grna_dt_local <- grna_dt_local[!duplicated(ID)]
  
  categories <- c("Both homoeologs", "Only one gene in pair", "No specific gRNA")
  total_pairs <- nrow(pairs_dt_local)
  thresholds <- sort(unique(thresholds))
  
  results <- vector("list", length(score_columns) * length(thresholds))
  idx <- 1L
  
  for (score_col in score_columns) {
    if (!score_col %in% names(grna_dt_local)) next
    
    score_values <- grna_dt_local[[score_col]]
    if (!is.numeric(score_values)) {
      score_values <- suppressWarnings(as.numeric(score_values))
      grna_dt_local[, (score_col) := score_values]
    }
    
    valid_idx <- !is.na(score_values)
    if (!any(valid_idx)) next
    
    grna_scores <- grna_dt_local[valid_idx, .(ID, gene_id, score_value = get(score_col))]
    
    for (th in thresholds) {
      genes_with_specific <- unique(grna_scores[score_value >= th, gene_id])
      
      if (length(genes_with_specific) == 0) {
        counts <- data.table(
          dataset = dataset_label,
          score_type = score_col,
          threshold = th,
          pair_category = categories,
          pair_count = c(0L, 0L, total_pairs)
        )
      } else {
        pair_cat_dt <- pairs_dt_local[, {
          geneA_has <- geneA %in% genes_with_specific
          geneB_has <- geneB %in% genes_with_specific
          category <- fcase(
            geneA_has & geneB_has, "Both homoeologs",
            geneA_has | geneB_has, "Only one gene in pair",
            default = "No specific gRNA"
          )
          .(pair_category = category)
        }]
        
        counts <- pair_cat_dt[, .(pair_count = .N), by = pair_category]
        counts <- merge(
          data.table(pair_category = categories),
          counts,
          by = "pair_category",
          all.x = TRUE
        )
        counts[is.na(pair_count), pair_count := 0L]
        counts[, `:=`(
          dataset = dataset_label,
          score_type = score_col,
          threshold = th
        )]
      }
      
      counts[, `:=`(
        pair_count = as.integer(pair_count),
        total_pairs = total_pairs,
        pair_percentage = if (total_pairs > 0) round(pair_count / total_pairs * 100, 2) else NA_real_
      )]
      
      setcolorder(counts, c(
        "dataset", "score_type", "threshold",
        "pair_category", "pair_count", "pair_percentage", "total_pairs"
      ))
      
      results[[idx]] <- counts
      idx <- idx + 1L
    }
  }
  
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(data.table())
  
  return(rbindlist(results, use.names = TRUE, fill = TRUE))
}

#' Calculate Retained gRNA Counts by Threshold
#'
#' @description Calculates how many gRNAs pass each threshold
#' @param grna_dt data.table of gRNA data
#' @param score_columns Character vector of score column names
#' @param thresholds Numeric vector of threshold values
#' @param dataset_label Label for this dataset
#' @return data.table with retained gRNA counts
calculate_retained_gRNA_counts <- function(
    grna_dt,
    score_columns,
    thresholds,
    dataset_label
) {
  
  if (length(score_columns) == 0) {
    return(data.table())
  }
  
  grna_dt_local <- as.data.table(copy(grna_dt))
  if (!"ID" %in% names(grna_dt_local)) {
    stop("gRNA data must contain 'ID' column")
  }
  
  grna_dt_local <- grna_dt_local[!duplicated(ID)]
  total_spacers <- uniqueN(grna_dt_local$ID)
  thresholds <- sort(unique(thresholds))
  
  results <- vector("list", length(score_columns) * length(thresholds))
  idx <- 1L
  
  for (score_col in score_columns) {
    if (!score_col %in% names(grna_dt_local)) next
    
    score_values <- grna_dt_local[[score_col]]
    if (!is.numeric(score_values)) {
      score_values <- suppressWarnings(as.numeric(score_values))
      grna_dt_local[, (score_col) := score_values]
    }
    
    valid_idx <- !is.na(score_values)
    if (!any(valid_idx)) next
    
    score_table <- grna_dt_local[valid_idx, .(ID, score_value = get(score_col))]
    
    for (th in thresholds) {
      retained <- score_table[score_value >= th, uniqueN(ID)]
      
      results[[idx]] <- data.table(
        dataset = dataset_label,
        score_type = score_col,
        threshold = th,
        retained_spacers = as.integer(retained),
        total_spacers = total_spacers,
        retained_percentage = if (total_spacers > 0) round(retained / total_spacers * 100, 2) else NA_real_
      )
      idx <- idx + 1L
    }
  }
  
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(data.table())
  
  return(rbindlist(results, use.names = TRUE, fill = TRUE))
}

#' Prepare Threshold Analysis
#'
#' @description Prepares data for threshold sweep analysis
#' @param pairs_dt data.table of homoeologous pairs
#' @param grna_all_dt data.table of all gRNAs
#' @param grna_filtered_dt data.table of filtered gRNAs
#' @param thresholds Numeric vector of thresholds to test
#' @param score_columns Character vector of score columns
#' @return List with threshold analysis data
#' @export
prepare_threshold_analysis <- function(
    pairs_dt,
    grna_all_dt,
    grna_filtered_dt,
    thresholds = seq(0.2, 1, by = 0.1),
    score_columns = c("score_cfd", "score_deephf", "score_deepspcas9")
) {
  
  if (is.null(pairs_dt) || nrow(pairs_dt) == 0) {
    message("No homoeologous pair data provided; skipping threshold sweep preparation")
    return(list(
      pair_counts = data.table(),
      spacer_counts = data.table(),
      dataset_summary = data.table(),
      thresholds = thresholds
    ))
  }
  
  pairs_dt_local <- as.data.table(copy(pairs_dt))
  
  dataset_definitions <- list(
    list(name = "Filtered", data = grna_filtered_dt),
    list(name = "All",      data = grna_all_dt)
  )
  
  pair_counts_list <- list()
  spacer_counts_list <- list()
  summary_list <- list()
  
  for (ds in dataset_definitions) {
    dataset_name <- ds$name
    grna_dt <- ds$data
    
    if (is.null(grna_dt) || nrow(grna_dt) == 0) {
      message("No gRNA data available for dataset '", dataset_name, "'. Skipping.")
      next
    }
    
    grna_dt <- as.data.table(copy(grna_dt))
    grna_dt <- grna_dt[!duplicated(ID)]
    
    if (!"gene_id" %in% names(grna_dt)) {
      message("gRNA data for dataset '", dataset_name, "' is missing 'gene_id'; skipping.")
      next
    }
    
    pair_genes <- unique(c(pairs_dt_local$geneA, pairs_dt_local$geneB))
    grna_dt <- grna_dt[gene_id %in% pair_genes]
    
    if (nrow(grna_dt) == 0) {
      message("No gRNAs linked to homoeologous genes for dataset '", dataset_name, "'.")
      next
    }
    
    score_candidates <- intersect(score_columns, names(grna_dt))
    if (length(score_candidates) == 0) {
      message("No requested score columns found for dataset '", dataset_name, "'.")
      next
    }
    
    available_scores <- character()
    for (score_col in score_candidates) {
      grna_dt[, (score_col) := suppressWarnings(as.numeric(get(score_col)))]
      if (any(!is.na(grna_dt[[score_col]]))) {
        available_scores <- c(available_scores, score_col)
      }
    }
    
    if (length(available_scores) == 0) {
      message("Score columns exist but contain only NA values for dataset '", dataset_name, "'.")
      next
    }
    
    pair_counts_dt <- calculate_pair_category_counts_by_threshold(
      grna_dt = grna_dt,
      pairs_dt = pairs_dt_local,
      score_columns = available_scores,
      thresholds = thresholds,
      dataset_label = dataset_name
    )
    
    spacer_counts_dt <- calculate_retained_gRNA_counts(
      grna_dt = grna_dt,
      score_columns = available_scores,
      thresholds = thresholds,
      dataset_label = dataset_name
    )
    
    if (nrow(pair_counts_dt) > 0) {
      pair_counts_list[[dataset_name]] <- pair_counts_dt
    }
    if (nrow(spacer_counts_dt) > 0) {
      spacer_counts_list[[dataset_name]] <- spacer_counts_dt
    }
    
    summary_list[[dataset_name]] <- data.table(
      dataset = dataset_name,
      total_pairs = nrow(pairs_dt_local),
      total_spacers = uniqueN(grna_dt$ID),
      available_scores = list(available_scores)
    )
  }
  
  return(list(
    pair_counts = if (length(pair_counts_list) > 0) rbindlist(pair_counts_list, use.names = TRUE, fill = TRUE) else data.table(),
    spacer_counts = if (length(spacer_counts_list) > 0) rbindlist(spacer_counts_list, use.names = TRUE, fill = TRUE) else data.table(),
    dataset_summary = if (length(summary_list) > 0) rbindlist(summary_list, use.names = TRUE, fill = TRUE) else data.table(),
    thresholds = thresholds
  ))
}

# ---- PLOT FUNCTIONS (ORIGINAL STYLING) -------------------------------------

#' Plot Global Category Distribution
#'
#' @description Creates bar plot of gRNA functional categories
#' @param classification_file Path to classification file
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_global_category_combo <- function(
    classification_file = "gRNA_classification_final.csv",
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_basic_data(classification_file)
  }
  
  # Calculate counts
  overall_counts <- data[, .(gRNA_count = .N), by = category_clean][order(-gRNA_count)]
  total_gRNAs <- sum(overall_counts$gRNA_count)
  overall_counts[, percentage := round(gRNA_count / total_gRNAs * 100, 2)]
  overall_counts[, category_clean := factor(category_clean, levels = category_clean)]
  overall_counts[, label_text := format_count_pct_label(gRNA_count, percentage, digits = 2)]
  
  # Check if rotation needed
  rotate_labels <- needs_rotation(levels(overall_counts$category_clean))
  bar_width <- calculate_bar_width(overall_counts$category_clean)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    overall_counts, 
    x = "category_clean", 
    y = "gRNA_count",
    fill = "category_clean",
    width = 0.9,
    alpha = 0.9,
    palette = get_palette(overall_counts$category_clean, PALETTE_CATEGORIES),
    xlab = FALSE,
    ylab = "Number of gRNAs",
    title = "Distribution of gRNA Functional Categories",
    subtitle = sprintf("Total unique gRNA: %s", scales::comma(total_gRNAs)),
    label = overall_counts$label_text,
    lab.pos = "out",
    lab.vjust = -0.6,
    lab.hjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.25))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "fig_all_grna_categories_combo.png"),
      p, 
      width = 6, 
      height = 5, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_all_grna_categories_combo.png"))
  }
  
  if (return_data) {
    return(list(plot = p, data = overall_counts))
  } else {
    return(p)
  }
}

#' Plot PAM Status Distribution
#'
#' @description Creates bar plot of PAM status
#' @param classification_file Path to classification file
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param prefix File prefix
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_pam_status <- function(
    classification_file = "gRNA_classification_final.csv",
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    prefix = "all_grna",
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_basic_data(classification_file)
  }
  
  if (!"pam_status" %in% names(data)) {
    message("No 'pam_status' column found in data")
    return(NULL)
  }
  
  # Filter and clean data
  pam_data <- data[!is.na(pam_status) & pam_status != ""]
  
  if (nrow(pam_data) == 0) {
    message("No PAM status data available")
    return(NULL)
  }
  
  pam_data[, pam_status_clean := trimws(pam_status)]
  
  # Standardize names
  replacements <- c(
    "non-canonical" = "Non-canonical",
    "Non_canonical" = "Non-canonical",
    "non_canonical" = "Non-canonical",
    "canonical" = "Canonical",
    "disrupted" = "Disrupted"
  )
  pam_data[, pam_status_clean := fcase(
    pam_status_clean %in% names(replacements), replacements[pam_status_clean],
    default = pam_status_clean
  )]
  
  # Calculate counts
  pam_counts <- pam_data[, .(count = .N), by = pam_status_clean][order(-count)]
  total_pam <- sum(pam_counts$count)
  pam_counts[, percentage := round(count / total_pam * 100, 2)]
  pam_counts[, label_text := format_count_pct_label(count, percentage, digits = 2)]
  
  # Order factor levels
  pam_order <- c("Canonical", "Non-canonical", "Disrupted")
  pam_counts[, pam_status_clean := factor(
    pam_status_clean, 
    levels = intersect(pam_order, pam_status_clean)
  )]
  
  bar_width <- calculate_bar_width(pam_counts$pam_status_clean)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    pam_counts, 
    x = "pam_status_clean", 
    y = "count",
    fill = "pam_status_clean",
    width = 0.9,
    alpha = 0.9,
    palette = get_palette(pam_counts$pam_status_clean, PALETTE_PAM),
    xlab = FALSE,
    ylab = "Count",
    title = "PAM Status Distribution",
    subtitle = sprintf("Total gRNAs with PAM status: %s", scales::comma(total_pam)),
    label = pam_counts$label_text,
    lab.pos = "out",
    lab.vjust = -0.5,
    lab.hjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.2))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, sprintf("fig_%s_pam_status.png", prefix)),
      p, 
      width = 5.5, 
      height = 5, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, sprintf("fig_%s_pam_status.png", prefix)))
  }
  
  if (return_data) {
    return(list(plot = p, data = pam_counts))
  } else {
    return(p)
  }
}

#' Plot Seed Status Distribution
#'
#' @description Creates bar plot of seed region status
#' @param classification_file Path to classification file
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param prefix File prefix
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_seed_status <- function(
    classification_file = "gRNA_classification_final.csv",
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    prefix = "all_grna",
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_basic_data(classification_file)
  }
  
  if (!"seed_status" %in% names(data)) {
    message("No 'seed_status' column found in data")
    return(NULL)
  }
  
  # Filter and clean data
  seed_data <- data[!is.na(seed_status) & seed_status != ""]
  
  if (nrow(seed_data) == 0) {
    message("No seed status data available")
    return(NULL)
  }
  
  seed_data[, seed_status_clean := trimws(seed_status)]
  
  # Standardize names
  replacements <- c("intact" = "Intact", "disrupted" = "Disrupted")
  seed_data[, seed_status_clean := fcase(
    seed_status_clean %in% names(replacements), replacements[seed_status_clean],
    default = seed_status_clean
  )]
  
  # Calculate counts
  seed_counts <- seed_data[, .(count = .N), by = seed_status_clean][order(-count)]
  total_seed <- sum(seed_counts$count)
  seed_counts[, percentage := round(count / total_seed * 100, 2)]
  seed_counts[, label_text := format_count_pct_label(count, percentage, digits = 2)]
  
  # Order factor levels
  seed_order <- c("Intact", "Disrupted")
  seed_counts[, seed_status_clean := factor(
    seed_status_clean, 
    levels = intersect(seed_order, seed_status_clean)
  )]
  
  bar_width <- calculate_bar_width(seed_counts$seed_status_clean)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    seed_counts, 
    x = "seed_status_clean", 
    y = "count",
    fill = "seed_status_clean",
    width = 0.9,
    alpha = 0.9,
    palette = get_palette(seed_counts$seed_status_clean, PALETTE_SEED),
    xlab = FALSE,
    ylab = "Count",
    title = "Seed Status Distribution",
    subtitle = sprintf("Total gRNAs with seed status: %s", scales::comma(total_seed)),
    label = seed_counts$label_text,
    lab.pos = "out",
    lab.vjust = -0.5,
    lab.hjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.2))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, sprintf("fig_%s_seed_status.png", prefix)),
      p, 
      width = 5, 
      height = 5, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, sprintf("fig_%s_seed_status.png", prefix)))
  }
  
  if (return_data) {
    return(list(plot = p, data = seed_counts))
  } else {
    return(p)
  }
}

#' Plot Top Details
#'
#' @description Creates horizontal bar plot of top detail types
#' @param classification_file Path to classification file
#' @param data Pre-loaded data (optional)
#' @param top_n Number of top details to show
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param prefix File prefix
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_top_details <- function(
    classification_file = "gRNA_classification_final.csv",
    data = NULL,
    top_n = 5,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    prefix = "all_grna",
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_basic_data(classification_file)
  }
  
  if (!"details" %in% names(data)) {
    message("No 'details' column found in data")
    return(NULL)
  }
  
  # Filter data
  details_data <- data[!is.na(details) & details != ""]
  
  if (nrow(details_data) == 0) {
    message("No details data available")
    return(NULL)
  }
  
  # Parse details
  all_details <- unlist(strsplit(details_data$details, ",", fixed = TRUE))
  all_details <- trimws(all_details)
  
  if (length(all_details) == 0) {
    message("No valid details found after parsing")
    return(NULL)
  }
  
  # Calculate counts
  details_counts <- data.table(detail = all_details)[, .(count = .N), by = detail][order(-count)]
  top_details <- head(details_counts, top_n)
  total_details <- sum(details_counts$count)
  top_details[, percentage := round(count / total_details * 100, 2)]
  top_details[, detail := factor(detail, levels = detail[order(count)])]
  top_details[, label_text := format_count_pct_label(count, percentage, digits = 2, multiline = FALSE)]
  
  bar_width <- calculate_bar_width(top_details$detail)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    top_details, 
    x = "detail", 
    y = "count",
    fill = "#8E6C8A",
    width = 0.9,
    alpha = 0.9,
    orientation = "horiz",
    xlab = "Detail Type",
    ylab = "Count",
    title = paste("Top", top_n, "Details"),
    subtitle = sprintf("Total details occurrences: %s", scales::comma(total_details)),
    label = top_details$label_text,
    lab.pos = "out",
    lab.hjust = -0.05,
    lab.vjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.5))
    ) +
    THEME_PUBLICATION +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, sprintf("fig_%s_top_details.png", prefix)),
      p, 
      width = 9, 
      height = 6, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, sprintf("fig_%s_top_details.png", prefix)))
  }
  
  if (return_data) {
    return(list(plot = p, data = top_details))
  } else {
    return(p)
  }
}

#' Plot Filtered Pairs Categories
#'
#' @description Creates bar plot of gRNA categories for filtered pairs
#' @param classification_file Path to classification file
#' @param grna_db Path to gRNA database
#' @param pairs_file Path to pairs file
#' @param cfd_threshold CFD threshold
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_filtered_pairs_categories_combo <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05,
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  }
  
  if (nrow(data$pairs_grna_filtered) == 0) {
    message("No data available for filtered pairs categories plot")
    return(NULL)
  }
  
  # Calculate counts
  category_counts <- data$pairs_grna_filtered[, .(gRNA_count = .N), by = category_clean][order(-gRNA_count)]
  total_gRNAs <- sum(category_counts$gRNA_count)
  category_counts[, percentage := round(gRNA_count / total_gRNAs * 100, 2)]
  category_counts[, category_clean := factor(category_clean, levels = category_clean)]
  category_counts[, label_text := format_count_pct_label(gRNA_count, percentage, digits = 2)]
  
  rotate_labels <- needs_rotation(levels(category_counts$category_clean))
  bar_width <- calculate_bar_width(category_counts$category_clean)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    category_counts, 
    x = "category_clean", 
    y = "gRNA_count",
    fill = "category_clean",
    width = 0.9,
    alpha = 0.9,
    palette = get_palette(category_counts$category_clean, PALETTE_CATEGORIES),
    xlab = FALSE,
    ylab = "Number of gRNAs",
    title = sprintf("gRNA Categories after CFD Filtering\n(Homoeologous Genes)"),
    subtitle = sprintf(
      "Unique gRNAs retained: %s | off-target CFD threshold: %.2f",
      scales::comma(total_gRNAs), cfd_threshold
    ),
    label = category_counts$label_text,
    lab.pos = "out",
    lab.vjust = -0.6,
    lab.hjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.25))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "fig_filtered_pairs_categories_combo.png"),
      p, 
      width = 6.5, 
      height = 5, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_filtered_pairs_categories_combo.png"))
  }
  
  if (return_data) {
    return(list(plot = p, data = category_counts))
  } else {
    return(p)
  }
}

#' Plot Pair Classification (Fully Functional Only)
#'
#' @description Creates bar plot of homoeologous pair classification
#' @param classification_file Path to classification file
#' @param grna_db Path to gRNA database
#' @param pairs_file Path to pairs file
#' @param cfd_threshold CFD threshold
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_pair_classification_combo <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05,
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  }
  
  if (is.null(data$pair_classification) || nrow(data$pair_classification) == 0) {
    message("No pair classification data available")
    return(NULL)
  }
  
  # Filter to only Fully_Functional gRNAs
  functional_genes <- data$pairs_grna_filtered[
    category_clean == "Fully_Functional", 
    unique(gene_id)
  ]
  
  pair_detail_filtered <- copy(data$pair_classification)
  pair_detail_filtered[, `:=`(
    geneA_has_functional = geneA %in% functional_genes,
    geneB_has_functional = geneB %in% functional_genes
  )]
  pair_detail_filtered[, pair_category := ifelse(
    geneA_has_functional & geneB_has_functional, "Both homoeologs",
    ifelse(geneA_has_functional | geneB_has_functional, "Only one gene in pair", "No specific gRNA")
  )]
  
  # Calculate counts
  pair_counts_filtered <- pair_detail_filtered[, .(pair_count = .N), by = pair_category]
  total_pairs <- nrow(pair_detail_filtered)
  if (total_pairs > 0) {
    pair_counts_filtered[, percentage := round(pair_count / total_pairs * 100, 2)]
  }
  
  pair_counts <- pair_counts_filtered
  pair_counts[, pair_category := factor(
    pair_category,
    levels = c("Both homoeologs", "Only one gene in pair", "No specific gRNA")
  )]
  pair_counts[, label_text := format_count_pct_label(pair_count, percentage, digits = 2)]
  
  bar_width <- calculate_bar_width(pair_counts$pair_category)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    pair_counts, 
    x = "pair_category", 
    y = "pair_count",
    fill = "pair_category",
    width = 0.9,
    alpha = 0.9,
    palette = get_palette(pair_counts$pair_category, PALETTE_PAIRS),
    xlab = FALSE,
    ylab = "Number of gene pairs",
    title = "Homoeologous Pairs Classification\n(Only Fully Functional gRNAs)",
    subtitle = sprintf(
      "Total pairs analyzed: %s | off-target CFD threshold: %.2f",
      scales::comma(total_pairs), cfd_threshold
    ),
    label = pair_counts$label_text,
    lab.pos = "out",
    lab.vjust = -0.6,
    lab.hjust = 0.5
  ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.2))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    font("xy.text", size = 12)
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "fig_pair_classification_combo.png"),
      p, 
      width = 6, 
      height = 5, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_pair_classification_combo.png"))
  }
  
  if (return_data) {
    return(list(plot = p, data = pair_counts))
  } else {
    return(p)
  }
}

#' Plot Scores and Identity Boxplots
#'
#' @description Creates boxplots for gRNA scores and sequence identity
#' @param classification_file Path to classification file
#' @param grna_db Path to gRNA database
#' @param pairs_file Path to pairs file
#' @param cfd_threshold CFD threshold
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_scores_and_identity_boxplots <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05,
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  }
  
  if (nrow(data$pairs_grna_filtered) == 0) {
    message("No data available for scores and identity boxplots")
    return(NULL)
  }
  
  plot_data <- data$pairs_grna_filtered
  
  # Identify available score columns
  possible_score_columns <- c(
    "score_deepspcas9", "score_deephf", "score_cfd",
    "deepspcas9_score", "deephf_score", "cfd_score"
  )
  score_columns <- intersect(possible_score_columns, names(plot_data))
  
  if (length(score_columns) == 0) {
    message("No score columns available for boxplots")
    message("Available columns: ", paste(names(plot_data), collapse = ", "))
    return(NULL)
  }
  
  # Reshape data for plotting
  score_long <- melt(
    plot_data,
    measure.vars = score_columns,
    variable.name = "score_type",
    value.name = "score_value"
  )[!is.na(score_value)]
  
  # Create readable labels
  score_labels <- c(
    score_deepspcas9 = "DeepSpCas9",
    score_deephf     = "DeepHF",
    score_cfd        = "Aggregate CFD",
    deepspcas9_score = "DeepSpCas9",
    deephf_score     = "DeepHF",
    cfd_score        = "Aggregate CFD"
  )
  
  score_long[, score_type := factor(
    score_labels[score_type], 
    levels = unique(score_labels[score_columns])
  )]
  
  # Create score boxplot (ORIGINAL STYLE)
  plot_scores <- ggboxplot(
    score_long, 
    x = "score_type", 
    y = "score_value",
    fill = "grey80",
    color = "black",
    xlab = FALSE,
    ylab = "Score",
    bxp.errorbar = TRUE,
    outlier.size = 0.1
  ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10), 
      limits = c(0, 1), 
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    THEME_PUBLICATION +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    font("xy.text", size = 12)
  
  # Check for identity data
  identity_available <- "perc_identity" %in% names(data$pairs_dt) &&
    any(!is.na(data$pairs_dt$perc_identity))
  
  if (identity_available) {
    identity_dt <- data$pair_classification[
      !is.na(perc_identity),
      .(pair_category, perc_identity)
    ]
    
    identity_dt[, pair_category := factor(
      pair_category,
      levels = c("No specific gRNA", "Only one gene in pair", "Both homoeologs")
    )]
    
    # Calculate sample sizes
    identity_counts <- identity_dt[, .(
      N = .N, 
      y = max(perc_identity, na.rm = TRUE) + 3
    ), by = pair_category]
    
    # Create identity boxplot (ORIGINAL STYLE)
    plot_identity <- ggboxplot(
      identity_dt, 
      x = "pair_category", 
      y = "perc_identity",
      fill = "grey80",
      color = "black",
      xlab = FALSE,
      ylab = "CDS nucleotide sequence identity, %",
      bxp.errorbar = TRUE,
      outlier.size = 0.1
    ) +
      geom_text(
        data = identity_counts,
        aes(x = pair_category, y = y, label = scales::comma(N)),
        vjust = 0,
        size = 3.5
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 10), 
        expand = expansion(mult = c(0.02, 0.08))
      ) +
      THEME_PUBLICATION +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank()
      ) +
      font("xy.text", size = 12)
    
    # Combine plots
    combined_plot <- ggarrange(
      plot_scores, plot_identity,
      ncol = 2, nrow = 1,
      common.legend = FALSE
    )
    
    print(combined_plot)
    
    if (save_plot) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      ggsave(
        file.path(output_dir, "fig_scores_and_identity_boxplots.png"),
        combined_plot, 
        width = 13, 
        height = 6, 
        dpi = plot_dpi, 
        bg = "white"
      )
      message("Plot saved to: ", file.path(output_dir, "fig_scores_and_identity_boxplots.png"))
    }
    
    if (return_data) {
      return(list(
        plot = combined_plot, 
        scores_data = score_long, 
        identity_data = identity_dt
      ))
    } else {
      return(combined_plot)
    }
    
  } else {
    # Return scores plot only
    print(plot_scores)
    
    if (save_plot) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      ggsave(
        file.path(output_dir, "fig_scores_boxplots.png"),
        plot_scores, 
        width = 6, 
        height = 5, 
        dpi = plot_dpi, 
        bg = "white"
      )
      message("Plot saved to: ", file.path(output_dir, "fig_scores_boxplots.png"))
    }
    
    if (return_data) {
      return(list(plot = plot_scores, scores_data = score_long, identity_data = NULL))
    } else {
      return(plot_scores)
    }
  }
}

#' Plot Homoeolog Pair All Categories (Swapped Layout)
#'
#' @description Creates grouped bar plot showing all categories for homoeologous pairs
#' @param classification_file Path to classification file
#' @param grna_db Path to gRNA database
#' @param pairs_file Path to pairs file
#' @param cfd_threshold CFD threshold
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_homoeolog_pair_all_categories_swapped <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05,
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    return_data = FALSE
) {
  
  if (is.null(data)) {
    data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  }
  
  if (nrow(data$pairs_grna_filtered) == 0) {
    message("No data available for homoeolog pair categories plot")
    return(NULL)
  }
  
  # Identify genes with each category
  functional_genes <- data$pairs_grna_filtered[
    , .(has_fully_functional = any(category_clean == "Fully_Functional"),
        has_reduced_efficiency = any(category_clean == "Reduced_Efficiency"),
        has_critical_failure = any(category_clean == "Critical_Failure")),
    by = gene_id
  ]
  
  pairs_dt_extended <- copy(data$pairs_dt)
  if (!"pair_id" %in% names(pairs_dt_extended)) {
    pairs_dt_extended[, pair_id := .I]
  }
  
  # Merge with functional genes
  pairs_with_scores <- merge(
    pairs_dt_extended, functional_genes,
    by.x = "geneA", by.y = "gene_id", all.x = TRUE
  )
  pairs_with_scores <- merge(
    pairs_with_scores, functional_genes,
    by.x = "geneB", by.y = "gene_id", all.x = TRUE,
    suffixes = c("_A", "_B")
  )
  
  # Fix NA values
  cols_to_fix <- c(
    "has_fully_functional_A", "has_reduced_efficiency_A", "has_critical_failure_A",
    "has_fully_functional_B", "has_reduced_efficiency_B", "has_critical_failure_B"
  )
  for (col in cols_to_fix) {
    pairs_with_scores[is.na(get(col)), (col) := FALSE]
  }
  
  # Classify pairs for each category
  pairs_with_scores[, category_fully_functional := fcase(
    has_fully_functional_A & has_fully_functional_B, "Both homoeologs",
    has_fully_functional_A | has_fully_functional_B, "Only one gene in pair",
    default = "No specific gRNA"
  )]
  pairs_with_scores[, category_reduced_efficiency := fcase(
    has_reduced_efficiency_A & has_reduced_efficiency_B, "Both homoeologs",
    has_reduced_efficiency_A | has_reduced_efficiency_B, "Only one gene in pair",
    default = "No specific gRNA"
  )]
  pairs_with_scores[, category_critical_failure := fcase(
    has_critical_failure_A & has_critical_failure_B, "Both homoeologs",
    has_critical_failure_A | has_critical_failure_B, "Only one gene in pair",
    default = "No specific gRNA"
  )]
  
  # Combine results
  fully_functional <- pairs_with_scores[, .(pair_id, category = category_fully_functional, grna_category = "Fully_Functional")]
  reduced_efficiency <- pairs_with_scores[, .(pair_id, category = category_reduced_efficiency, grna_category = "Reduced_Efficiency")]
  critical_failure <- pairs_with_scores[, .(pair_id, category = category_critical_failure, grna_category = "Critical_Failure")]
  
  all_categories_results <- rbindlist(list(fully_functional, reduced_efficiency, critical_failure))
  
  # Calculate statistics
  pair_stats_all_categories <- all_categories_results[, .(
    count = .N,
    percentage = round(.N / nrow(data$pairs_dt) * 100, 2)
  ), by = .(grna_category, category)]
  
  pair_stats_all_categories[, grna_category := factor(
    grna_category,
    levels = c("Fully_Functional", "Reduced_Efficiency", "Critical_Failure")
  )]
  pair_stats_all_categories[, category := factor(
    category,
    levels = c("Both homoeologs", "Only one gene in pair", "No specific gRNA")
  )]
  pair_stats_all_categories[, label_text := format_count_pct_label(count, percentage, digits = 1)]
  
  bar_width <- calculate_bar_width(pair_stats_all_categories$category)
  
  # Create plot (ORIGINAL STYLE)
  p <- ggbarplot(
    pair_stats_all_categories, 
    x = "category", 
    y = "count",
    fill = "grna_category",
    width = 0.9,
    alpha = 0.8,
    palette = get_palette(pair_stats_all_categories$grna_category, PALETTE_CATEGORIES),
    position = position_dodge(0.9),
    xlab = "Pair Category",
    ylab = "Number of Pairs",
    title = "Distribution of Homoeologous Pairs by gRNA Availability",
    subtitle = sprintf("Off-target CFD threshold: %.2f", cfd_threshold),
    legend = "bottom",
    legend.title = "gRNA Category"
  ) +
    geom_text(
      aes(label = label_text, group = grna_category),
      position = position_dodge(width = 0.9),
      vjust = -0.3,
      size = 3.2,
      fontface = "bold",
      lineheight = 0.9
    ) +
    scale_y_continuous(
      labels = scales::comma, 
      expand = expansion(mult = c(0, 0.2))
    ) +
    THEME_PUBLICATION +
    font("xy.text", size = 12) +
    theme(panel.grid.major.x = element_blank())
  
  print(p)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "fig_all_categories_analysis_swapped.png"),
      p, 
      width = 8, 
      height = 6, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_all_categories_analysis_swapped.png"))
  }
  
  if (return_data) {
    return(list(plot = p, data = pair_stats_all_categories))
  } else {
    return(p)
  }
}

#' Plot Pair Threshold Sweep
#'
#' @description Shows how pair classification changes with score thresholds
#' @param threshold_analysis Output from prepare_threshold_analysis()
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param dataset_titles Named vector of dataset titles
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_pair_threshold_sweep <- function(
    threshold_analysis,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    dataset_titles = c(
      Filtered = sprintf("Classification of Gene Pairs\n(Off-target CFD filtered)"),
      All = sprintf("Classification of Gene Pairs\n(Not filtered by off-target CFD)")
    ),
    return_data = FALSE
) {
  
  pair_counts <- threshold_analysis$pair_counts
  dataset_summary <- threshold_analysis$dataset_summary
  
  if (is.null(pair_counts) || nrow(pair_counts) == 0 || 
      is.null(dataset_summary) || nrow(dataset_summary) == 0) {
    message("No pair threshold sweep data available for plotting")
    return(NULL)
  }
  
  # Order factor levels
  pair_counts[, pair_category := factor(
    pair_category,
    levels = c("Both homoeologs", "Only one gene in pair", "No specific gRNA")
  )]
  
  dataset_order <- c("Filtered", "All")
  dataset_levels <- dataset_order[dataset_order %in% dataset_summary$dataset]
  dataset_plots <- list()
  
  for (ds in dataset_levels) {
    ds_counts <- copy(pair_counts[dataset == ds])
    if (nrow(ds_counts) == 0) next
    
    score_levels <- dataset_summary[dataset == ds, unlist(available_scores, use.names = FALSE)]
    ds_counts[, score_type := factor(score_type, levels = score_levels)]
    
    total_spacers <- dataset_summary[dataset == ds, unique(total_spacers)]
    threshold_values <- threshold_analysis$thresholds
    
    if (is.null(threshold_values) || length(threshold_values) == 0) {
      threshold_values <- sort(unique(ds_counts$threshold))
    }
    
    # Create plot for this dataset (ORIGINAL STYLE)
    ds_plot <- ggplot(ds_counts, aes(x = threshold, y = pair_count, color = pair_category)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.2) +
      facet_wrap(~score_type, nrow = 1, scales = "fixed") +
      scale_color_manual(values = PALETTE_PAIRS) +
      make_threshold_scale(threshold_values) +
      scale_y_continuous(
        labels = scales::comma, 
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = dataset_titles[ds],
        subtitle = sprintf("gRNA count: %s", scales::comma(total_spacers)),
        x = "Score Threshold",
        y = "Number of Pairs",
        color = "Category"
      ) +
      THEME_PUBLICATION +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 13),
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        panel.spacing = unit(0.5, "cm"),
        axis.title.y = element_text(margin = margin(r = 10))
      ) +
      guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
      font("xy.text", size = 12)
    
    dataset_plots[[ds]] <- ds_plot
  }
  
  dataset_plots <- Filter(Negate(is.null), dataset_plots)
  
  if (length(dataset_plots) == 0) {
    message("No dataset-specific plots could be generated for pair threshold sweep")
    return(NULL)
  }
  
  # Combine plots
  combined_plot <- ggarrange(
    plotlist = dataset_plots,
    ncol = length(dataset_plots),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  print(combined_plot)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = file.path(output_dir, "fig_pair_threshold_classification.png"),
      plot = combined_plot,
      width = 14,
      height = 6.5,
      dpi = plot_dpi,
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_pair_threshold_classification.png"))
  }
  
  if (return_data) {
    return(list(
      plot = combined_plot,
      data = pair_counts,
      dataset_summary = dataset_summary
    ))
  } else {
    return(combined_plot)
  }
}

#' Plot Spacer Retention Threshold Sweep
#'
#' @description Shows how many gRNAs are retained at different thresholds
#' @param threshold_analysis Output from prepare_threshold_analysis()
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @param dataset_titles Named vector of dataset titles
#' @param return_data Whether to return data
#' @return ggplot object or list
#' @export
plot_spacer_threshold_sweep <- function(
    threshold_analysis,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300,
    dataset_titles = c(
      Filtered = sprintf("Retained gRNAs in pairs at different Score Thresholds\n(Off-target CFD filtered)"),
      All = sprintf("Retained gRNAs in pairs at different Score Thresholds\n(Not filtered by off-target CFD)")
    ),
    return_data = FALSE
) {
  
  spacer_counts <- threshold_analysis$spacer_counts
  dataset_summary <- threshold_analysis$dataset_summary
  
  if (is.null(spacer_counts) || nrow(spacer_counts) == 0 || 
      is.null(dataset_summary) || nrow(dataset_summary) == 0) {
    message("No spacer retention data available for plotting")
    return(NULL)
  }
  
  dataset_order <- c("Filtered", "All")
  dataset_levels <- dataset_order[dataset_order %in% dataset_summary$dataset]
  dataset_plots <- list()
  
  for (ds in dataset_levels) {
    ds_counts <- copy(spacer_counts[dataset == ds])
    if (nrow(ds_counts) == 0) next
    
    score_levels <- dataset_summary[dataset == ds, unlist(available_scores, use.names = FALSE)]
    ds_counts[, score_type := factor(score_type, levels = score_levels)]
    
    total_spacers <- dataset_summary[dataset == ds, unique(total_spacers)]
    threshold_values <- threshold_analysis$thresholds
    
    if (is.null(threshold_values) || length(threshold_values) == 0) {
      threshold_values <- sort(unique(ds_counts$threshold))
    }
    
    # Create plot for this dataset (ORIGINAL STYLE)
    ds_plot <- ggplot(ds_counts, aes(x = threshold, y = retained_spacers, color = score_type)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.2) +
      scale_color_manual(values = PALETTE_SCORES) +
      make_threshold_scale(threshold_values) +
      scale_y_continuous(
        labels = scales::comma, 
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = dataset_titles[ds],
        subtitle = sprintf("gRNA count: %s", scales::comma(total_spacers)),
        x = "Score Threshold",
        y = "Retained gRNA count",
        color = "Score Type"
      ) +
      THEME_PUBLICATION +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10))
      ) +
      guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
      font("xy.text", size = 12)
    
    dataset_plots[[ds]] <- ds_plot
  }
  
  dataset_plots <- Filter(Negate(is.null), dataset_plots)
  
  if (length(dataset_plots) == 0) {
    message("No dataset-specific plots could be generated for spacer retention sweep")
    return(NULL)
  }
  
  # Combine plots
  combined_plot <- ggarrange(
    plotlist = dataset_plots,
    ncol = length(dataset_plots),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  print(combined_plot)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = file.path(output_dir, "fig_spacer_threshold_retention.png"),
      plot = combined_plot,
      width = 14,
      height = 6.5,
      dpi = plot_dpi,
      bg = "white"
    )
    message("Plot saved to: ", file.path(output_dir, "fig_spacer_threshold_retention.png"))
  }
  
  if (return_data) {
    return(list(
      plot = combined_plot,
      data = spacer_counts,
      dataset_summary = dataset_summary
    ))
  } else {
    return(combined_plot)
  }
}

#' Create Combined 2x2 Figure for All gRNAs
#'
#' @description Creates a 2x2 panel figure for publication
#' @param classification_file Path to classification file
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @return Combined ggplot object
#' @export
create_combined_all_grna_plots <- function(
    classification_file = "gRNA_classification_final.csv",
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300
) {
  
  if (is.null(data)) {
    data <- load_basic_data(classification_file)
  }
  
  message("Creating combined 2x2 plot for all gRNAs...")
  
  # Generate individual plots
  p1 <- plot_global_category_combo(data = data, save_plot = FALSE, return_data = FALSE)
  p2 <- plot_pam_status(data = data, save_plot = FALSE, prefix = "all_grna", return_data = FALSE)
  p3 <- plot_seed_status(data = data, save_plot = FALSE, prefix = "all_grna", return_data = FALSE)
  p4 <- plot_top_details(data = data, save_plot = FALSE, prefix = "all_grna", return_data = FALSE)
  
  # Combine into 2x2 grid
  combined_plot <- ggarrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    common.legend = FALSE
  )
  
  # Add overall title
  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob(
      "Comprehensive gRNA Analysis - All gRNAs", 
      face = "bold", 
      size = 16
    )
  )
  
  print(combined_plot)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "combined_all_grna_2x2.png"),
      combined_plot, 
      width = 16, 
      height = 12, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Combined plot saved to: ", file.path(output_dir, "combined_all_grna_2x2.png"))
  }
  
  return(combined_plot)
}

#' Create Combined 2x2 Figure for Filtered Pairs
#'
#' @description Creates a 2x2 panel figure for filtered pairs
#' @param classification_file Path to classification file
#' @param grna_db Path to gRNA database
#' @param pairs_file Path to pairs file
#' @param cfd_threshold CFD threshold
#' @param data Pre-loaded data (optional)
#' @param save_plot Whether to save plot
#' @param output_dir Output directory
#' @param plot_dpi Resolution
#' @return Combined ggplot object
#' @export
create_combined_filtered_pairs_plots <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    cfd_threshold = 0.05,
    data = NULL,
    save_plot = FALSE,
    output_dir = "plots",
    plot_dpi = 300
) {
  
  if (is.null(data)) {
    data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  }
  
  message("Creating combined 2x2 plot for filtered pairs...")
  
  # Generate individual plots
  p1 <- plot_filtered_pairs_categories_combo(data = data, save_plot = FALSE, return_data = FALSE)
  p2 <- plot_pam_status(data = data$pairs_grna_filtered, save_plot = FALSE, 
                        prefix = "filtered_pairs", return_data = FALSE)
  p3 <- plot_seed_status(data = data$pairs_grna_filtered, save_plot = FALSE, 
                         prefix = "filtered_pairs", return_data = FALSE)
  p4 <- plot_top_details(data = data$pairs_grna_filtered, save_plot = FALSE, 
                         prefix = "filtered_pairs", return_data = FALSE)
  
  # Combine into 2x2 grid
  combined_plot <- ggarrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    common.legend = FALSE
  )
  
  # Add overall title
  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob(
      "Comprehensive gRNA Analysis - Filtered Pairs", 
      face = "bold", 
      size = 16
    )
  )
  
  print(combined_plot)
  
  if (save_plot) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      file.path(output_dir, "combined_filtered_pairs_2x2.png"),
      combined_plot, 
      width = 16, 
      height = 12, 
      dpi = plot_dpi, 
      bg = "white"
    )
    message("Combined plot saved to: ", file.path(output_dir, "combined_filtered_pairs_2x2.png"))
  }
  
  return(combined_plot)
}

# ---- COMPLETE ANALYSIS PIPELINE ---------------------------------------------

#' Run Complete Analysis Pipeline
#'
#' @description Executes full analysis workflow with all plots and statistics
#' @param classification_file Path to classification CSV
#' @param grna_db Path to gRNA database TSV
#' @param pairs_file Path to homoeologous pairs file
#' @param output_dir Output directory for all results
#' @param cfd_threshold CFD score threshold for filtering
#' @param plot_dpi Resolution for saved plots
#' @param save_plots Whether to save plots
#' @param save_csv Whether to save CSV tables
#' @return List containing all analysis results
#' @export
run_complete_analysis <- function(
    classification_file = "gRNA_classification_final.csv",
    grna_db = "Cbp_msk_all_genes_grna.filtered.article.tsv",
    pairs_file = "pairs_perc_id_genes.txt",
    output_dir = "gRNA_functionality_summary",
    cfd_threshold = 0.05,
    plot_dpi = 300,
    save_plots = TRUE,
    save_csv = TRUE
) {
  
  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  plot_dir <- file.path(output_dir, "plots")
  table_dir <- file.path(output_dir, "tables")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Output directory: ", normalizePath(output_dir))
  message("Plots directory: ", normalizePath(plot_dir))
  message("Tables directory: ", normalizePath(table_dir))
  
  message("\n=== LOADING DATA ===")
  basic_data <- load_basic_data(classification_file)
  full_data <- load_full_data(classification_file, grna_db, pairs_file, cfd_threshold)
  
  # Prepare threshold analysis
  threshold_analysis <- prepare_threshold_analysis(
    pairs_dt = full_data$pairs_dt,
    grna_all_dt = full_data$pairs_grna_all,
    grna_filtered_dt = full_data$eligible_grnas,
    thresholds = seq(0.2, 1.0, by = 0.1),
    score_columns = c("score_cfd", "score_deephf", "score_deepspcas9")
  )
  
  if (save_csv) {
    dedup_file <- file.path(table_dir, "gRNA_classification_deduplicated.csv")
    fwrite(basic_data, dedup_file)
    message("Deduplicated classification saved to: ", dedup_file)
  }
  
  results <- list()
  
  message("\n=== GENERATING ALL PLOTS AND STATISTICS ===")
  
  message("1. Generating global category analysis...")
  results$global_categories <- plot_global_category_combo(
    data = basic_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi, return_data = TRUE
  )
  
  message("2. Generating PAM status analysis for all gRNAs...")
  results$pam_status_all <- plot_pam_status(
    data = basic_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "all_grna", return_data = TRUE
  )
  
  message("3. Generating seed status analysis for all gRNAs...")
  results$seed_status_all <- plot_seed_status(
    data = basic_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "all_grna", return_data = TRUE
  )
  
  message("4. Generating top details analysis for all gRNAs...")
  results$top_details_all <- plot_top_details(
    data = basic_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "all_grna", return_data = TRUE
  )
  
  message("5. Generating filtered pairs categories analysis...")
  results$filtered_categories <- plot_filtered_pairs_categories_combo(
    data = full_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi, return_data = TRUE
  )
  
  message("6. Generating pair classification analysis (Fully_Functional only)...")
  results$pair_classification <- plot_pair_classification_combo(
    data = full_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi, return_data = TRUE
  )
  
  message("7. Generating PAM status analysis for filtered pairs...")
  results$pam_status_filtered <- plot_pam_status(
    data = full_data$pairs_grna_filtered, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "filtered_pairs", return_data = TRUE
  )
  
  message("8. Generating seed status analysis for filtered pairs...")
  results$seed_status_filtered <- plot_seed_status(
    data = full_data$pairs_grna_filtered, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "filtered_pairs", return_data = TRUE
  )
  
  message("9. Generating top details analysis for filtered pairs...")
  results$top_details_filtered <- plot_top_details(
    data = full_data$pairs_grna_filtered, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi,
    prefix = "filtered_pairs", return_data = TRUE
  )
  
  message("10. Generating scores and identity boxplots...")
  results$scores_identity <- plot_scores_and_identity_boxplots(
    data = full_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi, return_data = TRUE
  )
  
  message("11. Generating homoeolog pair categories swapped analysis...")
  results$homoeolog_categories_swapped <- plot_homoeolog_pair_all_categories_swapped(
    data = full_data, save_plot = save_plots,
    output_dir = plot_dir, plot_dpi = plot_dpi, return_data = TRUE
  )
  
  message("12. Generating gene pair classification threshold sweep...")
  if (!is.null(threshold_analysis$pair_counts) && nrow(threshold_analysis$pair_counts) > 0) {
    results$pair_threshold_sweep <- plot_pair_threshold_sweep(
      threshold_analysis = threshold_analysis,
      save_plot = save_plots,
      output_dir = plot_dir,
      plot_dpi = plot_dpi,
      return_data = TRUE
    )
  } else {
    message("No pair threshold sweep data available; skipping plot.")
    results$pair_threshold_sweep <- NULL
  }
  
  message("13. Generating spacer retention threshold sweep...")
  if (!is.null(threshold_analysis$spacer_counts) && nrow(threshold_analysis$spacer_counts) > 0) {
    results$spacer_threshold_sweep <- plot_spacer_threshold_sweep(
      threshold_analysis = threshold_analysis,
      save_plot = save_plots,
      output_dir = plot_dir,
      plot_dpi = plot_dpi,
      return_data = TRUE
    )
  } else {
    message("No spacer retention data available; skipping plot.")
    results$spacer_threshold_sweep <- NULL
  }
  
  message("14. Generating combined 2x2 plot for all gRNAs...")
  results$combined_all_grna <- create_combined_all_grna_plots(
    data = basic_data, 
    save_plot = save_plots,
    output_dir = plot_dir, 
    plot_dpi = plot_dpi
  )
  
  message("15. Generating combined 2x2 plot for filtered pairs...")
  results$combined_filtered_pairs <- create_combined_filtered_pairs_plots(
    data = full_data, 
    save_plot = save_plots,
    output_dir = plot_dir, 
    plot_dpi = plot_dpi
  )
  
  results$threshold_analysis <- threshold_analysis
  
  # Save CSV files
  if (save_csv) {
    message("\n=== SAVING CSV FILES ===")
    
    if (!is.null(results$global_categories$data)) {
      fwrite(results$global_categories$data, file.path(table_dir, "stat_all_grna_categories.csv"))
    }
    if (!is.null(results$pam_status_all$data)) {
      fwrite(results$pam_status_all$data, file.path(table_dir, "stat_all_grna_pam_status.csv"))
    }
    if (!is.null(results$seed_status_all$data)) {
      fwrite(results$seed_status_all$data, file.path(table_dir, "stat_all_grna_seed_status.csv"))
    }
    if (!is.null(results$top_details_all$data)) {
      fwrite(results$top_details_all$data, file.path(table_dir, "stat_all_grna_top_details.csv"))
    }
    if (!is.null(results$filtered_categories$data)) {
      fwrite(results$filtered_categories$data, file.path(table_dir, "stat_filtered_pairs_categories.csv"))
    }
    if (!is.null(results$pair_classification$data)) {
      fwrite(results$pair_classification$data, file.path(table_dir, "stat_pair_classification.csv"))
    }
    if (!is.null(results$pam_status_filtered$data)) {
      fwrite(results$pam_status_filtered$data, file.path(table_dir, "stat_filtered_pairs_pam_status.csv"))
    }
    if (!is.null(results$seed_status_filtered$data)) {
      fwrite(results$seed_status_filtered$data, file.path(table_dir, "stat_filtered_pairs_seed_status.csv"))
    }
    if (!is.null(results$top_details_filtered$data)) {
      fwrite(results$top_details_filtered$data, file.path(table_dir, "stat_filtered_pairs_top_details.csv"))
    }
    if (!is.null(results$homoeolog_categories_swapped$data)) {
      fwrite(results$homoeolog_categories_swapped$data, file.path(table_dir, "stat_all_categories_swapped.csv"))
    }
    if (!is.null(results$scores_identity$scores_data)) {
      fwrite(results$scores_identity$scores_data, file.path(table_dir, "stat_scores_data.csv"))
    }
    if (!is.null(results$scores_identity$identity_data)) {
      fwrite(results$scores_identity$identity_data, file.path(table_dir, "stat_identity_data.csv"))
    }
    if (!is.null(results$pair_threshold_sweep$data)) {
      fwrite(results$pair_threshold_sweep$data, file.path(table_dir, "stat_pair_threshold_sweep.csv"))
    }
    if (!is.null(results$spacer_threshold_sweep$data)) {
      fwrite(results$spacer_threshold_sweep$data, file.path(table_dir, "stat_spacer_threshold_sweep.csv"))
    }
    if (!is.null(threshold_analysis$dataset_summary) && nrow(threshold_analysis$dataset_summary) > 0) {
      summary_dt <- copy(threshold_analysis$dataset_summary)
      summary_dt[, available_scores := vapply(
        available_scores,
        function(x) paste(x, collapse = ","),
        character(1)
      )]
      fwrite(summary_dt, file.path(table_dir, "stat_threshold_dataset_summary.csv"))
    }
    
    fwrite(full_data$pairs_grna_filtered, file.path(table_dir, "filtered_pairs_grna_data.csv"))
    fwrite(full_data$pair_classification, file.path(table_dir, "pair_classification_details.csv"))
    fwrite(full_data$pairs_dt, file.path(table_dir, "homoeologous_pairs.csv"))
    
    message("All CSV files saved to: ", table_dir)
  }
  
  # Generate summary report
  summary_file <- file.path(output_dir, "analysis_summary.txt")
  sink(summary_file)
  cat("gRNA ANALYSIS SUMMARY\n")
  cat("=====================\n\n")
  cat("Input files:\n")
  cat("- Classification file:", classification_file, "\n")
  cat("- gRNA database:", grna_db, "\n")
  cat("- Pairs file:", pairs_file, "\n\n")
  
  cat("Statistics:\n")
  cat("- Total gRNAs in classification:", scales::comma(nrow(basic_data)), "\n")
  cat("- gRNAs for homoeologous genes:", scales::comma(nrow(full_data$pairs_grna_filtered)), "\n")
  cat("- gRNAs failing CFD threshold:", scales::comma(length(full_data$bad_spacers)), "\n")
  cat("- Total homoeologous pairs:", scales::comma(nrow(full_data$pairs_dt)), "\n\n")
  
  cat("Pair classification:\n")
  if (!is.null(full_data$pair_counts)) {
    for (i in seq_len(nrow(full_data$pair_counts))) {
      cat("- ", full_data$pair_counts$pair_category[i], ": ",
          scales::comma(full_data$pair_counts$pair_count[i]), " (",
          full_data$pair_counts$percentage[i], "%)\n", sep = "")
    }
  }
  sink()
  
  message("\n=== ANALYSIS COMPLETE ===")
  message("Summary saved to: ", summary_file)
  if (save_plots) {
    message("All plots saved to: ", plot_dir)
  }
  if (save_csv) {
    message("All CSV files saved to: ", table_dir)
  }
  
  return(list(data = list(basic = basic_data, full = full_data), results = results))
}

# ---- COMMAND LINE INTERFACE -------------------------------------------------

if (!interactive()) {
  option_list <- list(
    make_option(
      c("-c", "--classification_file"), 
      type = "character",
      default = "gRNA_classification_final.csv",
      help = "Path to gRNA classification results CSV [default: %default]"
    ),
    make_option(
      c("-g", "--grna_db"), 
      type = "character",
      default = "Cbp_msk_all_genes_grna.filtered.article.tsv",
      help = "Path to gRNA database TSV with gene annotations [default: %default]"
    ),
    make_option(
      c("-p", "--pairs_file"), 
      type = "character",
      default = "pairs_perc_id_genes.txt",
      help = "Path to homoeologous gene pairs file (two columns) [default: %default]"
    ),
    make_option(
      c("-o", "--output_dir"), 
      type = "character",
      default = "gRNA_functionality_summary",
      help = "Output directory for tables and plots [default: %default]"
    ),
    make_option(
      "--cfd_threshold", 
      type = "double", 
      default = 0.05,
      help = "CFD off-target threshold used to filter gRNAs [default: %default]"
    ),
    make_option(
      "--plot_dpi", 
      type = "integer", 
      default = 300,
      help = "Resolution (dpi) for saved figures [default: %default]"
    ),
    make_option(
      "--no_csv", 
      action = "store_true", 
      default = FALSE,
      help = "Skip saving CSV files [default: %default]"
    )
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  results <- run_complete_analysis(
    classification_file = opt$classification_file,
    grna_db = opt$grna_db,
    pairs_file = opt$pairs_file,
    output_dir = opt$output_dir,
    cfd_threshold = opt$cfd_threshold,
    plot_dpi = opt$plot_dpi,
    save_plots = TRUE,
    save_csv = !opt$no_csv
  )
}

# ---- USAGE EXAMPLES ---------------------------------------------------------

# Example 1: Run complete analysis
# results <- run_complete_analysis(
#   classification_file = "data/gRNA_classification_final.csv",
#   grna_db = "data/grna_database.tsv",
#   pairs_file = "data/pairs.txt",
#   output_dir = "results/",
#   cfd_threshold = 0.05,
#   save_plots = TRUE,
#   save_csv = TRUE
# )

# Example 2: Generate individual plots interactively
# basic_data <- load_basic_data("data/gRNA_classification_final.csv")
# plot_global_category_combo(data = basic_data, save_plot = TRUE)
# plot_pam_status(data = basic_data, save_plot = TRUE)

# Example 3: Custom threshold analysis
# full_data <- load_full_data(
#   classification_file = "data/gRNA_classification_final.csv",
#   grna_db = "data/grna_database.tsv",
#   pairs_file = "data/pairs.txt"
# )
# threshold_analysis <- prepare_threshold_analysis(
#   pairs_dt = full_data$pairs_dt,
#   grna_all_dt = full_data$pairs_grna_all,
#   grna_filtered_dt = full_data$eligible_grnas,
#   thresholds = seq(0.1, 1.0, by = 0.05)
# )
# plot_pair_threshold_sweep(threshold_analysis, save_plot = TRUE)

# =============================================================================
# END OF SCRIPT
# =============================================================================
