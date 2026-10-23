# =====================================================================
# Two-group expression analysis for a user-defined gene panel
#
# Functions for loading an expression matrix, subsetting it to a panel
# of genes, testing each gene between two conditions, and plotting the
# result. See run_analysis.R for a command-line entry point.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})


# ---- Loading ---------------------------------------------------------

#' Read an expression matrix (genes as rows, samples as columns).
#'
#' The first column must hold gene IDs and the header must hold sample
#' names. Tab, comma and whitespace delimiters are handled. Some exported
#' tables contain the literal two-character string "/t" or "\t" instead
#' of real tab characters, so those are converted before parsing.
read_expression_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- gsub("/t", "\t", lines, fixed = TRUE)
  lines <- gsub("\\t", "\t", lines, fixed = TRUE)

  sep <- if (grepl("\t", lines[1])) "\t" else if (grepl(",", lines[1])) "," else ""

  mat <- read.table(text = lines, header = TRUE, sep = sep,
                    check.names = FALSE, stringsAsFactors = FALSE)

  # If the header names the ID column, read.table keeps it as a normal
  # column; move it into the row names. (If the header has one fewer
  # field than the data rows, read.table has already done this.)
  if (is.character(mat[[1]])) {
    rownames(mat) <- mat[[1]]
    mat <- mat[, -1, drop = FALSE]
  }

  if (!all(vapply(mat, is.numeric, logical(1)))) {
    stop("All sample columns must be numeric. Check the delimiter and header.")
  }
  mat
}


#' Remove Ensembl version suffixes, e.g. ENSG00000000001.12 -> ENSG00000000001.
strip_ensembl_version <- function(ids) {
  sub("\\.[0-9]+$", "", ids)
}


#' Read a gene panel: a two-column tab-separated file with a header
#' row, columns gene_id and gene_name.
read_gene_panel <- function(path) {
  panel <- read.delim(path, stringsAsFactors = FALSE)
  required <- c("gene_id", "gene_name")
  if (!all(required %in% names(panel))) {
    stop("Gene panel needs columns: ", paste(required, collapse = ", "))
  }
  panel$gene_id <- strip_ensembl_version(panel$gene_id)
  panel
}


# ---- Reshaping -------------------------------------------------------

#' Subset the matrix to the panel and convert to long format, assigning
#' each sample to a group by matching its name against a regex.
to_long_panel <- function(mat, panel, control_pattern,
                          control_label = "Control",
                          treatment_label = "Treatment") {
  rownames(mat) <- strip_ensembl_version(rownames(mat))

  missing <- setdiff(panel$gene_id, rownames(mat))
  if (length(missing) > 0) {
    warning("Not found in matrix: ", paste(missing, collapse = ", "))
  }

  mat[rownames(mat) %in% panel$gene_id, , drop = FALSE] |>
    rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "sample", values_to = "value") |>
    left_join(panel, by = "gene_id") |>
    mutate(
      group = ifelse(grepl(control_pattern, sample), control_label, treatment_label),
      group = factor(group, levels = c(control_label, treatment_label))
    )
}


# ---- Statistics ------------------------------------------------------

#' Significance label. Conventional thresholds:
#' **** p < 0.0001, *** p < 0.001, ** p < 0.01, * p < 0.05.
star_label <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "n.s."
  )
}


#' Two-tailed t-test per gene. Student's by default; set welch = TRUE
#' if group variances may differ. Benjamini-Hochberg adjusted p-values
#' are added because several genes are tested at once.
test_panel <- function(long_df, welch = FALSE) {
  long_df |>
    group_by(gene_id, gene_name) |>
    group_modify(~ {
      tt <- t.test(value ~ group, data = .x, var.equal = !welch)
      tibble(t_stat = unname(tt$statistic),
             df     = unname(tt$parameter),
             p_value = tt$p.value)
    }) |>
    ungroup() |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           significance = star_label(p_adj))
}


#' Per-gene, per-group summary (mean, SD, n) joined to the test results.
summarise_panel <- function(long_df, tests) {
  long_df |>
    group_by(gene_id, gene_name, group) |>
    summarise(mean = mean(value), SD = sd(value), n = n(), .groups = "drop") |>
    left_join(tests, by = c("gene_id", "gene_name"))
}


# ---- Plotting --------------------------------------------------------

#' Faceted bar plot: group means with SD error bars, individual
#' replicates jittered on top, and a significance bracket per gene
#' placed just above the highest point in that facet.
plot_panel <- function(long_df, tests, y_label = "Expression") {
  annot <- long_df |>
    group_by(gene_name, group) |>
    summarise(top = max(max(value), mean(value) + sd(value)), .groups = "drop") |>
    group_by(gene_name) |>
    summarise(y_max = max(top), .groups = "drop") |>
    left_join(select(tests, gene_name, significance), by = "gene_name") |>
    mutate(y_bracket = y_max * 1.08, y_text = y_max * 1.11)

  ggplot(long_df, aes(x = group, y = value, fill = group)) +
    stat_summary(fun = mean, geom = "bar", width = 0.6) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    geom_segment(data = annot,
                 aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = annot, aes(x = 1.5, y = y_text, label = significance),
              inherit.aes = FALSE, size = 4, vjust = 0) +
    facet_wrap(~ gene_name, scales = "free_y",
               ncol = min(length(unique(long_df$gene_name)), 4)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = y_label) +
    theme_bw() +
    theme(legend.position = "none")
}
