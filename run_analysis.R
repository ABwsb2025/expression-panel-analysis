#!/usr/bin/env Rscript
# =====================================================================
# Command-line entry point.
#
# Usage:
#   Rscript run_analysis.R --matrix <expression_matrix> \
#                          --panel <gene_panel.tsv> \
#                          --control-pattern "<regex>" \
#                          [--control-label Control] \
#                          [--treatment-label Treatment] \
#                          [--welch] [--outdir results]
# =====================================================================

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) == 0) return(".")
  dirname(normalizePath(sub("^--file=", "", arg)))
}
source(file.path(get_script_dir(), "R", "panel_analysis.R"))


# ---- Parse arguments -------------------------------------------------
parse_args <- function(args) {
  opts <- list(control_label = "Control", treatment_label = "Treatment",
               welch = FALSE, outdir = "results")
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--welch") {
      opts$welch <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for ", key)
    val <- args[i + 1]
    switch(key,
      "--matrix"          = opts$matrix <- val,
      "--panel"           = opts$panel <- val,
      "--control-pattern" = opts$control_pattern <- val,
      "--control-label"   = opts$control_label <- val,
      "--treatment-label" = opts$treatment_label <- val,
      "--outdir"          = opts$outdir <- val,
      stop("Unknown argument: ", key)
    )
    i <- i + 2
  }
  for (req in c("matrix", "panel", "control_pattern")) {
    if (is.null(opts[[req]])) stop("Required: --", gsub("_", "-", req))
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))


# ---- Run -------------------------------------------------------------
mat     <- read_expression_matrix(opts$matrix)
panel   <- read_gene_panel(opts$panel)
long_df <- to_long_panel(mat, panel, opts$control_pattern,
                         opts$control_label, opts$treatment_label)

counts <- table(long_df$group) / length(unique(long_df$gene_id))
message("Samples per group: ",
        paste(names(counts), counts, sep = " = ", collapse = ", "))
if (any(counts < 2)) stop("Each group needs at least 2 samples. Check --control-pattern.")

tests   <- test_panel(long_df, welch = opts$welch)
summary <- summarise_panel(long_df, tests)
plot    <- plot_panel(long_df, tests)

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)
write.table(summary, file.path(opts$outdir, "panel_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ggsave(file.path(opts$outdir, "panel_barplot.png"), plot,
       width = 2.7 * min(length(unique(long_df$gene_id)), 4) + 0.5,
       height = 3.5 * ceiling(length(unique(long_df$gene_id)) / 4), dpi = 300)

print(tests)
message("Wrote results to ", normalizePath(opts$outdir))
