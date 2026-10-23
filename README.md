# Two-group expression analysis for a gene panel

A small R pipeline for comparing expression of a chosen set of genes between two conditions (for example control vs knockdown) in an RNA-seq expression matrix. It handles input cleaning, per-gene statistical testing with multiple-testing correction, and publication-style plotting.

I wrote the original version during a summer research placement at UCL, then generalised it into a reusable tool. No data or results from that project are included here.

## What it does

1. **Loads** a genes × samples expression matrix (tab, comma or whitespace delimited). It also repairs files where tabs were exported as the literal string `/t` or `\t`.
2. **Normalises gene IDs** by stripping Ensembl version suffixes (`ENSG00000000001.12` → `ENSG00000000001`), so IDs match across annotation releases.
3. **Subsets** the matrix to a user-supplied gene panel and reshapes it to long format, assigning each sample to a group by regex on its column name.
4. **Tests** each gene with a two-tailed t-test (Student's by default, Welch's with `--welch`) and applies Benjamini-Hochberg correction across the panel.
5. **Summarises** mean, SD and replicate count per gene and group alongside the test statistics.
6. **Plots** a faceted bar chart with SD error bars, individual replicates overlaid, and a significance bracket per gene positioned automatically above the highest point.

## Requirements

R ≥ 4.1 with `dplyr`, `tidyr`, `tibble` and `ggplot2`:

```r
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2"))
```

## Usage

```bash
Rscript run_analysis.R \
  --matrix path/to/expression_matrix.tsv \
  --panel path/to/gene_panel.tsv \
  --control-pattern "^ctrl" \
  --control-label Control \
  --treatment-label Knockdown \
  --outdir results
```

| Argument | Required | Description |
|---|---|---|
| `--matrix` | yes | Expression matrix, genes as rows, samples as columns |
| `--panel` | yes | Gene panel file (see below) |
| `--control-pattern` | yes | Regex matching control sample names; all other samples form the second group |
| `--control-label` | no | Display name for the control group (default `Control`) |
| `--treatment-label` | no | Display name for the other group (default `Treatment`) |
| `--welch` | no | Use Welch's t-test instead of Student's |
| `--outdir` | no | Output directory (default `results`) |

### Input formats

**Expression matrix**: first column is the gene ID, header row is sample names.

```
gene_id          ctrl_1  ctrl_2  ctrl_3  kd_1  kd_2  kd_3
ENSG00000000001  ...     ...     ...     ...   ...   ...
```

**Gene panel**: tab-separated with columns `gene_id` and `gene_name`.

```
gene_id          gene_name
ENSG00000000001  GENE_A
ENSG00000000002  GENE_B
```

### Outputs

`panel_summary.tsv` contains one row per gene and group with mean, SD, n, t statistic, degrees of freedom, raw and adjusted p-values and a significance label. `panel_barplot.png` is the faceted plot.

Significance labels use adjusted p-values: `****` < 0.0001, `***` < 0.001, `**` < 0.01, `*` < 0.05.

## Using the functions directly

```r
source("R/panel_analysis.R")

mat     <- read_expression_matrix("expression_matrix.tsv")
panel   <- read_gene_panel("gene_panel.tsv")
long_df <- to_long_panel(mat, panel, control_pattern = "^ctrl")
tests   <- test_panel(long_df)
plot_panel(long_df, tests)
```

## Notes and limitations

This is a targeted test on a small, pre-selected panel, not a genome-wide differential expression analysis. For transcriptome-wide work on raw counts, use a count-based model such as DESeq2 or edgeR. With few replicates per group, t-test results should be read alongside effect size and replicate spread, which is why the plot shows every replicate.

## Licence

MIT
