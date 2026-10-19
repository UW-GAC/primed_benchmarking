# pgscompare

Combine many per-PGS score files on AnVIL/Terra, evaluate them with
[pgsmetrics](https://krillinor.github.io/pgsmetrics/), and compare them (against
PGS Catalog metadata, adjusted vs. raw, and across populations/subsets).

## Install (in each AnVIL workspace)

Run this in a **fresh R session** (restart the kernel first; installing after other
packages are loaded can trigger the "cannot be unloaded" error):

```r
install.packages("remotes")
remotes::install_github("UW-GAC/primed_benchmarking", subdir = "pgscompare")

# pgsmetrics is not on CRAN; install it separately (and keep upgrade = "never"):
remotes::install_github("krillinor/pgsmetrics", upgrade = "never")
```

Private repo? Create a GitHub token and use
`remotes::install_github("UW-GAC/primed_benchmarking", subdir = "pgscompare", auth_token = Sys.getenv("GITHUB_PAT"))`.
To pick up changes later, reinstall (add `force = TRUE`) and restart R.
`AnVIL`/`AnVILGCP` are only needed to read the Terra data table; they are already on
Terra images (or pass `pgs_table =` a data.frame).

## Workflow

`bucket_dir` defaults to `$WORKSPACE_BUCKET/pgs_comparison_results` on Terra, so nothing
bucket-specific is hard-coded. Pass `bucket_dir = NULL` to disable syncing.

```r
library(pgscompare)

# 1. Combine score files + compute metrics (one trait at a time is fine)
run_combine_and_compare(
  cohort_name   = "ARIC",
  pheno_path    = "gs://<bucket>/phenotypes.tsv",
  covariates    = c("age", "sex", paste0("PC", 1:5)),
  binary_traits = c(t2d = "t2d_status"),
  prevalence    = c(t2d = 0.1),
  keep_samplesets = "ARIC", sampleset_col = "sampleset",
  exclude_pgs_ids = "PGS002239_Top500",
  # subset_ids_path = "gs://<bucket>/k1k4_ids.txt",   # -> ARIC_k1k4_t2d_pgsmetrics.csv
  n_boot = 1000, n_cores = 4
)

# 2. PGS Catalog metadata (writes t2d_pgs_catalog_score_metadata.csv, t2d_pgs_catalog_performance_metrics.csv)
fetch_pgs_catalog_metadata(
  pgs_metadata_path = "pgs_comparison_results/ARIC_pgs_metadata.csv", trait = "t2d")

# 2b. Variant match rates (optional; feeds the match-rate plots)
combine_pgs_match_summaries("ARIC", sub_name = "t2d")   # -> ARIC_t2d_match_rate.csv

# 3. Join + all plots (add match rate file ARIC_t2d_match_rate.csv to the results folder/bucket)
run_join_and_plot("ARIC")                          # full cohort
run_join_and_plot("ARIC", subset_tag = "k1k4")     # a subset run

# 4. Population vs sub-population
run_compare_populations(
  "ARIC",
  list(list(label1 = "Population", tag1 = NULL,
            label2 = "K1K4 sub-population", tag2 = "k1k4")),
  color_var = "gwas_ancestry")                     # or "dev_ancestry" / NULL

# Housekeeping: delete local scratch (bucket copies untouched)
clear_pgs_local_files()
```

Every function documents its arguments (`?run_combine_and_compare`). Defaults match the
original scripts; the old global `CONFIG` list is now just function arguments, so a stale
copy in a notebook can no longer drift from the code.

Ready-to-run versions of these steps are in `inst/examples/` (`00_install.R` … `06_match_summaries.R`). `inst/python/combine_pgs_summaries.py` is the original Python version of the match-summary step, kept for reference.

## Output files (in `out_dir`)

`<cohort>[_<subset_tag>]_<trait>_pgsmetrics.csv`, `..._pgsmetrics_cov.csv`,
`..._sample_size.csv`, `..._or_predictions.csv`; `<cohort>_combined_pgs_wide.csv`,
`<cohort>_pgs_metadata.csv`; plots under `plots/<cohort>[/<subset_tag>]/<trait>/` and
`plots/<cohort>/population_comparisons/<trait>/`.

## Where it lives

This package is the `pgscompare/` folder of
[UW-GAC/primed_benchmarking](https://github.com/UW-GAC/primed_benchmarking), hence the
`subdir = "pgscompare"` in the install command. To install a branch, add `ref = "<branch>"`.
