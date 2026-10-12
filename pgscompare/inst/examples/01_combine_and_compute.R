## Step 1: combine the per-PGS score files and compute metrics for ONE trait.
## Reads the Terra data table `pgs_individual_file`; writes CSVs to
## pgs_comparison_results/ and pushes them to $WORKSPACE_BUCKET/pgs_comparison_results.
library(pgscompare)

COHORT <- "ARIC"
BUCKET <- Sys.getenv("WORKSPACE_BUCKET")          # gs://fc-... set by Terra

# --- full cohort ---
run_combine_and_compare(
  cohort_name   = COHORT,
  pheno_path    = file.path(BUCKET, "phenotypes/aric_pheno.tsv"),
  pheno_id_col  = "IID",
  covariates    = c("age", "sex", paste0("PC", 1:5)),
  binary_traits = c(t2d = "t2d_status"),          # output name = phenotype column
  prevalence    = c(t2d = 0.10),
  # table cleanup (all optional)
  keep_samplesets = COHORT, sampleset_col = "sampleset",
  exclude_pgs_ids = "PGS002239_Top500",
  n_boot = 1000, n_cores = 4
)

# --- same cohort restricted to a subset of people ---
# Output is named ARIC_k1k4_t2d_pgsmetrics.csv (tag = file name without extension).
run_combine_and_compute_subset <- function() {
  run_combine_and_compare(
    cohort_name   = COHORT,
    pheno_path    = file.path(BUCKET, "phenotypes/aric_pheno.tsv"),
    covariates    = c("age", "sex", paste0("PC", 1:5)),
    binary_traits = c(t2d = "t2d_status"), prevalence = c(t2d = 0.10),
    keep_samplesets = COHORT, exclude_pgs_ids = "PGS002239_Top500",
    subset_ids_path = file.path(BUCKET, "subsets/k1k4.txt")
  )
}
# run_combine_and_compute_subset()   # uncomment to run

# A continuous trait works the same way:
# run_combine_and_compare(cohort_name = COHORT, pheno_path = ..., covariates = ...,
#                         continuous_traits = c(bmi = "bmi"))
