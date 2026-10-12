## Step 2: PGS Catalog metadata + reported performance for the trait's PGS.
## Writes t2d_pgs_catalog_score_metadata.csv and t2d_pgs_catalog_performance_metrics.csv.
library(pgscompare)

fetch_pgs_catalog_metadata(
  pgs_metadata_path = "pgs_comparison_results/ARIC_pgs_metadata.csv",  # from step 1
  trait = "t2d"
)
