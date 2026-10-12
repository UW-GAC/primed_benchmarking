## Step 3: join everything and make the plots.
## Pulls inputs from the bucket, saves plots under
## pgs_comparison_results/plots/<cohort>[/<subset>]/<trait>/, pushes them back.
## Optional: put ARIC_t2d_match_rate.csv (cols pgs_model_id, percent, total_count, count)
## in the results folder / bucket to get the match-rate plots.
library(pgscompare)

run_join_and_plot("ARIC")                         # full cohort
run_join_and_plot("ARIC", subset_tag = "k1k4")    # subset run

# Narrow things down / customise:
run_join_and_plot(
  "ARIC",
  traits          = "t2d",
  metrics_to_plot = list(t2d = c("roc_auc", "r2_nagelkerke", "or")),
  default_type    = "partial",                    # or "full"
  metric_types    = list(or = "full"),
  metric_comparisons = list(
    list(x = "roc_auc", y = "r2_nagelkerke", color = "dev_ancestry"),
    list(x = "roc_auc", y = "or",            color = "gwas_ancestry")
  )
)
