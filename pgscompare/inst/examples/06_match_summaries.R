## Combine the per-PGS variant-match summaries for a cohort (R port of
## inst/python/combine_pgs_summaries.py) and write <cohort>_match_rate.csv,
## which run_join_and_plot() uses for the match-rate plots.
library(pgscompare)

combine_pgs_match_summaries("GeneSTAR")     # reads Terra table `pgs_scoring_file`
# Writes GeneSTAR_summary_combined.csv, ..._matched_flipped.csv, ..._matched_not_flipped.csv
# and GeneSTAR_match_rate.csv (pushed to $WORKSPACE_BUCKET/pgs_comparison_results).
