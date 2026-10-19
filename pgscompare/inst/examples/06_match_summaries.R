## Combine the per-PGS variant-match summaries for a cohort + trait (R port of the
## combine_pgs_summaries notebook in inst/python/) and write
## <cohort>_<trait>_match_rate.csv, which run_join_and_plot() uses for the
## match-rate plots. Both outputs are pushed to $WORKSPACE_BUCKET/pgs_comparison_results.
library(pgscompare)

combine_pgs_match_summaries("ARIC", sub_name = "t2d")   # reads Terra table `pgs_scoring_file_t2d`
# Writes ARIC_t2d_summary_combined.csv and ARIC_t2d_match_rate.csv
# (plus ARIC_t2d_summary_matched_[not_]flipped.csv, kept local).

# No sub_name: table `pgs_scoring_file`, files named ARIC_summary_combined.csv / ARIC_match_rate.csv
# combine_pgs_match_summaries("ARIC")
