## Step 4: compare two runs (e.g. whole cohort vs a sub-population).
## Needs both runs from step 1 (tag NULL = the no-subset run).
library(pgscompare)

run_compare_populations(
  "ARIC",
  population_comparisons = list(
    list(label1 = "Population",          tag1 = NULL,
         label2 = "K1K4 sub-population", tag2 = "k1k4")
  ),
  color_var = "gwas_ancestry"                     # "dev_ancestry", or NULL for no color
)
