## Run once per workspace, in a FRESH R session (restart the kernel first).
install.packages("remotes")
remotes::install_github("UW-GAC/primed_benchmarking", subdir = "pgscompare")
remotes::install_github("krillinor/pgsmetrics", upgrade = "never")  # not on CRAN
