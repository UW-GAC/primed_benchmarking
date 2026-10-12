## Offline / no-bucket use, and cleanup.
library(pgscompare)

# Skip all bucket syncing (work purely on local files):
# run_join_and_plot("ARIC", bucket_dir = NULL, sync_inputs = FALSE)

# Delete local scratch (downloaded score files + results incl. .rds caches).
# Bucket copies are untouched.
clear_pgs_local_files()
