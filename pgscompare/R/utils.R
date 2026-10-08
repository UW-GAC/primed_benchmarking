# Shared helpers -------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 ||
      (is.character(a) && length(a) == 1 && !is.na(a) && a == "")) b else a
}

#' Default results folder in the current Terra workspace bucket
#'
#' Terra sets the `WORKSPACE_BUCKET` environment variable inside notebooks and
#' RStudio. Returns `gs://<bucket>/<subdir>`, or `NULL` when it is not set
#' (i.e. not running on Terra), which every `bucket_dir` argument treats as
#' "skip bucket sync".
#' @param subdir folder under the bucket.
#' @export
default_bucket_dir <- function(subdir = "pgs_comparison_results") {
  b <- Sys.getenv("WORKSPACE_BUCKET", unset = "")
  if (!nzchar(b)) return(NULL)
  paste0(sub("/+$", "", b), "/", subdir)
}

# "<cohort>" or "<cohort>_<subset_tag>" (the prefix of every per-trait output).
name_prefix <- function(cohort_name, subset_tag = NULL) {
  if (is.null(subset_tag)) cohort_name else sprintf("%s_%s", cohort_name, subset_tag)
}

subset_tag_from_path <- function(subset_ids_path) {
  if (is.null(subset_ids_path)) NULL else tools::file_path_sans_ext(basename(subset_ids_path))
}

# Copy a gs:// file to local_dir (no-op for local paths / already downloaded).
localize_file <- function(gs_path, local_dir) {
  if (!startsWith(gs_path, "gs://")) return(gs_path)
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  local_path <- file.path(local_dir, basename(gs_path))
  if (!file.exists(local_path)) {
    status <- system2("gsutil", c("-q", "cp", shQuote(gs_path), shQuote(local_path)))
    if (status != 0 || !file.exists(local_path)) {
      stop(sprintf("gsutil cp failed for %s (exit status %d)", gs_path, status))
    }
  }
  local_path
}

#' Pull a bucket results folder down to a local folder
#' @param bucket_dir `gs://...` folder, or NULL to skip.
#' @param local_dir local destination.
#' @export
sync_from_bucket <- function(bucket_dir, local_dir) {
  if (is.null(bucket_dir)) return(invisible(NULL))
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("Syncing inputs from %s to %s ...\n", bucket_dir, local_dir))
  status <- system2("gsutil", c("-m", "cp", "-r",
                                shQuote(paste0(bucket_dir, "/*")), shQuote(local_dir)))
  if (status != 0) {
    warning("gsutil sync from bucket reported a non-zero exit status - check that files ",
            "exist under ", bucket_dir, " (e.g. `gsutil ls ", bucket_dir, "`).")
  }
  invisible(status)
}

#' Push CSVs and/or the plots folder to the bucket
#' @param bucket_dir `gs://...` folder, or NULL to skip.
#' @param local_dir local results folder.
#' @param csv_glob glob of CSVs inside `local_dir` to push (NULL = none), e.g.
#'   `"GeneSTAR*.csv"`.
#' @param plots push `<local_dir>/plots` too?
#' @export
sync_to_bucket <- function(bucket_dir, local_dir, csv_glob = "*.csv", plots = TRUE) {
  if (is.null(bucket_dir)) return(invisible(NULL))
  bad <- FALSE
  if (!is.null(csv_glob)) {
    cat(sprintf("Syncing CSVs (%s) from %s to %s ...\n", csv_glob, local_dir, bucket_dir))
    s <- system2("gsutil", c("-m", "cp", shQuote(file.path(local_dir, csv_glob)),
                             shQuote(paste0(bucket_dir, "/"))))
    bad <- bad || s != 0
  }
  plots_root <- file.path(local_dir, "plots")
  if (plots && dir.exists(plots_root)) {
    cat(sprintf("Syncing %s to %s ...\n", plots_root, bucket_dir))
    s <- system2("gsutil", c("-m", "cp", "-r", shQuote(plots_root),
                             shQuote(paste0(bucket_dir, "/"))))
    bad <- bad || s != 0
  }
  if (bad) warning("gsutil sync to bucket reported a non-zero exit status - check the output above.")
  else cat("Bucket sync complete.\n")
  invisible(!bad)
}

#' Delete the local scratch files the pipeline creates
#'
#' Removes the downloaded score files and the results folder (including the
#' per-chunk .rds caches). Anything already pushed to the bucket is untouched.
#' @param local_dir,out_dir the folders used in [run_combine_and_compare()].
#' @param confirm ask before deleting (interactive sessions only).
#' @export
clear_pgs_local_files <- function(local_dir = "pgs_sscore_files",
                                  out_dir = "pgs_comparison_results", confirm = TRUE) {
  targets <- c(local_dir, out_dir)
  targets <- targets[dir.exists(targets)]
  if (length(targets) == 0) { cat("Nothing to delete.\n"); return(invisible(NULL)) }
  if (confirm && interactive()) {
    ans <- readline(sprintf("Delete %s? [y/N] ", paste(targets, collapse = ", ")))
    if (!tolower(ans) %in% c("y", "yes")) { cat("Cancelled.\n"); return(invisible(NULL)) }
  }
  unlink(targets, recursive = TRUE)
  cat("Deleted:", paste(targets, collapse = ", "), "\n")
  invisible(targets)
}

# avtable() lives in AnVILGCP on newer installs, AnVIL on older ones.
resolve_avtable <- function() {
  for (pkg in c("AnVILGCP", "AnVIL")) {
    if (requireNamespace(pkg, quietly = TRUE) &&
        "avtable" %in% getNamespaceExports(pkg)) {
      return(getExportedValue(pkg, "avtable"))
    }
  }
  stop(
    "Could not find avtable() in AnVIL or AnVILGCP.\n",
    "Restart R, then in that CLEAN session run:\n",
    "    BiocManager::install(c('AnVIL', 'AnVILGCP'), update = FALSE, ask = FALSE)\n",
    "(Installing mid-session after other packages are loaded can fail with a\n",
    "namespace-unload error.) Alternatively pass `pgs_table =` a data.frame."
  )
}

chunk_vector <- function(x, size) unname(split(x, ceiling(seq_along(x) / size)))

# Which type variant (full/partial) a metric's columns are read from.
type_for_metric <- function(metric_name, opts) {
  t <- opts$metric_types[[metric_name]]
  if (is.null(t)) opts$default_type else t
}

metric_col <- function(metric_name, stat, opts) {
  sprintf("%s_%s_%s", metric_name, type_for_metric(metric_name, opts), stat)
}
metric_col_explicit <- function(metric_name, type_name, stat) {
  sprintf("%s_%s_%s", metric_name, type_name, stat)
}

# coord_equal() only when both axes cover similar ranges, else free axes.
axis_coord <- function(x_vals, y_vals, max_axis_ratio = 1.5) {
  x_range <- diff(range(x_vals, na.rm = TRUE))
  y_range <- diff(range(y_vals, na.rm = TRUE))
  if (!is.finite(x_range) || !is.finite(y_range) || min(x_range, y_range) <= 0) {
    return(ggplot2::coord_cartesian())
  }
  if (max(x_range, y_range) / min(x_range, y_range) <= max_axis_ratio) {
    ggplot2::coord_equal()
  } else {
    ggplot2::coord_cartesian()
  }
}

# "<cohort>_<trait>_<description>.csv" first, then "<cohort>_<description>.csv".
resolve_optional <- function(results_dir, cohort_name, description, trait_name) {
  p <- file.path(results_dir, sprintf("%s_%s_%s.csv", cohort_name, trait_name, description))
  if (file.exists(p)) return(p)
  p <- file.path(results_dir, sprintf("%s_%s.csv", cohort_name, description))
  if (file.exists(p)) return(p)
  NA_character_
}
resolve_trait_or_global <- function(results_dir, trait_specific_name, global_name) {
  p <- file.path(results_dir, trait_specific_name)
  if (file.exists(p)) return(p)
  p <- file.path(results_dir, global_name)
  if (file.exists(p)) return(p)
  NA_character_
}

save_plot <- function(p, path, width, height, ...) {
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 150, ...)
  cat(sprintf("  saved %s\n", path))
}

sanitize_for_filename <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# Traits available for a run. File names alone are ambiguous (the no-subset run's
# pattern "<cohort>_*_pgsmetrics.csv" also matches "<cohort>_<tag>_<trait>_..."),
# so use the subset_tag recorded in the "_sample_size.csv" files when present and
# fall back to the filename pattern otherwise.
detect_traits <- function(results_dir, cohort_name, tag = NULL) {
  prefix <- name_prefix(cohort_name, tag)
  ss_pattern <- paste0("^", prefix, "_(.*)_sample_size\\.csv$")
  ss_files <- list.files(results_dir, pattern = ss_pattern, full.names = TRUE)
  traits <- character(0)
  for (f in ss_files) {
    ss <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(ss) || !all(c("trait", "subset_tag") %in% names(ss))) next
    file_tag <- ss$subset_tag[1]
    if ((is.null(tag) && is.na(file_tag)) || (!is.null(tag) && identical(as.character(file_tag), tag))) {
      traits <- c(traits, as.character(ss$trait[1]))
    }
  }
  if (length(traits) == 0) {
    pattern <- paste0("^", prefix, "_(.*)_pgsmetrics\\.csv$")
    traits <- sub(pattern, "\\1", list.files(results_dir, pattern = pattern))
  }
  traits[file.exists(file.path(results_dir, sprintf("%s_%s_pgsmetrics.csv", prefix, traits)))]
}
