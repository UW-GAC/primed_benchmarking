#' Combine the per-PGS variant-match summary files for a cohort
#'
#' R port of `inst/python/combine_pgs_summaries.py`. Reads the Terra data table
#' of scoring files (default `pgs_scoring_file`), keeps the rows whose
#' `match_summary` path ends in `<cohort_name>_summary.csv`, reads each summary
#' file from the bucket, and stacks them. Writes (to `out_dir`):
#' * `<cohort>_summary_combined.csv`
#' * `<cohort>_summary_<status>_flipped.csv` and `..._not_flipped.csv`
#' * `<cohort>_match_rate.csv` (if `write_match_rate`): one row per PGS with
#'   `cohort, pgs_model_id, total_count, count, percent, pgs_scoring_file_id`
#'   (`count` = variants with `match_status == match_status_filter`,
#'   `percent` = 100 * count / total_count). [run_join_and_plot()] picks this up
#'   for the match-rate plots.
#'
#' @param cohort_name cohort name, e.g. "GeneSTAR".
#' @param entity_type Terra data table with the scoring files.
#' @param summary_col column holding the `gs://` path of each summary file.
#' @param id_col column with the PGS id (matches the file's `accession`).
#' @param file_id_col column with the scoring-file entity id; NULL = the
#'   table's first column (Terra's `<table>_id`).
#' @param match_status_filter status used for the flipped / not-flipped split.
#' @param out_dir local output folder.
#' @param bucket_dir where the combined/flipped CSVs are pushed (NULL = skip).
#'   Defaults to `<WORKSPACE_BUCKET>/pgs_combined_summaries`.
#' @param write_match_rate also write `<cohort>_match_rate.csv`.
#' @param match_rate_bucket_dir where that file is pushed; defaults to the
#'   shared `pgs_comparison_results` folder that [run_join_and_plot()] reads.
#' @param scoring_table optional data.frame to use instead of the Terra table
#'   (handy for testing; paths may then be local files).
#' @return Invisibly, a list with `combined`, `matched_flipped`,
#'   `matched_not_flipped`, `match_rate` and `failures`.
#' @export
combine_pgs_match_summaries <- function(
    cohort_name, entity_type = "pgs_scoring_file", summary_col = "match_summary",
    id_col = "pgs_model_id", file_id_col = NULL, match_status_filter = "matched",
    out_dir = "pgs_comparison_results",
    bucket_dir = default_bucket_dir("pgs_combined_summaries"),
    write_match_rate = TRUE, match_rate_bucket_dir = default_bucket_dir(),
    scoring_table = NULL) {

  summary_filename <- paste0(cohort_name, "_summary.csv")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  full <- if (is.null(scoring_table)) as.data.table(resolve_avtable()(entity_type)) else as.data.table(scoring_table)
  cat(sprintf("Pulled %d rows from data table '%s'\n", nrow(full), entity_type))
  for (col in c(summary_col, id_col)) {
    if (!(col %in% names(full))) {
      stop(sprintf("Column '%s' not found in '%s'. Columns present: %s",
                   col, entity_type, paste(names(full), collapse = ", ")))
    }
  }
  if (is.null(file_id_col)) file_id_col <- names(full)[1]

  # a shared table can span cohorts: keep rows whose summary file is this cohort's
  tbl <- full[which(endsWith(as.character(full[[summary_col]]), summary_filename))]
  cat(sprintf("Filtered to %d rows matching '%s'\n", nrow(tbl), summary_filename))
  if (nrow(tbl) == 0) {
    stop(sprintf("No rows found with %s ending in '%s'. Check cohort_name, or inspect the '%s' column to see the naming convention this workspace uses.",
                 summary_col, summary_filename, summary_col))
  }

  frames <- list(); failures <- list()
  for (i in seq_len(nrow(tbl))) {
    pgs_id  <- as.character(tbl[[id_col]][i])
    gs_path <- as.character(tbl[[summary_col]][i])
    df <- tryCatch(read_gs_csv(gs_path), error = function(e) e)
    if (inherits(df, "error")) {
      failures[[length(failures) + 1]] <- list(pgs_id = pgs_id, path = gs_path, error = conditionMessage(df))
      next
    }
    if ("accession" %in% names(df) && !all(df$accession == pgs_id)) {
      cat(sprintf("WARNING: accession mismatch in file for %s\n", pgs_id))
    }
    df[, pgs_scoring_file_id := as.character(tbl[[file_id_col]][i])]
    frames[[length(frames) + 1]] <- df
  }
  if (length(failures) > 0) {
    cat(sprintf("\n%d file(s) failed to download:\n", length(failures)))
    for (f in failures) cat(sprintf("  %s: %s\n    %s\n", f$pgs_id, f$path, f$error))
  }
  if (length(frames) == 0) stop("No summary files could be read.")

  combined <- rbindlist(frames, use.names = TRUE, fill = TRUE)
  cat(sprintf("\nCombined %d files -> %d rows, %d unique PGS accessions.\n",
              length(frames), nrow(combined), uniqueN(combined$accession)))
  f_combined <- file.path(out_dir, sprintf("%s_summary_combined.csv", cohort_name))
  fwrite(combined, f_combined)
  cat(sprintf("Saved combined table to %s\n", f_combined))

  cat("\nVariant counts by match_status across all PGS scores:\n")
  print(combined[, list(count = sum(count)), by = match_status])

  combined[, total_count := sum(count), by = accession]
  # match_flipped is blank for unmatched rows -> real logical (NA when blank)
  combined[, match_flipped_bool := unname(c(true = TRUE, false = FALSE)[tolower(as.character(match_flipped))])]

  pick <- function(flipped) combined[which(combined$match_status == match_status_filter &
                                             combined$match_flipped_bool %in% flipped)]
  matched_flipped     <- pick(TRUE)
  matched_not_flipped <- pick(FALSE)
  cat(sprintf("\n%s + flipped: %d rows, %s total variants across %d PGS scores\n",
              match_status_filter, nrow(matched_flipped), format(sum(matched_flipped$count), big.mark = ","),
              uniqueN(matched_flipped$accession)))
  cat(sprintf("%s + not flipped: %d rows, %s total variants\n",
              match_status_filter, nrow(matched_not_flipped), format(sum(matched_not_flipped$count), big.mark = ",")))
  f_flip <- file.path(out_dir, sprintf("%s_summary_%s_flipped.csv", cohort_name, match_status_filter))
  f_nof  <- file.path(out_dir, sprintf("%s_summary_%s_not_flipped.csv", cohort_name, match_status_filter))
  fwrite(matched_flipped, f_flip); fwrite(matched_not_flipped, f_nof)

  match_rate <- NULL
  if (write_match_rate) {
    match_rate <- combined[, list(
      total_count = sum(count),
      count       = sum(count[match_status == match_status_filter]),
      pgs_scoring_file_id = pgs_scoring_file_id[1]), by = list(pgs_model_id = accession)]
    match_rate[, percent := 100 * count / total_count]
    match_rate[, cohort := cohort_name]
    setcolorder(match_rate, c("cohort", "pgs_model_id", "total_count", "count", "percent", "pgs_scoring_file_id"))
    f_mr <- file.path(out_dir, sprintf("%s_match_rate.csv", cohort_name))
    fwrite(match_rate, f_mr)
    cat(sprintf("Wrote match rates for %d PGS to %s\n", nrow(match_rate), f_mr))
  }

  if (!is.null(bucket_dir)) {
    for (f in c(f_combined, f_flip, f_nof)) {
      system2("gsutil", c("cp", shQuote(f), shQuote(paste0(bucket_dir, "/"))))
    }
  }
  if (write_match_rate && !is.null(match_rate_bucket_dir)) {
    system2("gsutil", c("cp", shQuote(f_mr), shQuote(paste0(match_rate_bucket_dir, "/"))))
  }
  invisible(list(combined = combined, matched_flipped = matched_flipped,
                 matched_not_flipped = matched_not_flipped, match_rate = match_rate,
                 failures = failures))
}

# read a CSV from gs:// (via gsutil cat) or a local path
read_gs_csv <- function(path) {
  if (!startsWith(path, "gs://")) return(fread(path))
  out <- suppressWarnings(system2("gsutil", c("cat", shQuote(path)), stdout = TRUE, stderr = FALSE))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop(sprintf("gsutil cat failed for %s (exit status %d)", path, st))
  fread(text = paste(out, collapse = "\n"))
}
