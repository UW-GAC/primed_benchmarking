#' Fetch PGS Catalog score metadata and reported performance metrics
#'
#' Looks up each PGS id in the PGS Catalog REST API and writes
#' `[<trait>_]pgs_catalog_score_metadata.csv` (one row per PGS: method, number of
#' variants, GWAS and development ancestry/sample size) and
#' `[<trait>_]pgs_catalog_performance_metrics.csv` (one row per PGS x
#' evaluation x metric).
#'
#' @param pgs_ids character vector of PGS Catalog ids. If NULL they are read
#'   from `pgs_metadata_path`.
#' @param pgs_metadata_path CSV with a `pgs_id` column, typically
#'   `<out_dir>/<cohort_name>_pgs_metadata.csv` written by
#'   [run_combine_and_compare()].
#' @param trait optional trait name; when given, output files are prefixed
#'   `<trait>_` (the naming [run_join_and_plot()] looks for first).
#' @param out_dir where to write.
#' @param base_url PGS Catalog REST base URL.
#' @param majority_threshold an ancestry is called "majority" at this fraction
#'   of the sample size, otherwise "Multi-ancestry".
#' @param request_delay seconds between API calls.
#' @param bucket_dir `gs://` folder to push the two CSVs to; NULL = skip.
#' @return Invisibly, list(score_metadata, performance).
#' @export
fetch_pgs_catalog_metadata <- function(
    pgs_ids = NULL, pgs_metadata_path = NULL, trait = NULL,
    out_dir = "pgs_comparison_results",
    base_url = "https://www.pgscatalog.org/rest",
    majority_threshold = 0.5, request_delay = 0.34,
    bucket_dir = default_bucket_dir()) {

  if (is.null(pgs_ids)) {
    if (is.null(pgs_metadata_path)) stop("Give `pgs_ids` or `pgs_metadata_path`.")
    meta_in <- fread(pgs_metadata_path)
    stopifnot("pgs_id" %in% names(meta_in))
    pgs_ids <- meta_in$pgs_id
  }
  pgs_ids <- unique(as.character(pgs_ids))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("Fetching PGS Catalog metadata for %d unique PGS IDs\n", length(pgs_ids)))
  file_prefix <- if (is.null(trait)) "" else paste0(trait, "_")

  ## ---- score-level metadata ----
  score_meta_list <- vector("list", length(pgs_ids))
  for (i in seq_along(pgs_ids)) {
    pgs_id <- pgs_ids[i]
    cat(sprintf("[%d/%d] %s\n", i, length(pgs_ids), pgs_id))
    sc <- get_score_metadata(pgs_id, base_url)
    Sys.sleep(request_delay)
    if (is.null(sc)) { score_meta_list[[i]] <- data.table(pgs_id = pgs_id); next }

    gwas <- summarize_ancestry(samples_list_to_dt(sc$samples_variants), majority_threshold)
    dev  <- summarize_ancestry(samples_list_to_dt(sc$samples_training), majority_threshold)
    score_meta_list[[i]] <- data.table(
      pgs_id = pgs_id,
      pgs_name = sc$name %||% NA_character_,
      trait_reported = sc$trait_reported %||% NA_character_,
      method_name = sc$method_name %||% NA_character_,
      variants_number = as.integer(sc$variants_number %||% NA),
      gwas_ancestry = gwas$ancestry, gwas_sample_size = gwas$n_total,
      dev_ancestry = dev$ancestry, dev_sample_size = dev$n_total
    )
  }
  score_meta <- rbindlist(score_meta_list, fill = TRUE)
  meta_file <- file.path(out_dir, paste0(file_prefix, "pgs_catalog_score_metadata.csv"))
  fwrite(score_meta, meta_file)
  cat(sprintf("Wrote score metadata for %d PGS to %s\n", nrow(score_meta), meta_file))

  ## ---- reported performance ----
  perf_list <- vector("list", length(pgs_ids))
  for (i in seq_along(pgs_ids)) {
    pgs_id <- pgs_ids[i]
    cat(sprintf("[perf %d/%d] %s\n", i, length(pgs_ids), pgs_id))
    entries <- get_all_performance_metrics(pgs_id, base_url, request_delay)
    if (length(entries) == 0) next
    perf_list[[i]] <- rbindlist(
      lapply(entries, extract_performance_rows, pgs_id = pgs_id, threshold = majority_threshold),
      fill = TRUE)
  }
  perf <- rbindlist(perf_list, fill = TRUE)
  perf_file <- file.path(out_dir, paste0(file_prefix, "pgs_catalog_performance_metrics.csv"))
  fwrite(perf, perf_file)
  cat(sprintf("Wrote %d catalog-reported performance entries (across %d PGS) to %s\n",
              nrow(perf), length(unique(perf$pgs_id)), perf_file))

  if (!is.null(bucket_dir)) {
    system2("gsutil", c("-m", "cp", shQuote(meta_file), shQuote(perf_file),
                        shQuote(paste0(bucket_dir, "/"))))
  }
  cat("Done.\n")
  invisible(list(score_metadata = score_meta, performance = perf))
}

# Parse sample-group objects by hand (simplifyVector = FALSE) so differently
# shaped nested fields don't break rbindlist.
samples_list_to_dt <- function(samples) {
  if (is.null(samples) || length(samples) == 0) return(data.table())
  rbindlist(lapply(samples, function(s) {
    data.table(sample_number  = suppressWarnings(as.numeric(s$sample_number %||% NA)),
               ancestry_broad = as.character(s$ancestry_broad %||% NA))
  }), fill = TRUE)
}

# Total sample size + majority ancestry (or "Multi-ancestry").
summarize_ancestry <- function(samples_dt, threshold) {
  if (is.null(samples_dt) || nrow(samples_dt) == 0) {
    return(list(ancestry = NA_character_, n_total = NA_integer_))
  }
  n_total <- sum(samples_dt$sample_number, na.rm = TRUE)

  # "European, East Asian" style rows: split and divide the count evenly; drop
  # blank pieces so a trailing comma doesn't create a phantom "" category.
  by_anc <- samples_dt[, {
    cats <- trimws(unlist(strsplit(ancestry_broad, ",")))
    cats <- cats[!is.na(cats) & cats != ""]
    if (length(cats) == 0) cats <- NA_character_
    list(ancestry_broad = cats, sample_number = sample_number / length(cats))
  }, by = seq_len(nrow(samples_dt))]

  totals <- by_anc[, list(n = sum(sample_number, na.rm = TRUE)), by = ancestry_broad]
  totals <- totals[order(-n)]
  if (nrow(totals) == 0 || n_total == 0 || is.na(n_total)) {
    return(list(ancestry = NA_character_, n_total = n_total))
  }
  ancestry <- if (totals$n[1] / n_total >= threshold) totals$ancestry_broad[1] else "Multi-ancestry"
  list(ancestry = ancestry, n_total = round(n_total))
}

get_score_metadata <- function(pgs_id, base_url) {
  resp <- httr::GET(sprintf("%s/score/%s", base_url, pgs_id))
  if (httr::status_code(resp) != 200) {
    warning(sprintf("Score lookup failed for %s (HTTP %d)", pgs_id, httr::status_code(resp)))
    return(NULL)
  }
  jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
}

get_all_performance_metrics <- function(pgs_id, base_url, request_delay) {
  url <- sprintf("%s/performance/search?pgs_id=%s", base_url, pgs_id)
  all_results <- list()
  while (!is.null(url)) {
    resp <- httr::GET(url)
    if (httr::status_code(resp) != 200) {
      warning(sprintf("Performance lookup failed for %s (HTTP %d)", pgs_id, httr::status_code(resp)))
      break
    }
    parsed <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
    if (!is.null(parsed$results) && length(parsed$results) > 0) {
      all_results <- c(all_results, parsed$results)
    }
    url <- parsed$`next`
    Sys.sleep(request_delay)
  }
  all_results
}

metrics_list_to_dt <- function(metrics, category) {
  if (is.null(metrics) || length(metrics) == 0) return(data.table())
  rbindlist(lapply(metrics, function(m) {
    data.table(metric_category = category,
               name_short = m$name_short %||% NA_character_,
               name_long  = m$name_long %||% NA_character_,
               estimate   = suppressWarnings(as.numeric(m$estimate %||% NA)),
               ci_lower   = suppressWarnings(as.numeric(m$ci_lower %||% NA)),
               ci_upper   = suppressWarnings(as.numeric(m$ci_upper %||% NA)))
  }), fill = TRUE)
}

extract_performance_rows <- function(pgs_id, entry, threshold) {
  eval_summary <- summarize_ancestry(samples_list_to_dt(entry$sampleset$samples), threshold)
  pm <- entry$performance_metrics
  metrics_dt <- rbindlist(list(
    metrics_list_to_dt(pm$effect_sizes, "effect_size"),
    metrics_list_to_dt(pm$class_acc,    "classification"),
    metrics_list_to_dt(pm$othermetrics, "other")
  ), fill = TRUE)
  if (nrow(metrics_dt) == 0) {
    metrics_dt <- data.table(metric_category = NA_character_, name_short = NA_character_,
                             name_long = NA_character_, estimate = NA_real_,
                             ci_lower = NA_real_, ci_upper = NA_real_)
  }
  covariates <- if (is.null(entry$covariates)) NA_character_ else paste(unlist(entry$covariates), collapse = ", ")
  base_row <- list(pgs_id = pgs_id, ppm_id = entry$id %||% NA_character_,
                   phenotyping_reported = entry$phenotyping_reported %||% NA_character_,
                   publication = entry$publication$id %||% NA_character_,
                   eval_ancestry = eval_summary$ancestry, eval_sample_size = eval_summary$n_total,
                   covariates = covariates)
  for (nm in names(base_row)) metrics_dt[, (nm) := base_row[[nm]]]
  setcolorder(metrics_dt, c(names(base_row), setdiff(names(metrics_dt), names(base_row))))
  metrics_dt
}
