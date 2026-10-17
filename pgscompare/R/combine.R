#' Combine per-PGS score files and evaluate them with pgsmetrics
#'
#' Reads the Terra data table of per-PGS score files, downloads and merges them
#' into one wide table, merges phenotypes, optionally subsets individuals, and
#' runs `pgsmetrics::pgsmetrics()` (plus a per-SD odds ratio for binary traits)
#' in chunks with per-chunk `.rds` checkpoints.
#'
#' Outputs in `out_dir` (all prefixed `<cohort_name>[_<subset_tag>]`):
#' `_<trait>_pgsmetrics.csv`, `_<trait>_pgsmetrics_cov.csv`,
#' `_<trait>_sample_size.csv`, `_<trait>_or_predictions.csv` (binary), plus
#' `<cohort_name>_combined_pgs_wide.csv` and `<cohort_name>_pgs_metadata.csv`.
#'
#' @param cohort_name short cohort name; prefixes every output file.
#' @param pheno_path phenotype CSV/TSV (`gs://` or local).
#' @param covariates covariate column names in the phenotype file.
#' @param binary_traits,continuous_traits named character vectors:
#'   output name -> phenotype column, e.g. `c(t2d = "t2d_status")`.
#' @param prevalence named numeric vector (names match `binary_traits`) of
#'   population prevalence, used for liability-scale R2.
#' @param pheno_id_col ID column in the phenotype file.
#' @param table_name Terra data table listing the score files.
#' @param pgs_id_col,file_path_col,ancestry_adjusted_col columns of that table.
#' @param sscore_id_col ID column inside each score file.
#' @param score_col_suffix suffix of the score column in RAW files (`"_AVG"`).
#' @param sampleset_col,keep_samplesets optional filter on a sampleset column.
#' @param exclude_file_ids values of `pgs_individual_file_id` to drop.
#' @param exclude_pgs_ids PGS ids to skip entirely (both adjusted and raw).
#' @param reuse_combined_scores reuse cached combined score/metadata CSVs.
#' @param subset_ids_path file of IIDs to restrict to (one per line, or CSV/TSV
#'   whose first column is IID). The file's basename becomes the `subset_tag`.
#' @param chunk_size max PGS per pgsmetrics call.
#' @param skip_existing reuse per-chunk `.rds` results if present.
#' @param n_boot,n_cores bootstrap replicates / cores for pgsmetrics.
#' @param n_sample_individuals,sample_seed individuals sampled for the OR
#'   individual-prediction output.
#' @param local_dir scratch folder for downloaded files.
#' @param out_dir results folder.
#' @param bucket_dir `gs://` folder to push this cohort's CSVs to; NULL = skip.
#'   Defaults to `<WORKSPACE_BUCKET>/pgs_comparison_results` on Terra.
#' @param pgs_table optional data.frame to use instead of reading `table_name`
#'   with `avtable()` (handy for testing outside Terra).
#' @return Invisibly, a list with `files` (paths written) and `results`.
#' @export
run_combine_and_compare <- function(
    cohort_name, pheno_path, covariates,
    binary_traits = NULL, continuous_traits = NULL, prevalence = NULL,
    pheno_id_col = "IID",
    table_name = "pgs_individual_file",
    pgs_id_col = "pgs_id", file_path_col = "sscore_path",
    sscore_id_col = "IID", score_col_suffix = "_AVG",
    ancestry_adjusted_col = "ancestry_adjusted",
    sampleset_col = "sampleset", keep_samplesets = NULL,
    exclude_file_ids = NULL, exclude_pgs_ids = NULL,
    reuse_combined_scores = TRUE, subset_ids_path = NULL,
    chunk_size = 50, skip_existing = TRUE, n_boot = 1000, n_cores = 4,
    n_sample_individuals = 6, sample_seed = 42,
    local_dir = "pgs_sscore_files", out_dir = "pgs_comparison_results",
    bucket_dir = default_bucket_dir(), pgs_table = NULL) {

  if (!requireNamespace("pgsmetrics", quietly = TRUE)) {
    stop("Package 'pgsmetrics' is required. Install it in a clean session with:\n",
         "  remotes::install_github('krillinor/pgsmetrics', upgrade = 'never')")
  }
  if (length(binary_traits) + length(continuous_traits) == 0) {
    stop("Provide at least one of `binary_traits` / `continuous_traits`.")
  }
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  subset_tag <- subset_tag_from_path(subset_ids_path)
  prefix     <- name_prefix(cohort_name, subset_tag)
  out_suffix <- if (is.null(subset_tag)) "" else paste0("_", subset_tag)

  ## ---- data table ----
  if (is.null(pgs_table)) {
    pgs_table <- as.data.table(resolve_avtable()(table_name))
  } else {
    pgs_table <- as.data.table(copy(pgs_table))
  }
  pgs_table <- prepare_pgs_table(pgs_table, table_name, pgs_id_col, file_path_col,
                                 ancestry_adjusted_col, sampleset_col, keep_samplesets,
                                 exclude_file_ids, exclude_pgs_ids)

  ## ---- combine score files (or reuse cache) ----
  combined_path <- file.path(out_dir, sprintf("%s_combined_pgs_wide.csv", cohort_name))
  meta_path     <- file.path(out_dir, sprintf("%s_pgs_metadata.csv", cohort_name))

  if (reuse_combined_scores && file.exists(combined_path) && file.exists(meta_path)) {
    cat(sprintf("Found existing combined scores at %s - loading instead of re-combining %d score files. Use reuse_combined_scores = FALSE to force.\n",
                combined_path, nrow(pgs_table)))
    pgs_wide <- fread(combined_path)
    pgs_meta <- fread(meta_path)
  } else {
    comb <- combine_score_files(pgs_table, pgs_id_col, file_path_col, sscore_id_col,
                                score_col_suffix, local_dir)
    pgs_wide <- comb$pgs_wide
    pgs_meta <- comb$pgs_meta
    fwrite(pgs_wide, combined_path)
    fwrite(pgs_meta, meta_path)
  }
  cat(sprintf("%d PGS (%d adjusted, %d raw) across %d individuals\n",
              nrow(pgs_meta), sum(pgs_meta$ancestry_adjusted),
              sum(!pgs_meta$ancestry_adjusted), nrow(pgs_wide)))

  ## ---- phenotypes ----
  pheno <- fread(localize_file(pheno_path, local_dir))
  if (!(pheno_id_col %in% names(pheno))) {
    stop(sprintf("ID column '%s' not found in phenotype file. Columns present: %s",
                 pheno_id_col, paste(names(pheno), collapse = ", ")))
  }
  setnames(pheno, pheno_id_col, "IID")
  analysis_df <- merge(pgs_wide, pheno, by = "IID")
  cat(sprintf("Merged analysis dataset: %d individuals with both scores and phenotypes\n",
              nrow(analysis_df)))
  missing_covars <- setdiff(covariates, names(analysis_df))
  if (length(missing_covars) > 0) {
    stop("Missing covariate column(s) in phenotype file: ", paste(missing_covars, collapse = ", "))
  }

  ## ---- optional subset ----
  if (!is.null(subset_ids_path)) {
    local_subset <- localize_file(subset_ids_path, local_dir)
    first_line <- readLines(local_subset, n = 1)
    subset_ids <- if (grepl(",|\t", first_line)) fread(local_subset)[[1]] else readLines(local_subset)
    subset_ids <- as.character(subset_ids)
    n_before   <- nrow(analysis_df)
    n_not_found <- length(subset_ids) - sum(analysis_df$IID %in% subset_ids)
    if (n_not_found > 0) {
      cat(sprintf("Note: %d of %d subset IDs were not found in the scores/phenotype data\n",
                  n_not_found, length(subset_ids)))
    }
    analysis_df <- analysis_df[analysis_df$IID %in% subset_ids]
    cat(sprintf("Subset to %d/%d individuals via subset_ids_path (%s)\n",
                nrow(analysis_df), n_before, subset_ids_path))
  }

  ## ---- make sure everything pgsmetrics needs is numeric ----
  analysis_df <- check_analysis_types(
    analysis_df,
    pgs_cols   = setdiff(names(pgs_wide), "IID"),
    covariates = covariates,
    dep_cols   = c(binary_traits, continuous_traits),
    binary_dep_cols = binary_traits)

  ## ---- run pgsmetrics ----
  ctx <- list(analysis_df = analysis_df, covariates = covariates, n_boot = n_boot,
              n_cores = n_cores, skip_existing = skip_existing, out_dir = out_dir,
              out_suffix = out_suffix, chunk_size = chunk_size,
              pgs_groups = {
                g <- list(adjusted = pgs_meta$unique_name[pgs_meta$ancestry_adjusted],
                          raw      = pgs_meta$unique_name[!pgs_meta$ancestry_adjusted])
                g[lengths(g) > 0]
              })

  results <- list()
  or_predictions <- list()
  set.seed(sample_seed)

  for (trait_name in names(binary_traits)) {
    dep_col <- binary_traits[[trait_name]]
    K <- if (trait_name %in% names(prevalence)) prevalence[[trait_name]] else NULL

    base_ok  <- stats::complete.cases(analysis_df[, c(dep_col, covariates), with = FALSE])
    base_df  <- data.table(IID = analysis_df$IID[base_ok], outcome = analysis_df[[dep_col]][base_ok])
    n_per    <- ceiling(n_sample_individuals / 2)
    cases    <- base_df$IID[base_df$outcome == 1]
    controls <- base_df$IID[base_df$outcome == 0]
    sample_ids <- c(sample_safely(cases, n_per), sample_safely(controls, n_per))
    cat(sprintf("Sampled %d individuals for '%s' (%d cases, %d controls) for OR prediction profiles\n",
                length(sample_ids), trait_name, sum(sample_ids %in% cases), sum(sample_ids %in% controls)))

    cat(sprintf("Running pgsmetrics for binary trait '%s' (%s)\n", trait_name, dep_col))
    out <- run_pgsmetrics_for_trait(ctx, trait_name, dep_col, is_binary = TRUE, K = K, sample_ids = sample_ids)
    results[[trait_name]] <- out$metrics
    or_predictions[[trait_name]] <- out$predictions$all
  }
  for (trait_name in names(continuous_traits)) {
    dep_col <- continuous_traits[[trait_name]]
    cat(sprintf("Running pgsmetrics for continuous trait '%s' (%s)\n", trait_name, dep_col))
    out <- run_pgsmetrics_for_trait(ctx, trait_name, dep_col, is_binary = FALSE)
    results[[trait_name]] <- out$metrics
  }

  ## ---- write outputs ----
  files <- c(combined_path, meta_path)
  all_traits <- c(binary_traits, continuous_traits)
  for (trait_name in names(all_traits)) {
    w <- widen_metrics(results[[trait_name]][["all"]])
    f <- file.path(out_dir, sprintf("%s_%s_pgsmetrics.csv", prefix, trait_name))
    fwrite(w$main, f); files <- c(files, f)
    if (!is.null(w$cov)) {
      f <- file.path(out_dir, sprintf("%s_%s_pgsmetrics_cov.csv", prefix, trait_name))
      fwrite(w$cov, f); files <- c(files, f)
    }

    dep_col <- all_traits[[trait_name]]
    usable  <- stats::complete.cases(analysis_df[, c(dep_col, covariates), with = FALSE])
    ss <- data.table(cohort = cohort_name,
                     subset_tag = if (is.null(subset_tag)) NA_character_ else subset_tag,
                     trait = trait_name, n_individuals = nrow(analysis_df), n_usable = sum(usable))
    if (trait_name %in% names(binary_traits)) {
      ss[, n_cases    := sum(analysis_df[[dep_col]][usable] == 1)]
      ss[, n_controls := sum(analysis_df[[dep_col]][usable] == 0)]
    }
    f <- file.path(out_dir, sprintf("%s_%s_sample_size.csv", prefix, trait_name))
    fwrite(ss, f); files <- c(files, f)
    cat(sprintf("Sample size for '%s' (%s): %d individuals, %d usable%s\n",
                trait_name, prefix, ss$n_individuals, ss$n_usable,
                if ("n_cases" %in% names(ss)) sprintf(" (%d cases / %d controls)", ss$n_cases, ss$n_controls) else ""))
  }
  for (trait_name in names(binary_traits)) {
    preds <- or_predictions[[trait_name]]
    if (is.null(preds) || nrow(preds) == 0) next
    f <- file.path(out_dir, sprintf("%s_%s_or_predictions.csv", prefix, trait_name))
    fwrite(preds, f); files <- c(files, f)
    cat(sprintf("Wrote %d individual OR-prediction rows for '%s' to %s\n", nrow(preds), trait_name, f))
  }

  cat(sprintf("Done. Results written to %s/\n", out_dir))
  sync_to_bucket(bucket_dir, out_dir, csv_glob = paste0(cohort_name, "*.csv"), plots = FALSE)
  invisible(list(files = files, results = results))
}

sample_safely <- function(x, n) if (length(x) == 0) x else x[sample.int(length(x), min(n, length(x)))]

## ---------------------------------------------------------------------------
## Data-table preparation and score-file combination
## ---------------------------------------------------------------------------

prepare_pgs_table <- function(pgs_table, table_name, pgs_id_col, file_path_col,
                              ancestry_adjusted_col, sampleset_col, keep_samplesets,
                              exclude_file_ids, exclude_pgs_ids) {
  for (col in c(pgs_id_col, file_path_col, ancestry_adjusted_col)) {
    if (!(col %in% names(pgs_table))) {
      stop(sprintf("Column '%s' not found in data table. Columns present: %s",
                   col, paste(names(pgs_table), collapse = ", ")))
    }
  }
  raw_flag <- pgs_table[[ancestry_adjusted_col]]
  adj <- if (is.logical(raw_flag)) raw_flag else
    tolower(as.character(raw_flag)) %in% c("true", "yes", "1", "adjusted")
  pgs_table[, ancestry_adjusted := adj]
  cat(sprintf("Loaded %d rows from data table '%s'\n", nrow(pgs_table), table_name))

  if (!is.null(keep_samplesets) && !is.null(sampleset_col) && sampleset_col %in% names(pgs_table)) {
    n_before <- nrow(pgs_table)
    keep <- pgs_table[[sampleset_col]] %in% keep_samplesets
    if (any(!keep)) {
      cat(sprintf("Dropping %d row(s) outside kept sampleset(s) %s:\n",
                  sum(!keep), paste(keep_samplesets, collapse = ", ")))
      print(data.table(pgs_id = pgs_table[[pgs_id_col]][!keep], sampleset = pgs_table[[sampleset_col]][!keep]))
    }
    pgs_table <- pgs_table[which(keep)]
    cat(sprintf("Kept %d/%d rows after sampleset filter\n", nrow(pgs_table), n_before))
  } else if (!is.null(keep_samplesets)) {
    warning("keep_samplesets is set (", paste(keep_samplesets, collapse = ", "),
            ") but the sampleset filter did NOT run - sampleset_col is ",
            if (is.null(sampleset_col)) "NULL" else sprintf("'%s', not found in the table", sampleset_col),
            ". Available columns: ", paste(names(pgs_table), collapse = ", "), ".")
  }

  if (!is.null(exclude_file_ids) && "pgs_individual_file_id" %in% names(pgs_table)) {
    n_before <- nrow(pgs_table)
    pgs_table <- pgs_table[which(!(pgs_table[["pgs_individual_file_id"]] %in% exclude_file_ids))]
    cat(sprintf("Excluded %d row(s) by pgs_individual_file_id; %d/%d rows remain\n",
                n_before - nrow(pgs_table), nrow(pgs_table), n_before))
  }
  if (!is.null(exclude_pgs_ids)) {
    n_before <- nrow(pgs_table)
    pgs_table <- pgs_table[which(!(pgs_table[[pgs_id_col]] %in% exclude_pgs_ids))]
    cat(sprintf("Excluded PGS %s entirely; %d/%d rows remain\n",
                paste(exclude_pgs_ids, collapse = ", "), nrow(pgs_table), n_before))
  }
  pgs_table
}

# RAW files: "#IID", ..., "<PGS_ID>_AVG". ADJUSTED files: "IID" + bare "<PGS_ID>".
read_one_score_file <- function(path, pgs_name, id_col, score_col_suffix, is_adjusted) {
  dt <- fread(path)
  names(dt) <- sub("^#", "", names(dt))
  score_col <- if (is_adjusted) pgs_name else paste0(pgs_name, score_col_suffix)
  for (col in c(id_col, score_col)) {
    if (!(col %in% names(dt))) {
      stop(sprintf("Column '%s' not found in %s. Columns present: %s",
                   col, path, paste(names(dt), collapse = ", ")))
    }
  }
  out <- data.table(IID = dt[[id_col]], score = dt[[score_col]])
  setnames(out, "score", pgs_name)
  out
}

combine_score_files <- function(pgs_table, pgs_id_col, file_path_col, sscore_id_col,
                                score_col_suffix, local_dir) {
  # Unique column per score: "<pgs_id>_adj" / "<pgs_id>_raw".
  pgs_table[, unique_name := paste0(pgs_table[[pgs_id_col]], ifelse(ancestry_adjusted, "_adj", "_raw"))]

  if (anyDuplicated(pgs_table$unique_name) > 0) {
    dup_names <- unique(pgs_table$unique_name[
      duplicated(pgs_table$unique_name) | duplicated(pgs_table$unique_name, fromLast = TRUE)])
    cat("Duplicate (pgs_id, ancestry_adjusted) rows found - full colliding rows:\n")
    print(pgs_table[pgs_table$unique_name %in% dup_names])
    stop("Duplicate (pgs_id, ancestry_adjusted) combinations in the data table (printed above). ",
         "If driven by a 'sampleset' column set `keep_samplesets`; if it is a genuine duplicate ",
         "entry add the unwanted row's pgs_individual_file_id to `exclude_file_ids`, or skip the ",
         "PGS entirely via `exclude_pgs_ids`, then rerun.")
  }

  pgs_meta <- data.table(pgs_id = as.character(pgs_table[[pgs_id_col]]),
                         ancestry_adjusted = pgs_table$ancestry_adjusted,
                         unique_name = pgs_table$unique_name)

  score_list <- vector("list", nrow(pgs_table))
  for (i in seq_len(nrow(pgs_table))) {
    pgs_id      <- as.character(pgs_table[[pgs_id_col]][i])
    unique_name <- pgs_table$unique_name[i]
    local_path  <- localize_file(as.character(pgs_table[[file_path_col]][i]), local_dir)
    cat(sprintf("[%d/%d] Reading %s (%s, %s)\n", i, nrow(pgs_table), unique_name,
                ifelse(pgs_table$ancestry_adjusted[i], "ancestry-adjusted", "raw"), local_path))
    score_list[[i]] <- read_one_score_file(local_path, pgs_id, sscore_id_col, score_col_suffix,
                                           pgs_table$ancestry_adjusted[i])
    setnames(score_list[[i]], pgs_id, unique_name)
  }
  list(pgs_wide = Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), score_list),
       pgs_meta = pgs_meta)
}

## ---------------------------------------------------------------------------
## Metrics
## ---------------------------------------------------------------------------

# pgsmetrics() has no odds ratio, so fit one logistic model per PGS
# (dep ~ scale(pgs) + covars) and report OR per 1 SD with a Wald CI, in the same
# long shape pgsmetrics uses (metric = "or", type = "full"). If `sample_ids` is
# given, also return those individuals' predicted probabilities from the same fits.
compute_or_metrics <- function(df, chunk_pgs, dep_col, covars, sample_ids = NULL) {
  pred_rows <- list()
  metric_rows <- lapply(chunk_pgs, function(pgs_name) {
    d <- copy(df)
    d[, (pgs_name) := as.numeric(scale(d[[pgs_name]]))]
    fmla <- stats::as.formula(paste(dep_col, "~", paste(c(pgs_name, covars), collapse = " + ")))
    fit <- tryCatch(stats::glm(fmla, data = d, family = stats::binomial()), error = function(e) NULL)

    if (!is.null(fit) && !is.null(sample_ids)) {
      present <- d[d$IID %in% sample_ids]
      if (nrow(present) > 0) {
        pred_rows[[pgs_name]] <<- data.table(
          IID = present$IID, pgs = pgs_name,
          predicted_prob = stats::predict(fit, newdata = present, type = "response"),
          actual_outcome = present[[dep_col]])
      }
    }
    if (is.null(fit) || !(pgs_name %in% names(stats::coef(fit)))) {
      return(data.table(pgs = pgs_name, metric = "or", type = "full",
                        observed = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
    }
    est <- stats::coef(fit)[[pgs_name]]
    se  <- summary(fit)$coefficients[pgs_name, "Std. Error"]
    data.table(pgs = pgs_name, metric = "or", type = "full", observed = exp(est),
               ci_lower = exp(est - 1.96 * se), ci_upper = exp(est + 1.96 * se))
  })
  list(metrics = rbindlist(metric_rows, fill = TRUE),
       predictions = if (length(pred_rows) > 0) rbindlist(pred_rows, fill = TRUE) else NULL)
}

run_one_chunk <- function(ctx, trait_name, group_name, chunk_idx, n_chunks, chunk_pgs,
                          dep_col, is_binary, K, sample_ids = NULL) {
  out_rds <- file.path(ctx$out_dir, sprintf("pgsmetrics_%s_%s%s_chunk%02d.rds",
                                            trait_name, group_name, ctx$out_suffix, chunk_idx))
  df <- ctx$analysis_df
  needed_cols <- c(dep_col, ctx$covariates, chunk_pgs)

  if (ctx$skip_existing && file.exists(out_rds)) {
    cached <- readRDS(out_rds)
    if (is_binary && !("or" %in% cached$metrics$metric)) {
      cat(sprintf("  %s chunk %d/%d: cached result has no OR yet - backfilling OR only\n",
                  group_name, chunk_idx, n_chunks))
      chunk_df <- df[stats::complete.cases(df[, needed_cols, with = FALSE])]
      or_result <- compute_or_metrics(chunk_df, chunk_pgs, dep_col, ctx$covariates, sample_ids)
      cached$metrics <- rbindlist(list(cached$metrics, or_result$metrics), fill = TRUE)
      cached$or_predictions <- or_result$predictions
      saveRDS(cached, out_rds)
      rm(chunk_df); gc(full = TRUE, verbose = FALSE)
    } else {
      cat(sprintf("  %s chunk %d/%d: found existing result at %s, skipping\n",
                  group_name, chunk_idx, n_chunks, out_rds))
    }
    return(list(metrics = cached$metrics, predictions = cached$or_predictions))
  }

  complete <- stats::complete.cases(df[, needed_cols, with = FALSE])
  if (any(!complete)) {
    cat(sprintf("  %s chunk %d/%d: dropping %d/%d rows with missing values\n",
                group_name, chunk_idx, n_chunks, sum(!complete), nrow(df)))
  }
  chunk_df <- df[which(complete)]
  cat(sprintf("  %s chunk %d/%d (%d PGS, n=%d) starting at %s\n", group_name, chunk_idx, n_chunks,
              length(chunk_pgs), nrow(chunk_df), format(Sys.time())))

  attach_pgsmetrics_deps()
  call_pgsmetrics <- function(pgs, n_boot, n_cores) {
    pgsmetrics::pgsmetrics(chunk_df, pgs = pgs, dep = dep_col, covars = ctx$covariates,
                           K = if (is_binary) K else NULL, n_boot = n_boot, n_cores = n_cores)
  }
  # mclapply workers run silently, so a failure inside them only shows up as
  # "all scheduled cores encountered errors in user code". Catch that warning and
  # re-run a tiny serial job (n_cores = 1) so the REAL error message is shown.
  res <- withCallingHandlers(
    call_pgsmetrics(chunk_pgs, ctx$n_boot, ctx$n_cores),
    warning = function(w) {
      if (grepl("encountered errors in user code", conditionMessage(w), fixed = TRUE)) {
        cat("\n*** pgsmetrics bootstrap workers failed; re-running 2 PGS serially to find the real error ***\n")
        err <- tryCatch({ call_pgsmetrics(chunk_pgs[seq_len(min(2, length(chunk_pgs)))], 2, 1); NULL },
                        error = function(e) e)
        if (!is.null(err)) {
          stop("pgsmetrics failed for ", trait_name, " (", group_name, " chunk ", chunk_idx,
               "). Underlying error: ", conditionMessage(err), call. = FALSE)
        }
        stop("pgsmetrics bootstrap workers failed in parallel but a small serial run succeeded - ",
             "most likely out of memory in the forked workers. Lower n_cores (try 1-2) and/or ",
             "chunk_size, then rerun (finished chunks are cached).", call. = FALSE)
      }
    })

  or_predictions <- NULL
  if (is_binary) {
    or_result <- compute_or_metrics(chunk_df, chunk_pgs, dep_col, ctx$covariates, sample_ids)
    res$metrics <- rbindlist(list(res$metrics, or_result$metrics), fill = TRUE)
    or_predictions <- or_result$predictions
  }
  res$or_predictions <- or_predictions

  saveRDS(res, out_rds)  # checkpoint so a later OOM crash doesn't lose this chunk
  cat(sprintf("  %s chunk %d/%d: done at %s, saved to %s\n",
              group_name, chunk_idx, n_chunks, format(Sys.time()), out_rds))
  metrics <- res$metrics
  rm(res, chunk_df); gc(full = TRUE, verbose = FALSE)
  list(metrics = metrics, predictions = or_predictions)
}

run_pgsmetrics_for_trait <- function(ctx, trait_name, dep_col, is_binary, K = NULL, sample_ids = NULL) {
  group_metrics <- list(); group_predictions <- list()
  for (group_name in names(ctx$pgs_groups)) {
    chunks <- chunk_vector(ctx$pgs_groups[[group_name]], ctx$chunk_size)
    cat(sprintf("Trait '%s', group '%s': %d PGS split into %d chunk(s) of up to %d\n",
                trait_name, group_name, length(ctx$pgs_groups[[group_name]]), length(chunks), ctx$chunk_size))
    chunk_results <- lapply(seq_along(chunks), function(j) {
      run_one_chunk(ctx, trait_name, group_name, j, length(chunks), chunks[[j]],
                    dep_col, is_binary, K, sample_ids)
    })
    group_metrics[[group_name]] <- rbindlist(lapply(chunk_results, `[[`, "metrics"), use.names = TRUE, fill = TRUE)
    group_predictions[[group_name]] <- rbindlist(lapply(chunk_results, `[[`, "predictions"), use.names = TRUE, fill = TRUE)
  }
  group_metrics$all     <- rbindlist(group_metrics, use.names = TRUE, fill = TRUE)
  group_predictions$all <- rbindlist(group_predictions, use.names = TRUE, fill = TRUE)
  list(metrics = group_metrics, predictions = group_predictions)
}

## ---------------------------------------------------------------------------
## Long -> wide
## ---------------------------------------------------------------------------

relabel_value_suffix <- function(x) {
  for (prefix in c("observed", "ci_lower", "ci_upper")) {
    hit <- grepl(paste0("^", prefix, "_"), x)
    x[hit] <- paste0(sub(paste0("^", prefix, "_"), "", x[hit]), "_", prefix)
  }
  x
}

# main: one row per PGS, <metric>_<full|partial>_<observed|ci_lower|ci_upper>.
# cov: covariate-only model metrics, consolidated to a single row.
widen_metrics <- function(dt) {
  cov_dt  <- dt[dt$type == "cov"]
  main_dt <- dt[dt$type != "cov"]

  cov_wide <- NULL
  if (nrow(cov_dt) > 0) {
    cov_full <- dcast(cov_dt, pgs ~ metric, value.var = c("observed", "ci_lower", "ci_upper"))
    value_cols <- setdiff(names(cov_full), "pgs")
    distinct_rows <- unique(cov_full[, value_cols, with = FALSE])
    if (nrow(distinct_rows) > 1) {
      cat(sprintf("Note: covariate-only metrics differ across %d value combinations (different chunk complete-case n?) - keeping the first.\n",
                  nrow(distinct_rows)))
    }
    cov_wide <- distinct_rows[1]
    setnames(cov_wide, names(cov_wide), relabel_value_suffix(names(cov_wide)))
    setcolorder(cov_wide, sort(names(cov_wide)))
  }

  main_wide <- dcast(main_dt, pgs ~ metric + type, value.var = c("observed", "ci_lower", "ci_upper"))
  value_cols <- setdiff(names(main_wide), "pgs")
  setnames(main_wide, value_cols, relabel_value_suffix(value_cols))
  setcolorder(main_wide, c("pgs", sort(setdiff(names(main_wide), "pgs"))))
  list(main = main_wide, cov = cov_wide)
}

# The original working scripts attached dplyr and pgsmetrics; some pgsmetrics code
# paths may rely on that (e.g. in forked bootstrap workers), so attach them when
# installed instead of assuming a clean search path.
attach_pgsmetrics_deps <- function() {
  for (pkg in c("dplyr", "pgsmetrics")) {
    if (requireNamespace(pkg, quietly = TRUE) && !(paste0("package:", pkg) %in% search())) {
      suppressPackageStartupMessages(attachNamespace(pkg))
    }
  }
  invisible(NULL)
}

# pgsmetrics hands the data to compiled code that requires numeric columns; a
# single character column (e.g. sex coded "Female"/"Male") fails with
# "Not compatible with requested type: [type=character; target=double]" -
# inside the workers, so it surfaces as a vague parallel-bootstrap error.
# Check up front and say exactly which column is the problem.
check_analysis_types <- function(df, pgs_cols, covariates, dep_cols, binary_dep_cols = NULL) {
  df <- copy(df)
  is_num <- function(x) is.numeric(x) || is.logical(x)

  # score columns: coerce character -> numeric (e.g. "NA"/"nan" text), reporting NAs created
  for (col in pgs_cols) {
    if (!is_num(df[[col]])) {
      before <- sum(is.na(df[[col]]))
      set(df, j = col, value = suppressWarnings(as.numeric(as.character(df[[col]]))))
      cat(sprintf("Note: score column '%s' was %s; coerced to numeric (%d new NA)\n",
                  col, "non-numeric", sum(is.na(df[[col]])) - before))
    }
  }

  # covariates and outcomes: don't guess a coding - stop and say what to fix
  bad <- character(0)
  for (col in unique(c(covariates, dep_cols))) {
    x <- df[[col]]
    if (is.logical(x)) { set(df, j = col, value = as.integer(x)); next }
    if (!is.numeric(x)) {
      vals <- utils::head(sort(unique(as.character(x[!is.na(x)]))), 5)
      bad <- c(bad, sprintf("  - '%s' is %s (e.g. %s)", col, class(x)[1], paste(shQuote(vals), collapse = ", ")))
    }
  }
  if (length(bad) > 0) {
    stop("These covariate/outcome columns are not numeric, which pgsmetrics requires:\n",
         paste(bad, collapse = "\n"),
         "\nRecode them as numbers in the phenotype file (e.g. sex 0/1, binary outcomes 0/1) and rerun.",
         call. = FALSE)
  }

  for (nm in names(binary_dep_cols)) {
    col <- binary_dep_cols[[nm]]
    vals <- unique(df[[col]][!is.na(df[[col]])])
    if (!all(vals %in% c(0, 1))) {
      warning(sprintf("Binary outcome '%s' (column '%s') has values other than 0/1: %s. Cases/controls ",
                      nm, col, paste(utils::head(sort(vals), 5), collapse = ", ")),
              "are assumed to be 1/0 for sample sizes and OR predictions; recode (e.g. 1/2 -> 0/1) if needed.",
              call. = FALSE)
    }
  }
  df
}
