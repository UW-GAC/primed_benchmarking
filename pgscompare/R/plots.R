#' Default metric-vs-metric comparisons for [run_join_and_plot()]
#' @return list of `list(x, y, color)`.
#' @export
default_metric_comparisons <- function() {
  list(list(x = "roc_auc", y = "r2_nagelkerke", color = "dev_ancestry"))
}

#' Default mapping of your metrics to PGS Catalog `name_short` families
#'
#' Each entry: `own_metric`, `own_type` ("full"/"partial"), a
#' `catalog_pattern` regex matched against the catalog's `name_short`, and a
#' `label`.
#' @export
default_catalog_comparisons <- function() {
  list(
    list(own_metric = "roc_auc", own_type = "partial",
         catalog_pattern = "Incremental AUROC|∆AUC", label = "Incremental AUC"),
    list(own_metric = "roc_auc", own_type = "full",
         catalog_pattern = "^AUROC$|PGS AUROC \\(no covariates\\)|C-index",
         label = "Full-model AUROC / C-index"),
    list(own_metric = "r2_nagelkerke", own_type = "partial",
         catalog_pattern = "Incremental R2|Nagelkerke R2|Covariate-adjusted pseudo-R2",
         label = "Incremental / Nagelkerke R2"),
    list(own_metric = "r2_liability", own_type = "partial",
         catalog_pattern = "liability R2", label = "Liability R2"),
    list(own_metric = "or", own_type = "full",
         # plain / per-allele / per-SD ORs only; tertile/percentile ORs are on a different scale
         catalog_pattern = "^OR$|per risk allele|per SD|per standard deviation",
         label = "Odds Ratio (per SD)")
  )
}

#' Join your pgsmetrics results with PGS Catalog data and make comparison plots
#'
#' Pulls results from the bucket (if `bucket_dir` is set), joins each trait's
#' `<prefix>_<trait>_pgsmetrics.csv` with the cohort's PGS metadata, match rate
#' and the PGS Catalog score metadata, then plots: metrics vs ancestry / sample
#' size / variants / match rate; metric vs metric; adjusted vs raw; yours vs
#' PGS Catalog; and odds-ratio summaries. Plots go to
#' `<plot_dir>/<cohort>[/<subset_tag>]/<trait>/` and are pushed to the bucket.
#'
#' @param cohort_name cohort name used in [run_combine_and_compare()].
#' @param subset_tag the run's subset tag (basename of its `subset_ids_path`
#'   without extension), or NULL for the full-cohort run.
#' @param results_dir local results folder.
#' @param traits traits to plot; NULL auto-detects.
#' @param metric_types named list overriding the full/partial variant per
#'   metric; others use `default_type`.
#' @param default_type `"partial"` or `"full"`.
#' @param max_axis_ratio use equal x/y scaling only when axis ranges are within
#'   this ratio.
#' @param metrics_to_plot named list: trait -> metrics; NULL/absent = all.
#' @param plot_dir root plots folder (default `<results_dir>/plots`).
#' @param match_rate_path force a single match-rate CSV for all traits.
#' @param metric_comparisons see [default_metric_comparisons()].
#' @param catalog_comparisons see [default_catalog_comparisons()].
#' @param bucket_dir `gs://` folder to pull inputs from / push outputs to.
#' @param sync_inputs pull from the bucket first?
#' @return Invisibly, the joined master table.
#' @export
run_join_and_plot <- function(
    cohort_name, subset_tag = NULL, results_dir = "pgs_comparison_results",
    traits = NULL, metric_types = list(or = "full"), default_type = "partial",
    max_axis_ratio = 1.5, metrics_to_plot = list(),
    plot_dir = file.path(results_dir, "plots"), match_rate_path = NULL,
    metric_comparisons = default_metric_comparisons(),
    catalog_comparisons = default_catalog_comparisons(),
    bucket_dir = default_bucket_dir(), sync_inputs = TRUE) {

  opts <- list(cohort_name = cohort_name, metric_types = metric_types,
               default_type = default_type, max_axis_ratio = max_axis_ratio)
  prefix <- name_prefix(cohort_name, subset_tag)
  plot_dir <- file.path(plot_dir, cohort_name)
  if (!is.null(subset_tag)) plot_dir <- file.path(plot_dir, subset_tag)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  if (sync_inputs) sync_from_bucket(bucket_dir, results_dir)

  ld <- load_master(cohort_name, subset_tag, results_dir, traits, match_rate_path)
  master <- ld$master; traits <- ld$traits
  metric_names_all <- unique(sub("_(full|partial)_observed$", "",
                                 grep("_(full|partial)_observed$", names(master), value = TRUE)))
  cat("Available metrics:", paste(sort(metric_names_all), collapse = ", "), "\n")

  fwrite(master, file.path(results_dir, sprintf("%s_pgs_comparison_master_table_wide.csv", prefix)))

  metrics_for <- function(trait_name) {
    m <- metrics_to_plot[[trait_name]]
    if (is.null(m)) sort(metric_names_all) else m
  }
  numeric_x_vars     <- c("gwas_sample_size", "dev_sample_size", "variants_number", "match_rate_percent")
  categorical_x_vars <- c("gwas_ancestry", "dev_ancestry")

  for (trait_name in traits) {
    tdir <- file.path(plot_dir, trait_name)
    dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
    trait_dt <- master[master$trait == trait_name]

    ## metric vs ancestry / sample size / variants / match rate
    for (metric_name in metrics_for(trait_name)) {
      if (!(metric_col(metric_name, "observed", opts) %in% names(trait_dt))) {
        cat(sprintf("  [%s/%s] column not found - skipping\n", trait_name, metric_name)); next
      }
      for (x_var in intersect(numeric_x_vars, names(trait_dt))) {
        p <- plot_metric_vs_numeric(trait_dt, x_var, metric_name, trait_name, opts)
        if (!is.null(p)) save_plot(p, file.path(tdir, sprintf("%s_vs_%s.png", metric_name, x_var)), 7, 5)
      }
      for (x_var in intersect(categorical_x_vars, names(trait_dt))) {
        p <- plot_metric_vs_categorical(trait_dt, x_var, metric_name, trait_name, opts)
        if (!is.null(p)) save_plot(p, file.path(tdir, sprintf("%s_vs_%s.png", metric_name, x_var)), 7, 5)
      }
    }

    ## metric vs metric
    for (cmp in metric_comparisons) {
      color_var <- if (!is.null(cmp$color)) cmp$color else "ancestry_adjusted"
      p <- plot_metric_vs_metric(master, trait_name, cmp$x, cmp$y, color_var, opts)
      if (!is.null(p)) {
        save_plot(p, file.path(tdir, sprintf("%s_vs_%s_by_%s.png", cmp$x, cmp$y, color_var)), 7, 5)
      }
    }

    ## adjusted vs raw
    for (metric_name in metrics_for(trait_name)) {
      type_name <- type_for_metric(metric_name, opts)
      plots <- plot_adjusted_vs_raw(master, trait_name, metric_name, opts)
      for (kind in names(plots)) {
        save_plot(plots[[kind]], file.path(tdir, sprintf("adj_vs_raw_%s_%s_%s.png", metric_name, type_name, kind)), 7, 6)
      }
    }

    ## yours vs PGS Catalog
    for (mapping in catalog_comparisons) {
      p <- plot_own_vs_catalog(master, ld$catalog_perf[[trait_name]], trait_name, mapping, opts)
      if (!is.null(p)) {
        save_plot(p, file.path(tdir, sprintf("vs_catalog_%s_%s.png", mapping$own_metric, mapping$own_type)), 9, 5)
      }
    }

    ## odds ratios
    p <- plot_or_pooled(master, trait_name, opts)
    if (!is.null(p)) save_plot(p, file.path(tdir, "or_pooled_summary.png"), 6, 5)

    p <- plot_or_forest(master, trait_name, opts)
    if (!is.null(p)) {
      obs_col <- metric_col("or", "observed", opts)
      n_pgs <- length(unique(master$pgs[master$trait == trait_name & !is.na(master[[obs_col]])]))
      save_plot(p, file.path(tdir, "or_forest.png"), 8, max(6, 0.12 * n_pgs), limitsize = FALSE)
    }
    p <- plot_or_individual_predictions(master, trait_name, results_dir, prefix, opts)
    if (!is.null(p)) save_plot(p, file.path(tdir, "or_individual_predictions.png"), 8, 5)
  }

  cat("Done generating plots.\n")
  sync_to_bucket(bucket_dir, results_dir, csv_glob = "*.csv", plots = TRUE)
  invisible(master)
}

## ---------------------------------------------------------------------------
## Loading / joining
## ---------------------------------------------------------------------------

load_master <- function(cohort_name, subset_tag, results_dir, traits, match_rate_path) {
  prefix <- name_prefix(cohort_name, subset_tag)
  if (length(traits) == 0) {
    traits <- detect_traits(results_dir, cohort_name, subset_tag)
  }
  if (length(traits) == 0) stop("No '", prefix, "_<trait>_pgsmetrics.csv' files found in ", results_dir)
  cat(sprintf("Loading pgsmetrics results for cohort '%s', traits: %s\n", cohort_name,
              paste(traits, collapse = ", ")))

  catalog_perf <- list()
  pieces <- list()

  for (trait_name in traits) {
    main_path <- file.path(results_dir, sprintf("%s_%s_pgsmetrics.csv", prefix, trait_name))
    if (!file.exists(main_path)) {
      warning(sprintf("Main metrics file not found for trait '%s': %s - skipping", trait_name, main_path))
      next
    }
    dt <- fread(main_path)
    dt[, trait := trait_name]

    meta_path <- resolve_optional(results_dir, cohort_name, "pgs_metadata", trait_name)
    if (!is.na(meta_path)) {
      meta <- fread(meta_path)
      stopifnot(all(c("pgs_id", "ancestry_adjusted", "unique_name") %in% names(meta)))
      dt <- merge(dt, meta, by.x = "pgs", by.y = "unique_name", all.x = TRUE)
      n_unmatched <- sum(is.na(dt$pgs_id))
      if (n_unmatched > 0) {
        cat(sprintf("  [%s] Warning: %d rows did not match a pgs_id in %s (check naming)\n",
                    trait_name, n_unmatched, meta_path))
      }
    } else {
      warning(sprintf("pgs_metadata not found for trait '%s' (looked for '%s_%s_pgs_metadata.csv' and '%s_pgs_metadata.csv')",
                      trait_name, cohort_name, trait_name, cohort_name))
    }

    # match rate: same value for the "_adj" and "_raw" version of a score
    mpath <- if (!is.null(match_rate_path)) match_rate_path else
      resolve_optional(results_dir, cohort_name, "match_rate", trait_name)
    if (!is.na(mpath) && file.exists(mpath)) {
      mr <- fread(mpath)
      setnames(mr, "pgs_model_id", "pgs_id", skip_absent = TRUE)
      setnames(mr, "percent", "match_rate_percent", skip_absent = TRUE)
      mr <- unique(mr[, list(pgs_id, match_rate_percent, total_count, count)])
      if (anyDuplicated(mr$pgs_id) > 0) {
        mr <- mr[, list(match_rate_percent = mean(match_rate_percent),
                        total_count = sum(total_count), count = sum(count)), by = pgs_id]
      }
      dt <- merge(dt, mr, by = "pgs_id", all.x = TRUE)
      cat(sprintf("  [%s] joined match-rate data for %d PGS from %s\n", trait_name, nrow(mr), mpath))
    } else {
      cat(sprintf("  [%s] no match-rate file found - match-rate plots will be skipped\n", trait_name))
    }

    cmeta_path <- resolve_trait_or_global(results_dir, sprintf("%s_pgs_catalog_score_metadata.csv", trait_name),
                                          "pgs_catalog_score_metadata.csv")
    if (!is.na(cmeta_path)) {
      cmeta <- fread(cmeta_path)
      stopifnot("pgs_id" %in% names(cmeta))
      for (col in names(cmeta)[vapply(cmeta, is.character, logical(1))]) {
        cmeta[which(cmeta[[col]] == ""), (col) := NA_character_]  # "" -> NA everywhere
      }
      dt <- merge(dt, cmeta, by = "pgs_id", all.x = TRUE)
      cat(sprintf("  [%s] joined PGS Catalog score metadata for %d PGS from %s\n", trait_name, nrow(cmeta), cmeta_path))
    } else {
      cat(sprintf("  [%s] no PGS Catalog score metadata found - ancestry/sample-size/variant plots will be skipped\n", trait_name))
    }

    cperf_path <- resolve_trait_or_global(results_dir, sprintf("%s_pgs_catalog_performance_metrics.csv", trait_name),
                                          "pgs_catalog_performance_metrics.csv")
    catalog_perf[[trait_name]] <- if (!is.na(cperf_path)) fread(cperf_path) else {
      cat(sprintf("  [%s] no PGS Catalog performance metrics found - yours-vs-catalog plots will be skipped\n", trait_name))
      data.table(pgs_id = character(0), name_short = character(0))
    }

    cov_path <- file.path(results_dir, sprintf("%s_%s_pgsmetrics_cov.csv", prefix, trait_name))
    if (file.exists(cov_path)) {
      cat(sprintf("Covariate-only model metrics for '%s' (%s):\n", trait_name, cov_path))
      print(fread(cov_path))
    }
    pieces[[trait_name]] <- dt
  }

  master <- rbindlist(pieces, fill = TRUE)
  cat(sprintf("Loaded %d rows across %d traits\n", nrow(master), length(pieces)))
  list(master = master, catalog_perf = catalog_perf, traits = names(pieces))
}

## ---------------------------------------------------------------------------
## Plot builders
## ---------------------------------------------------------------------------

plot_metric_vs_numeric <- function(d, x_var, metric_name, trait_name, opts) {
  obs_col <- metric_col(metric_name, "observed", opts)
  lo_col  <- metric_col(metric_name, "ci_lower", opts)
  hi_col  <- metric_col(metric_name, "ci_upper", opts)
  d <- d[which(!is.na(d[[x_var]]) & !is.na(d[[obs_col]]))]
  if (nrow(d) < 3) return(NULL)

  ggplot(d, aes(x = .data[[x_var]], y = .data[[obs_col]], color = ancestry_adjusted)) +
    geom_errorbar(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), width = 0, alpha = 0.4) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(title = sprintf("%s — %s (%s): %s vs %s", opts$cohort_name, trait_name,
                         type_for_metric(metric_name, opts), metric_name, x_var),
         x = x_var, y = sprintf("%s (with 95%% CI)", metric_name), color = "Ancestry-adjusted") +
    theme_minimal()
}

plot_metric_vs_categorical <- function(d, x_var, metric_name, trait_name, opts) {
  obs_col <- metric_col(metric_name, "observed", opts)
  d <- d[which(!is.na(d[[x_var]]) & !is.na(d[[obs_col]]))]
  if (nrow(d) < 3) return(NULL)

  ggplot(d, aes(x = .data[[x_var]], y = .data[[obs_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(aes(color = ancestry_adjusted), width = 0.15, size = 2, alpha = 0.8) +
    labs(title = sprintf("%s — %s (%s): %s by %s", opts$cohort_name, trait_name,
                         type_for_metric(metric_name, opts), metric_name, x_var),
         x = x_var, y = metric_name, color = "Ancestry-adjusted") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_metric_vs_metric <- function(master, trait_name, x_metric, y_metric, color_var, opts) {
  x_col <- metric_col(x_metric, "observed", opts)
  y_col <- metric_col(y_metric, "observed", opts)
  if (!(x_col %in% names(master)) || !(y_col %in% names(master))) return(NULL)
  d <- master[which(master$trait == trait_name & !is.na(master[[x_col]]) & !is.na(master[[y_col]]))]
  if (nrow(d) < 3) return(NULL)
  # the column can exist (filled in from another trait) yet be all-NA for this trait
  if (!(color_var %in% names(d)) || all(is.na(d[[color_var]]))) {
    cat(sprintf("  [%s] no values for color variable '%s' - falling back to ancestry_adjusted\n", trait_name, color_var))
    color_var <- "ancestry_adjusted"
  }
  corr <- stats::cor(d[[x_col]], d[[y_col]], use = "complete.obs")

  ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_var]])) +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = sprintf("%s — %s: %s vs %s", opts$cohort_name, trait_name, x_metric, y_metric),
         subtitle = sprintf("Pearson r = %.2f, colored by %s", corr, color_var),
         x = x_metric, y = y_metric, color = color_var) +
    theme_minimal()
}

plot_adjusted_vs_raw <- function(master, trait_name, metric_name, opts) {
  type_name <- type_for_metric(metric_name, opts)
  obs_col <- metric_col(metric_name, "observed", opts)
  lo_col  <- metric_col(metric_name, "ci_lower", opts)
  hi_col  <- metric_col(metric_name, "ci_upper", opts)
  if (!(obs_col %in% names(master)) || !("ancestry_adjusted" %in% names(master))) return(list())

  d <- master[which(master$trait == trait_name & !is.na(master$ancestry_adjusted) & !is.na(master[[obs_col]]))]
  if (nrow(d) == 0) return(list())
  d[, group := ifelse(ancestry_adjusted, "adjusted", "raw")]

  wide <- dcast(d, pgs_id ~ group, value.var = c(obs_col, lo_col, hi_col))
  setnames(wide,
           old = c(paste0(obs_col, "_adjusted"), paste0(obs_col, "_raw"),
                   paste0(lo_col, "_adjusted"), paste0(lo_col, "_raw"),
                   paste0(hi_col, "_adjusted"), paste0(hi_col, "_raw")),
           new = c("observed_adjusted", "observed_raw", "ci_lower_adjusted", "ci_lower_raw",
                   "ci_upper_adjusted", "ci_upper_raw"), skip_absent = TRUE)
  if (!all(c("observed_adjusted", "observed_raw") %in% names(wide))) return(list())
  wide <- wide[which(!is.na(wide$observed_adjusted) & !is.na(wide$observed_raw))]
  if (nrow(wide) < 2) {
    cat(sprintf("  [%s/%s/%s] fewer than 2 PGS with both adjusted+raw versions - skipping\n",
                trait_name, metric_name, type_name))
    return(list())
  }
  wide[, delta := observed_adjusted - observed_raw]
  tt <- tryCatch(stats::t.test(wide$observed_adjusted, wide$observed_raw, paired = TRUE), error = function(e) NULL)
  cat(sprintf("  [%s/%s/%s] n=%d pairs, mean delta (adjusted-raw) = %.4f%s\n", trait_name, metric_name,
              type_name, nrow(wide), mean(wide$delta),
              if (!is.null(tt)) sprintf(", paired t-test p = %.4g", tt$p.value) else ""))

  plots <- list()
  plots$scatter <- ggplot(wide, aes(x = observed_raw, y = observed_adjusted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = ci_lower_adjusted, ymax = ci_upper_adjusted), width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = ci_lower_raw, xmax = ci_upper_raw), height = 0, alpha = 0.3) +
    geom_point(size = 2) +
    labs(title = sprintf("%s — %s: %s (%s) — adjusted vs. raw", opts$cohort_name, trait_name, metric_name, type_name),
         x = sprintf("%s, raw", metric_name), y = sprintf("%s, ancestry-adjusted", metric_name)) +
    axis_coord(c(wide$observed_raw, wide$ci_lower_raw, wide$ci_upper_raw),
               c(wide$observed_adjusted, wide$ci_lower_adjusted, wide$ci_upper_adjusted),
               opts$max_axis_ratio) +
    theme_minimal()

  plots$slope <- ggplot(d, aes(x = group, y = .data[[obs_col]], group = pgs_id)) +
    geom_line(alpha = 0.3) + geom_point(aes(color = group), size = 2) +
    scale_x_discrete(limits = c("raw", "adjusted")) +
    labs(title = sprintf("%s — %s: %s (%s) — per-PGS change, raw → adjusted", opts$cohort_name, trait_name, metric_name, type_name),
         x = NULL, y = metric_name, color = NULL) +
    theme_minimal() + theme(legend.position = "none")

  plots$delta <- ggplot(wide, aes(x = delta)) +
    geom_histogram(bins = min(20, nrow(wide)), fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = mean(wide$delta), color = "firebrick") +
    labs(title = sprintf("%s — %s: %s (%s) — distribution of (adjusted − raw)", opts$cohort_name, trait_name, metric_name, type_name),
         x = sprintf("%s, adjusted − raw (red line = mean)", metric_name), y = "Number of PGS") +
    theme_minimal()
  plots
}

plot_own_vs_catalog <- function(master, catalog_perf, trait_name, mapping, opts) {
  obs_col <- metric_col_explicit(mapping$own_metric, mapping$own_type, "observed")
  lo_col  <- metric_col_explicit(mapping$own_metric, mapping$own_type, "ci_lower")
  hi_col  <- metric_col_explicit(mapping$own_metric, mapping$own_type, "ci_upper")
  if (!(obs_col %in% names(master)) || is.null(catalog_perf) || nrow(catalog_perf) == 0 ||
      !("ancestry_adjusted" %in% names(master))) return(NULL)

  sel <- which(master$trait == trait_name & !is.na(master[[obs_col]]))
  if (length(sel) == 0) return(NULL)
  own_d <- data.table(pgs = master$pgs[sel], pgs_id = master$pgs_id[sel],
                      ancestry_adjusted = master$ancestry_adjusted[sel],
                      estimate = master[[obs_col]][sel], ci_lower = master[[lo_col]][sel],
                      ci_upper = master[[hi_col]][sel], source = "Yours")

  csel <- which(grepl(mapping$catalog_pattern, catalog_perf$name_short, ignore.case = TRUE) &
                  catalog_perf$pgs_id %in% own_d$pgs_id)
  if (length(csel) == 0) {
    cat(sprintf("  [%s] no catalog matches for pattern '%s' among your PGS - skipping\n",
                mapping$label, mapping$catalog_pattern))
    return(NULL)
  }
  cat_d <- catalog_perf[csel, list(pgs_id, publication, estimate, ci_lower, ci_upper)]
  cat_d[, source := paste0("Catalog (", publication, ")")]
  cat_d <- merge(cat_d, unique(own_d[, list(pgs_id, pgs, ancestry_adjusted)]),
                 by = "pgs_id", allow.cartesian = TRUE)

  cols <- c("pgs", "pgs_id", "ancestry_adjusted", "source", "estimate", "ci_lower", "ci_upper")
  plot_dt <- rbindlist(list(own_d[, cols, with = FALSE], cat_d[, cols, with = FALSE]), fill = TRUE)
  plot_dt[, pgs := factor(pgs, levels = own_d$pgs[order(own_d$estimate)])]

  ggplot(plot_dt, aes(x = pgs, y = estimate, color = source, shape = ancestry_adjusted)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width = 0.6), width = 0, alpha = 0.6) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    labs(title = sprintf("%s — %s: %s — yours vs. PGS Catalog", opts$cohort_name, trait_name, mapping$label),
         x = "PGS", y = mapping$label, color = "Source", shape = "Ancestry-adjusted") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), legend.position = "bottom")
}

# Pooled OR per group (adjusted/raw): fixed-effect meta-analysis on the log-OR
# scale; SE recovered from each PGS's stored 95% CI; I^2 reported.
plot_or_pooled <- function(master, trait_name, opts) {
  obs_col <- metric_col("or", "observed", opts)
  lo_col  <- metric_col("or", "ci_lower", opts)
  hi_col  <- metric_col("or", "ci_upper", opts)
  if (!(obs_col %in% names(master)) || !("ancestry_adjusted" %in% names(master))) return(NULL)

  d <- master[which(master$trait == trait_name & !is.na(master[[obs_col]]) & !is.na(master$ancestry_adjusted))]
  if (nrow(d) == 0) return(NULL)
  d <- copy(d)
  d[, group := ifelse(ancestry_adjusted, "adjusted", "raw")]
  d[, log_or := log(d[[obs_col]])]
  d[, se := (log(d[[hi_col]]) - log(d[[lo_col]])) / (2 * 1.96)]
  d <- d[which(is.finite(d$log_or) & is.finite(d$se) & d$se > 0)]
  if (nrow(d) == 0) return(NULL)

  pooled <- d[, {
    w <- 1 / se^2
    pooled_log_or <- sum(w * log_or) / sum(w)
    pooled_se <- sqrt(1 / sum(w))
    k <- .N
    Q <- sum(w * (log_or - pooled_log_or)^2)
    i2 <- if (k > 1) max(0, (Q - (k - 1)) / Q) * 100 else NA_real_
    list(pooled_or = exp(pooled_log_or), ci_lower = exp(pooled_log_or - 1.96 * pooled_se),
         ci_upper = exp(pooled_log_or + 1.96 * pooled_se), k = k, i2 = i2)
  }, by = group]

  caption_txt <- paste(sprintf("%s: k=%d PGS, I²=%s", pooled$group, pooled$k,
                               ifelse(is.na(pooled$i2), "NA", sprintf("%.0f%%", pooled$i2))), collapse = "    ")
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_jitter(data = d, aes(x = group, y = .data[[obs_col]]), width = 0.15, size = 1, alpha = 0.25, color = "grey40") +
    geom_errorbar(data = pooled, aes(x = group, ymin = ci_lower, ymax = ci_upper), width = 0.1, linewidth = 0.8, color = "firebrick") +
    geom_point(data = pooled, aes(x = group, y = pooled_or), size = 4, shape = 18, color = "firebrick") +
    scale_y_log10() +
    labs(title = sprintf("%s — %s: pooled Odds Ratio per 1-SD PRS increase", opts$cohort_name, trait_name),
         subtitle = "Fixed-effect meta-analysis across all PGS per group (red diamond); grey points = individual PGS",
         x = NULL, y = "Odds Ratio (log scale)", caption = caption_txt) +
    theme_minimal()
}

plot_or_forest <- function(master, trait_name, opts) {
  obs_col <- metric_col("or", "observed", opts)
  lo_col  <- metric_col("or", "ci_lower", opts)
  hi_col  <- metric_col("or", "ci_upper", opts)
  if (!(obs_col %in% names(master)) || !("ancestry_adjusted" %in% names(master))) return(NULL)
  d <- master[which(master$trait == trait_name & !is.na(master[[obs_col]]))]
  if (nrow(d) == 0) return(NULL)
  mean_or <- mean(d[[obs_col]], na.rm = TRUE)
  d <- copy(d)
  d[, pgs := factor(pgs, levels = pgs[order(d[[obs_col]])])]

  ggplot(d, aes(x = pgs, y = .data[[obs_col]], color = ancestry_adjusted)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = mean_or, color = "firebrick", linetype = "dotted") +
    geom_errorbar(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), width = 0, alpha = 0.4) +
    geom_point(size = 1.3) + coord_flip() +
    labs(title = sprintf("%s — %s: Odds Ratio per PGS (mean = %.2f, dotted line)", opts$cohort_name, trait_name, mean_or),
         x = NULL, y = "Odds Ratio (per 1-SD increase, 95% CI)", color = "Ancestry-adjusted") +
    theme_minimal() + theme(axis.text.y = element_text(size = 5))
}

plot_or_individual_predictions <- function(master, trait_name, results_dir, prefix, opts) {
  pred_path <- file.path(results_dir, sprintf("%s_%s_or_predictions.csv", prefix, trait_name))
  if (!file.exists(pred_path)) {
    cat(sprintf("  [%s] no individual OR-prediction file at %s - skipping\n", trait_name, pred_path))
    return(NULL)
  }
  preds <- fread(pred_path)
  sel <- which(master$trait == trait_name)
  lookup <- unique(data.table(pgs = master$pgs[sel], ancestry_adjusted = master$ancestry_adjusted[sel]))
  preds <- merge(preds, lookup, by = "pgs", all.x = TRUE)

  ind_order <- preds[, list(m = mean(predicted_prob, na.rm = TRUE)), by = IID][order(m)]$IID
  preds[, IID := factor(IID, levels = ind_order)]
  preds[, outcome_label := ifelse(actual_outcome == 1, "Case", "Control")]

  ggplot(preds, aes(x = IID, y = predicted_prob)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.15, color = "grey40") +
    geom_jitter(aes(color = ancestry_adjusted), width = 0.15, size = 1, alpha = 0.5) +
    facet_wrap(~ outcome_label, scales = "free_x") +
    labs(title = sprintf("%s — %s: predicted risk per individual, across all %d PGS",
                         opts$cohort_name, trait_name, length(unique(preds$pgs))),
         x = "Sampled individual", y = "Predicted probability", color = "Ancestry-adjusted") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
