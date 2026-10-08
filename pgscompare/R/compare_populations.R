#' Compare pgsmetrics results between two runs on different individuals
#'
#' For each entry in `population_comparisons`, reads
#' `<cohort>[_<tag>]_<trait>_pgsmetrics.csv` for both sides, matches PGS by
#' name and makes a paired scatter, slope plot and delta histogram per metric
#' (only PGS present in both runs), optionally colored by GWAS or development
#' ancestry, with sample sizes in the subtitle.
#'
#' @param cohort_name cohort name used in [run_combine_and_compare()].
#' @param population_comparisons list of `list(label1, tag1, label2, tag2,
#'   traits = NULL)`; a `tag` is a run's subset tag (NULL = no-subset run).
#' @param results_dir local results folder.
#' @param metric_types,default_type,max_axis_ratio,metrics_to_plot see
#'   [run_join_and_plot()].
#' @param color_var `"gwas_ancestry"`, `"dev_ancestry"` or NULL.
#' @param plot_dir root plots folder (default `<results_dir>/plots`).
#' @param bucket_dir `gs://` folder to pull inputs from / push plots to.
#' @param sync_inputs pull from the bucket first?
#' @export
run_compare_populations <- function(
    cohort_name, population_comparisons, results_dir = "pgs_comparison_results",
    metric_types = list(or = "full"), default_type = "partial",
    metrics_to_plot = list(), color_var = "gwas_ancestry", max_axis_ratio = 1.5,
    plot_dir = file.path(results_dir, "plots"),
    bucket_dir = default_bucket_dir(), sync_inputs = TRUE) {

  opts <- list(cohort_name = cohort_name, metric_types = metric_types,
               default_type = default_type, max_axis_ratio = max_axis_ratio)
  plot_dir <- file.path(plot_dir, cohort_name, "population_comparisons")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  if (sync_inputs) sync_from_bucket(bucket_dir, results_dir)

  for (cmp in population_comparisons) {
    traits <- cmp$traits
    if (is.null(traits)) {
      traits <- intersect(detect_traits(results_dir, cohort_name, cmp$tag1),
                          detect_traits(results_dir, cohort_name, cmp$tag2))
      cat(sprintf("Auto-detected shared traits for '%s' vs '%s': %s\n",
                  cmp$label1, cmp$label2, paste(traits, collapse = ", ")))
    }

    for (trait_name in traits) {
      dt1 <- load_population_metrics(results_dir, cohort_name, trait_name, cmp$tag1)
      dt2 <- load_population_metrics(results_dir, cohort_name, trait_name, cmp$tag2)
      if (is.null(dt1) || is.null(dt2)) {
        cat(sprintf("  [%s] could not load both populations for '%s' vs '%s' - skipping\n",
                    trait_name, cmp$label1, cmp$label2))
        next
      }
      tdir <- file.path(plot_dir, trait_name)
      dir.create(tdir, showWarnings = FALSE, recursive = TRUE)

      ancestry_lookup <- load_ancestry_lookup(results_dir, cohort_name, trait_name)

      n1 <- format_n(load_sample_size(results_dir, cohort_name, trait_name, cmp$tag1))
      n2 <- format_n(load_sample_size(results_dir, cohort_name, trait_name, cmp$tag2))
      if (is.na(n1) || is.na(n2)) {
        cat(sprintf("  [%s] sample-size file missing for %s - rerun run_combine_and_compare() for that run (cached chunks make it quick)\n",
                    trait_name, paste(c(if (is.na(n1)) cmp$label1, if (is.na(n2)) cmp$label2), collapse = " and ")))
      }
      subtitle_text <- if (is.na(n1) && is.na(n2)) NULL else
        sprintf("%s: %s   |   %s: %s", cmp$label1, ifelse(is.na(n1), "n unavailable", n1),
                cmp$label2, ifelse(is.na(n2), "n unavailable", n2))
      if (!is.null(subtitle_text)) cat(sprintf("  [%s] %s\n", trait_name, subtitle_text))

      metrics_here <- metrics_to_plot[[trait_name]]
      if (is.null(metrics_here)) {
        metrics_here <- unique(sub("_(full|partial)_observed$", "",
                                   grep("_(full|partial)_observed$", names(dt1), value = TRUE)))
      }
      color_tag <- ancestry_filename_tag(ancestry_lookup, color_var)

      for (metric_name in metrics_here) {
        type_name <- type_for_metric(metric_name, opts)
        plots <- plot_population_comparison(trait_name, metric_name, type_name, dt1, cmp$label1,
                                            dt2, cmp$label2, opts, ancestry_lookup, color_var, subtitle_text)
        for (kind in names(plots)) {
          fname <- sprintf("%s_vs_%s_%s_%s_%s%s.png", sanitize_for_filename(cmp$label1),
                           sanitize_for_filename(cmp$label2), metric_name, type_name, kind, color_tag)
          save_plot(plots[[kind]], file.path(tdir, fname), 7, 6)
        }
      }
    }
  }
  cat("Done comparing populations.\n")
  sync_to_bucket(bucket_dir, results_dir, csv_glob = NULL, plots = TRUE)
  invisible(NULL)
}

load_population_metrics <- function(results_dir, cohort_name, trait_name, tag) {
  path <- file.path(results_dir, sprintf("%s_%s_pgsmetrics.csv", name_prefix(cohort_name, tag), trait_name))
  if (!file.exists(path)) { cat(sprintf("  Population file not found: %s\n", path)); return(NULL) }
  fread(path)
}

load_sample_size <- function(results_dir, cohort_name, trait_name, tag) {
  path <- file.path(results_dir, sprintf("%s_%s_sample_size.csv", name_prefix(cohort_name, tag), trait_name))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

format_n <- function(ss) {
  if (is.null(ss)) return(NA_character_)
  txt <- sprintf("n = %s", format(ss$n_usable, big.mark = ","))
  if (all(c("n_cases", "n_controls") %in% names(ss)) && !is.na(ss$n_cases)) {
    txt <- sprintf("%s (%s cases / %s controls)", txt, format(ss$n_cases, big.mark = ","),
                   format(ss$n_controls, big.mark = ","))
  }
  txt
}

# pgs (unique_name) -> gwas_ancestry / dev_ancestry, from the default run's
# pgs_metadata + PGS Catalog score metadata (ancestry doesn't depend on subset).
load_ancestry_lookup <- function(results_dir, cohort_name, trait_name) {
  meta_path <- resolve_optional(results_dir, cohort_name, "pgs_metadata", trait_name)
  if (is.na(meta_path)) {
    cat(sprintf("  [%s] no pgs_metadata found - ancestry coloring unavailable\n", trait_name))
    return(NULL)
  }
  meta <- fread(meta_path)
  stopifnot(all(c("pgs_id", "unique_name") %in% names(meta)))
  lookup <- data.table(pgs = meta$unique_name, pgs_id = meta$pgs_id)

  cpath <- resolve_trait_or_global(results_dir, sprintf("%s_pgs_catalog_score_metadata.csv", trait_name),
                                   "pgs_catalog_score_metadata.csv")
  if (!is.na(cpath)) {
    cmeta <- fread(cpath)
    keep_cols <- intersect(c("pgs_id", "gwas_ancestry", "dev_ancestry"), names(cmeta))
    lookup <- merge(lookup, cmeta[, keep_cols, with = FALSE], by = "pgs_id", all.x = TRUE)
    for (col in intersect(c("gwas_ancestry", "dev_ancestry"), names(lookup))) {
      lookup[which(lookup[[col]] == ""), (col) := NA_character_]
    }
  } else {
    cat(sprintf("  [%s] no PGS Catalog score metadata found - ancestry coloring unavailable\n", trait_name))
  }
  lookup
}

ancestry_color_available <- function(ancestry_lookup, color_var) {
  !is.null(color_var) && !is.null(ancestry_lookup) && color_var %in% names(ancestry_lookup)
}
ancestry_filename_tag <- function(ancestry_lookup, color_var) {
  if (!ancestry_color_available(ancestry_lookup, color_var)) return("")
  paste0("_by_", sub("_ancestry$", "", color_var))
}

plot_population_comparison <- function(trait_name, metric_name, type_name, dt1, label1, dt2, label2,
                                       opts, ancestry_lookup = NULL, color_var = NULL, subtitle_text = NULL) {
  obs_col <- sprintf("%s_%s_observed", metric_name, type_name)
  lo_col  <- sprintf("%s_%s_ci_lower", metric_name, type_name)
  hi_col  <- sprintf("%s_%s_ci_upper", metric_name, type_name)
  if (!(obs_col %in% names(dt1)) || !(obs_col %in% names(dt2))) return(list())

  pick <- function(dt) {
    ok <- which(!is.na(dt[[obs_col]]))
    data.table(pgs = dt$pgs[ok], observed = dt[[obs_col]][ok],
               ci_lower = dt[[lo_col]][ok], ci_upper = dt[[hi_col]][ok])
  }
  d1 <- pick(dt1); d2 <- pick(dt2)
  wide <- merge(d1, d2, by = "pgs", suffixes = c("_1", "_2"))
  if (nrow(wide) < 2) {
    cat(sprintf("  [%s/%s/%s] fewer than 2 PGS shared between '%s' and '%s' - skipping\n",
                trait_name, metric_name, type_name, label1, label2))
    return(list())
  }
  wide[, delta := observed_2 - observed_1]

  has_color <- ancestry_color_available(ancestry_lookup, color_var)
  if (has_color) {
    wide <- merge(wide, data.table(pgs = ancestry_lookup$pgs, color = ancestry_lookup[[color_var]]),
                  by = "pgs", all.x = TRUE)
  }
  tt <- tryCatch(stats::t.test(wide$observed_2, wide$observed_1, paired = TRUE), error = function(e) NULL)
  cat(sprintf("  [%s/%s/%s] n=%d shared PGS, mean delta (%s − %s) = %.4f%s\n", trait_name, metric_name,
              type_name, nrow(wide), label2, label1, mean(wide$delta),
              if (!is.null(tt)) sprintf(", paired t-test p = %.4g", tt$p.value) else ""))

  plots <- list()
  scatter_aes <- if (has_color) aes(x = observed_1, y = observed_2, color = color) else aes(x = observed_1, y = observed_2)
  plots$scatter <- ggplot(wide, scatter_aes) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = ci_lower_2, ymax = ci_upper_2), width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = ci_lower_1, xmax = ci_upper_1), height = 0, alpha = 0.3) +
    geom_point(size = 2) +
    labs(title = sprintf("%s — %s: %s (%s) — %s vs. %s", opts$cohort_name, trait_name, metric_name, type_name, label1, label2),
         subtitle = subtitle_text, x = sprintf("%s, %s", metric_name, label1),
         y = sprintf("%s, %s", metric_name, label2), color = if (has_color) color_var else NULL) +
    axis_coord(c(wide$observed_1, wide$ci_lower_1, wide$ci_upper_1),
               c(wide$observed_2, wide$ci_lower_2, wide$ci_upper_2), opts$max_axis_ratio) +
    theme_minimal()

  long <- rbindlist(list(
    d1[d1$pgs %in% wide$pgs, list(pgs, population = label1, observed)],
    d2[d2$pgs %in% wide$pgs, list(pgs, population = label2, observed)]))
  if (has_color) long <- merge(long, wide[, list(pgs, color)], by = "pgs", all.x = TRUE)
  slope_aes   <- if (has_color) aes(x = population, y = observed, group = pgs, color = color) else aes(x = population, y = observed, group = pgs)
  slope_line  <- if (has_color) geom_line(alpha = 0.4) else geom_line(alpha = 0.3, color = "grey40")
  slope_point <- if (has_color) geom_point(size = 2, alpha = 0.8) else geom_point(size = 2, color = "grey30")
  plots$slope <- ggplot(long, slope_aes) + slope_line + slope_point +
    scale_x_discrete(limits = c(label1, label2)) +
    labs(title = sprintf("%s — %s: %s (%s) — per-PGS change, %s → %s", opts$cohort_name, trait_name, metric_name, type_name, label1, label2),
         subtitle = subtitle_text, x = NULL, y = metric_name, color = if (has_color) color_var else NULL) +
    theme_minimal()

  delta_aes  <- if (has_color) aes(x = delta, fill = color) else aes(x = delta)
  delta_hist <- if (has_color) geom_histogram(bins = min(20, nrow(wide)), alpha = 0.85) else
    geom_histogram(bins = min(20, nrow(wide)), fill = "steelblue", alpha = 0.7)
  plots$delta <- ggplot(wide, delta_aes) + delta_hist +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = mean(wide$delta), color = "firebrick") +
    labs(title = sprintf("%s — %s: %s (%s) — distribution of (%s − %s)", opts$cohort_name, trait_name, metric_name, type_name, label2, label1),
         subtitle = subtitle_text, x = sprintf("%s, %s − %s (red line = mean)", metric_name, label2, label1),
         y = "Number of shared PGS", fill = if (has_color) color_var else NULL) +
    theme_minimal()
  plots
}
