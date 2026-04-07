#' Compute two-way ancestry counts from admixture proportion data
#'
#' Given a data frame or tibble of individual-level admixture proportions with
#' columns whose names begin with \code{"K"} (one column per reference
#' population), this function counts the number of individuals that satisfy
#' two-way ancestry criteria for every pair of reference populations.
#'
#' @param admixture_anc_prop_list A data frame or tibble containing at least
#'   two numeric columns whose names begin with \code{"K"} (e.g. \code{K1},
#'   \code{K2}, \code{KAFR}, \code{KEUR}). Each row represents one individual
#'   and each \code{K} column holds that individual's estimated proportion of
#'   ancestry from the corresponding reference population. Other columns are
#'   ignored.
#' @param cohort_name Character scalar. Name of the cohort; used to populate
#'   the \code{Cohort} column in the returned tibble.
#' @param threshold Numeric scalar in \eqn{[0, 1]}. Minimum combined ancestry
#'   proportion (\code{x1 + x2}) required for an individual to be counted in
#'   any category for a given pair. Default is \code{0.9}.
#' @param min_prop Numeric scalar in \eqn{[0, 1]}. Minimum individual ancestry
#'   proportion required for an individual to be classified as admixed (i.e.
#'   having meaningful ancestry from both reference populations). Default is
#'   \code{0.10}.
#'
#' @return A tibble with one row per pair of reference populations and the
#'   following columns:
#'   \describe{
#'     \item{\code{Cohort}}{Cohort name (from \code{cohort_name}).}
#'     \item{\code{Ref_Pop1}}{Name of the first reference population column.}
#'     \item{\code{Ref_Pop2}}{Name of the second reference population column.}
#'     \item{\code{Count_two_way}}{Number of individuals with
#'       \code{x1 >= min_prop}, \code{x2 >= min_prop}, and
#'       \code{x1 + x2 >= threshold} (genuinely admixed between the two
#'       populations).}
#'     \item{\code{Excluded_Ref1_lt10_and_Ref2_lt90}}{Number of individuals
#'       with \code{x1 < min_prop}, \code{x2 < threshold}, and
#'       \code{x1 + x2 >= threshold}.}
#'     \item{\code{Excluded_Ref2_lt10_and_Ref1_lt90}}{Number of individuals
#'       with \code{x2 < min_prop}, \code{x1 < threshold}, and
#'       \code{x1 + x2 >= threshold}.}
#'     \item{\code{Excluded_Ref1_lt10_and_Ref2_gt90}}{Number of individuals
#'       with \code{x1 < min_prop} and \code{x2 >= threshold} (predominantly
#'       Ref2 ancestry).}
#'     \item{\code{Excluded_Ref2_lt10_and_Ref1_gt90}}{Number of individuals
#'       with \code{x2 < min_prop} and \code{x1 >= threshold} (predominantly
#'       Ref1 ancestry).}
#'     \item{\code{n}}{Total number of individuals with
#'       \code{x1 + x2 >= threshold}.}
#'   }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   sample_id = 1:5,
#'   KAFR      = c(0.95, 0.05, 0.50, 0.02, 0.80),
#'   KEUR      = c(0.03, 0.93, 0.48, 0.03, 0.15)
#' )
#' get_two_way_ancestry(df, cohort_name = "my_cohort")
#' }
#'
#' @importFrom dplyr select starts_with bind_rows
#' @importFrom tibble tibble
#' @export
get_two_way_ancestry <- function(admixture_anc_prop_list,
                                 cohort_name,
                                 threshold = 0.9,
                                 min_prop  = 0.10) {
    if (!is.character(cohort_name) || length(cohort_name) != 1 ||
            is.na(cohort_name)) {
        stop("'cohort_name' must be a scalar character string.")
    }
    if (!is.numeric(threshold) || length(threshold) != 1 ||
            is.na(threshold) || threshold < 0 || threshold > 1) {
        stop("'threshold' must be a numeric value in [0, 1].")
    }
    if (!is.numeric(min_prop) || length(min_prop) != 1 ||
            is.na(min_prop) || min_prop < 0 || min_prop > 1) {
        stop("'min_prop' must be a numeric value in [0, 1].")
    }

    adm <- dplyr::select(admixture_anc_prop_list, dplyr::starts_with("K"))

    if (ncol(adm) < 2) {
        stop(
            "'admixture_anc_prop_list' must contain at least two columns ",
            "whose names begin with 'K'."
        )
    }

    pairs <- as.data.frame(
        t(combn(colnames(adm), 2)),
        stringsAsFactors = FALSE
    )
    colnames(pairs) <- c("Ref1", "Ref2")

    out_list <- vector("list", nrow(pairs))

    for (i in seq_len(nrow(pairs))) {
        x1 <- adm[[pairs$Ref1[i]]]
        x2 <- adm[[pairs$Ref2[i]]]

        two_way      <- (x1 >= min_prop) & (x2 >= min_prop) &
                            (x1 + x2 >= threshold)

        if(sum(two_way)>100) {
          print(paste("Count of two_way for",pairs$Ref1[i],"and",pairs$Ref2[i]," is ",sum(two_way)))
          print(paste0("Writing file: ", cohort,"_",pairs$Ref1[i],"_",pairs$Ref2[i],"_two_way_anc.prop_IDs.txt"))
          write_delim(admixture_anc_prop_list %>% filter(two_way) %>% select(ID),
                      file=paste0("Writing file: ", cohort,"_",pairs$Ref1[i],"_",pairs$Ref2[i],"_two_way_anc.prop_IDs.txt"))
          }
      
        # x1 < min_prop and x2 < threshold, but their sum still reaches
        # threshold (intermediate state: Ref1 very low, Ref2 moderate-high)
        exclude_ref1 <- (x1 < min_prop)  & (x2 < threshold) &
                            (x1 + x2 >= threshold)
        # x2 < min_prop and x1 < threshold, but their sum still reaches
        # threshold (intermediate state: Ref2 very low, Ref1 moderate-high)
        exclude_ref2 <- (x1 < threshold) & (x2 < min_prop)  &
                            (x1 + x2 >= threshold)
        # x1 < min_prop and x2 >= threshold (predominantly Ref2 ancestry);
        # the sum check is redundant here but kept for consistency with spec
        mono_ref1    <- (x1 < min_prop)  & (x2 >= threshold) &
                            (x1 + x2 >= threshold)
        # x2 < min_prop and x1 >= threshold (predominantly Ref1 ancestry);
        # the sum check is redundant here but kept for consistency with spec
        mono_ref2    <- (x1 >= threshold) & (x2 < min_prop)  &
                            (x1 + x2 >= threshold)
        total        <- (x1 + x2 >= threshold)
        if(sum(total)>100) {
          print(paste("Count of total n for",pairs$Ref1[i],"and",pairs$Ref2[i]," is ",sum(total)))
          print(paste0("Writing file: ", cohort,"_",pairs$Ref1[i],"_",pairs$Ref2[i],"_total_n_anc.prop_IDs.txt"))
          write_delim(admixture_anc_prop_list %>% filter(total) %>% select(ID),
                      file=paste0("Writing file: ", cohort,"_",pairs$Ref1[i],"_",pairs$Ref2[i],"_total_n_anc.prop_IDs.txt"))
          }

        out_list[[i]] <- tibble::tibble(
            Cohort                           = cohort_name,
            Ref_Pop1                         = pairs$Ref1[i],
            Ref_Pop2                         = pairs$Ref2[i],
            Count_two_way                    = sum(two_way),
            Excluded_Ref1_lt10_and_Ref2_lt90 = sum(exclude_ref1),
            Excluded_Ref2_lt10_and_Ref1_lt90 = sum(exclude_ref2),
            Excluded_Ref1_lt10_and_Ref2_gt90 = sum(mono_ref1),
            Excluded_Ref2_lt10_and_Ref1_gt90 = sum(mono_ref2),
            n                                = sum(total)
        )
    }

    dplyr::bind_rows(out_list)
}
