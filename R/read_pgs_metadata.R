#' Read and reformat the PGS Catalog "all metadata scores" table
#'
#' Downloads and cleans the PGS Catalog "all metadata scores" CSV file from the
#' EBI FTP server. The raw CSV is read as all-character, then column names are
#' standardized to \code{snake_case}, whitespace is trimmed, obviously numeric
#' or date-like columns are converted to their native types, and optionally
#' columns whose values are pipe/semicolon/comma-separated lists are converted
#' to list-columns.
#'
#' @param url Character. URL of the CSV file.
#'   Defaults to the canonical EBI FTP location.
#' @param split_list_columns Logical. If \code{TRUE} (default), columns whose
#'   values frequently contain pipe (\code{|}), semicolon (\code{;}), or
#'   comma-space (\code{,\\ }) delimiters are split into list-columns using
#'   \code{\link[stringr]{str_split}}.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} with cleaned columns ready for
#'   analysis and joins. When \code{split_list_columns = TRUE}, some columns
#'   will be list-columns (each element a character vector).
#'
#' @examples
#' \dontrun{
#' pgs <- read_pgs_all_metadata_scores()
#' dplyr::glimpse(pgs)
#' dplyr::count(pgs, pgs_id, sort = TRUE)
#' }
#'
#' @importFrom readr read_csv cols col_character parse_double parse_date
#' @importFrom dplyr mutate across coalesce any_of all_of if_else
#' @importFrom stringr str_squish str_detect str_replace str_split str_to_lower
#' @importFrom janitor clean_names
#' @importFrom purrr map_lgl
#' @export
read_pgs_all_metadata_scores <- function(
    url = paste0(
        "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/",
        "pgs_all_metadata_scores.csv"
    ),
    split_list_columns = TRUE
) {
    # Read everything as character to avoid accidental parsing surprises
    raw <- readr::read_csv(
        url,
        col_types = readr::cols(.default = readr::col_character()),
        na = c("", "NA", "N/A", "NULL")
    )
    raw <- janitor::clean_names(raw)
    raw <- dplyr::mutate(
        raw,
        dplyr::across(dplyr::everything(),
                      ~ dplyr::if_else(is.na(.x), NA_character_,
                                       stringr::str_squish(.x)))
    )

    # Normalise the primary PGS identifier if present.
    # Prefer pgs_id if it exists; fall back to pgs or score_id.
    id_candidates <- intersect(c("pgs_id", "pgs", "score_id"), names(raw))
    if (length(id_candidates) > 0) {
        raw <- dplyr::mutate(
            raw,
            pgs_id = dplyr::coalesce(!!!lapply(id_candidates,
                                               function(n) raw[[n]])),
            pgs_id = dplyr::if_else(
                !is.na(.data$pgs_id),
                stringr::str_replace(.data$pgs_id, "^(?!PGS)(.+)$", "\\1"),
                .data$pgs_id
            )
        )
    }

    # Convert columns that look fully numeric or fully date (YYYY-MM-DD).
    # We only attempt conversion on character columns to avoid double-conversion.
    char_cols <- names(raw)[vapply(raw, is.character, logical(1))]

    numeric_cols <- char_cols[vapply(char_cols, function(cn) {
        .is_numeric_col(raw[[cn]])
    }, logical(1))]

    date_cols <- setdiff(char_cols, numeric_cols)[vapply(
        setdiff(char_cols, numeric_cols), function(cn) {
            .is_date_col(raw[[cn]])
        }, logical(1))]

    if (length(numeric_cols) > 0) {
        raw <- dplyr::mutate(
            raw,
            dplyr::across(
                dplyr::all_of(numeric_cols),
                ~ readr::parse_double(.x, na = c("", "NA", "N/A", "NULL"))
            )
        )
    }

    if (length(date_cols) > 0) {
        raw <- dplyr::mutate(
            raw,
            dplyr::across(
                dplyr::all_of(date_cols),
                ~ readr::parse_date(.x, format = "%Y-%m-%d",
                                    na = c("", "NA", "N/A", "NULL"))
            )
        )
    }

    # Convenience normalizations for common columns — done while they are still
    # character (before optional list-col splitting).
    raw <- dplyr::mutate(
        raw,
        dplyr::across(
            dplyr::any_of(c("trait_reported", "trait_mapped", "trait_efo")),
            stringr::str_squish
        ),
        dplyr::across(
            dplyr::any_of(c("ancestry_distribution", "ancestry_broad")),
            stringr::str_to_lower
        ),
        dplyr::across(
            dplyr::any_of(c("license", "use_restriction")),
            stringr::str_squish
        )
    )

    # Optional: split list-like character columns into list-columns.
    if (isTRUE(split_list_columns)) {
        remaining_char <- names(raw)[vapply(raw, is.character, logical(1))]
        list_cols <- remaining_char[purrr::map_lgl(
            remaining_char,
            ~ .is_listlike_col(raw[[.x]])
        )]

        if (length(list_cols) > 0) {
            raw <- dplyr::mutate(
                raw,
                dplyr::across(
                    dplyr::all_of(list_cols),
                    ~ stringr::str_split(.x, pattern = "\\s*(\\||;|,)\\s*")
                )
            )
        }
    }

    raw
}


# ---- Internal helpers --------------------------------------------------------

#' @keywords internal
.is_numeric_col <- function(x) {
    v <- x[!is.na(x)]
    if (length(v) == 0L) return(FALSE)
    all(stringr::str_detect(
        v,
        "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"
    ))
}

#' @keywords internal
.is_date_col <- function(x) {
    v <- x[!is.na(x)]
    if (length(v) == 0L) return(FALSE)
    all(stringr::str_detect(v, "^\\d{4}-\\d{2}-\\d{2}$"))
}

#' @keywords internal
.is_listlike_col <- function(x) {
    v <- x[!is.na(x)]
    if (length(v) == 0L) return(FALSE)
    mean(stringr::str_detect(v, "\\||;|,\\s")) >= 0.10
}
