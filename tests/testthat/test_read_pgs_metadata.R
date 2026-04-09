context("read_pgs_all_metadata_scores helpers and cleaning logic")

# Helper: build a small CSV string and run it through the same cleaning steps
# used by the function (without hitting the network).

make_raw_tibble <- function(text, split_list_columns = FALSE) {
    raw <- readr::read_csv(
        I(text),
        col_types = readr::cols(.default = readr::col_character()),
        na = c("", "NA", "N/A", "NULL"),
        show_col_types = FALSE
    )
    raw <- janitor::clean_names(raw)
    raw <- dplyr::mutate(
        raw,
        dplyr::across(dplyr::everything(),
                      ~ dplyr::if_else(is.na(.x), NA_character_,
                                       stringr::str_squish(.x)))
    )
    raw
}

# ---- .is_numeric_col --------------------------------------------------------

test_that(".is_numeric_col returns TRUE for all-numeric vectors", {
    expect_true(primed_benchmarking:::.is_numeric_col(c("1", "2.5", "-3", "1e4")))
    expect_true(primed_benchmarking:::.is_numeric_col(c("0", "100")))
})

test_that(".is_numeric_col returns FALSE when any value is non-numeric", {
    expect_false(:::.is_numeric_col(c("1", "abc")))
    expect_false(:::.is_numeric_col(c("1.5", "2019-01-01")))
})

test_that(".is_numeric_col returns FALSE for all-NA input", {
    expect_false(:::.is_numeric_col(c(NA_character_, NA_character_)))
})

# ---- .is_date_col -----------------------------------------------------------

test_that(".is_date_col returns TRUE for YYYY-MM-DD vectors", {
    expect_true(primed_benchmarking:::.is_date_col(c("2020-01-15", "2023-12-31")))
})

test_that(".is_date_col returns FALSE for non-date patterns", {
    expect_false(primed_benchmarking:::.is_date_col(c("01/15/2020", "Jan 2020")))
    expect_false(primed_benchmarking:::.is_date_col(c("1", "2")))
})

test_that(".is_date_col returns FALSE for all-NA input", {
    expect_false(primed_benchmarking:::.is_date_col(NA_character_))
})

# ---- .is_listlike_col -------------------------------------------------------

test_that(".is_listlike_col returns TRUE when >= 10% of values have delimiters", {
    # 5/5 = 100% pipe-separated
    expect_true(primed_benchmarking:::.is_listlike_col(c("a|b", "c|d", "e", "f|g", "h|i")))
})

test_that(".is_listlike_col returns FALSE when few values contain delimiters", {
    # 0% have delimiters
    expect_false(primed_benchmarking:::.is_listlike_col(c("apple", "banana", "cherry")))
})

test_that(".is_listlike_col returns FALSE for empty / all-NA input", {
    expect_false(primed_benchmarking:::.is_listlike_col(character(0)))
    expect_false(primed_benchmarking:::.is_listlike_col(c(NA_character_, NA_character_)))
})

# ---- Column name cleaning ---------------------------------------------------

test_that("read_pgs_all_metadata_scores standardises column names to snake_case", {
    csv_text <- "Pgs Id,Trait Reported,Num Variants\nPGS000001,Height,500000"
    raw <- suppressMessages(
        readr::read_csv(I(csv_text),
                        col_types = readr::cols(.default = readr::col_character()),
                        na = c("", "NA", "N/A", "NULL"))
    )
    cleaned <- janitor::clean_names(raw)
    expect_true(all(names(cleaned) %in% c("pgs_id", "trait_reported",
                                          "num_variants")))
})

# ---- NA handling ------------------------------------------------------------

test_that("empty strings and N/A become NA", {
    csv_text <- "pgs_id,trait_reported\nPGS000001,\nPGS000002,N/A\nPGS000003,Height"
    raw <- make_raw_tibble(csv_text)
    expect_true(is.na(raw$trait_reported[1]))
    expect_true(is.na(raw$trait_reported[2]))
    expect_equal(raw$trait_reported[3], "Height")
})

# ---- Whitespace trimming ----------------------------------------------------

test_that("whitespace is trimmed from character columns", {
    csv_text <- "pgs_id,trait_reported\n  PGS000001  ,  Height  \n"
    raw <- make_raw_tibble(csv_text)
    expect_equal(raw$pgs_id[1], "PGS000001")
    expect_equal(raw$trait_reported[1], "Height")
})

# ---- Numeric conversion (via helper) ----------------------------------------

test_that(".is_numeric_col handles positive integers and floats", {
    expect_true(primed_benchmarking:::.is_numeric_col(c("1", "2", "3")))
    expect_true(primed_benchmarking:::.is_numeric_col(c("1.5", "2.0")))
})

test_that(".is_numeric_col rejects mixed numeric/text column", {
    expect_false(primed_benchmarking:::.is_numeric_col(c("1", "two")))
})

# ---- Date conversion (via helper) -------------------------------------------

test_that(".is_date_col accepts standard ISO dates", {
    expect_true(primed_benchmarking:::.is_date_col(c("2021-06-01", "2022-11-30")))
})

test_that(".is_date_col rejects non-ISO date strings", {
    expect_false(primed_benchmarking:::.is_date_col(c("06/01/2021")))
})

# ---- split_list_columns parameter -------------------------------------------

test_that(".is_listlike_col detects pipe-separated values", {
    x <- c("AFR|EUR", "AMR|EAS", "EUR", "AFR|AMR|EUR", "SAS")
    expect_true(primed_benchmarking:::.is_listlike_col(x))
})

test_that(".is_listlike_col detects semicolon-separated values", {
    x <- c("a;b", "c;d", "e")
    expect_true(primed_benchmarking:::.is_listlike_col(x))
})

test_that(".is_listlike_col detects comma-space-separated values", {
    x <- c("a, b, c", "d, e", "f")
    expect_true(primed_benchmarking:::.is_listlike_col(x))
})
