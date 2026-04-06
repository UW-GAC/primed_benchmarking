context("Two-way ancestry summary")

# Helper: build a minimal admixture data frame --------------------------------

make_adm <- function(k1, k2) {
    data.frame(KAFR = k1, KEUR = k2, stringsAsFactors = FALSE)
}

# Tests for input validation ---------------------------------------------------

test_that("get_two_way_ancestry errors on non-character cohort_name", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, cohort_name = 123),
        "'cohort_name' must be a scalar character string"
    )
})

test_that("get_two_way_ancestry errors on length-2 cohort_name", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, cohort_name = c("a", "b")),
        "'cohort_name' must be a scalar character string"
    )
})

test_that("get_two_way_ancestry errors on NA cohort_name", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, cohort_name = NA_character_),
        "'cohort_name' must be a scalar character string"
    )
})

test_that("get_two_way_ancestry errors when threshold is out of range", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, "cohort", threshold = 1.5),
        "'threshold' must be a numeric value in"
    )
    expect_error(
        get_two_way_ancestry(df, "cohort", threshold = -0.1),
        "'threshold' must be a numeric value in"
    )
})

test_that("get_two_way_ancestry errors when threshold is non-numeric", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, "cohort", threshold = "high"),
        "'threshold' must be a numeric value in"
    )
})

test_that("get_two_way_ancestry errors when min_prop is out of range", {
    df <- make_adm(0.9, 0.05)
    expect_error(
        get_two_way_ancestry(df, "cohort", min_prop = 1.1),
        "'min_prop' must be a numeric value in"
    )
    expect_error(
        get_two_way_ancestry(df, "cohort", min_prop = -0.05),
        "'min_prop' must be a numeric value in"
    )
})

test_that("get_two_way_ancestry errors when fewer than two K columns present", {
    df <- data.frame(KAFR = c(0.9, 0.1), other = c(0.1, 0.9))
    expect_error(
        get_two_way_ancestry(df, "cohort"),
        "at least two columns"
    )
})

test_that("get_two_way_ancestry errors when no K columns present", {
    df <- data.frame(pop1 = c(0.9, 0.1), pop2 = c(0.1, 0.9))
    expect_error(
        get_two_way_ancestry(df, "cohort"),
        "at least two columns"
    )
})

# Tests for correct output structure -------------------------------------------

test_that("get_two_way_ancestry returns a tibble with expected columns", {
    df <- make_adm(c(0.5, 0.9, 0.05), c(0.5, 0.05, 0.93))
    result <- get_two_way_ancestry(df, cohort_name = "test_cohort")

    expected_cols <- c(
        "Cohort", "Ref_Pop1", "Ref_Pop2",
        "Count_two_way",
        "Excluded_Ref1_lt10_and_Ref2_lt90",
        "Excluded_Ref2_lt10_and_Ref1_lt90",
        "Excluded_Ref1_lt10_and_Ref2_gt90",
        "Excluded_Ref2_lt10_and_Ref1_gt90",
        "n"
    )
    expect_true(all(expected_cols %in% colnames(result)))
    expect_true(inherits(result, "tbl_df"))
})

test_that("get_two_way_ancestry returns one row per pair for two K columns", {
    df <- make_adm(c(0.5, 0.9), c(0.5, 0.05))
    result <- get_two_way_ancestry(df, cohort_name = "cohort")
    expect_equal(nrow(result), 1L)
})

test_that("get_two_way_ancestry returns three rows for three K columns", {
    df <- data.frame(
        K1 = c(0.8, 0.1, 0.1),
        K2 = c(0.1, 0.8, 0.1),
        K3 = c(0.1, 0.1, 0.8)
    )
    result <- get_two_way_ancestry(df, cohort_name = "cohort")
    expect_equal(nrow(result), 3L)
    expect_equal(sort(result$Ref_Pop1), c("K1", "K1", "K2"))
    expect_equal(sort(result$Ref_Pop2), c("K2", "K3", "K3"))
})

# Tests for correct counts -----------------------------------------------------

test_that("get_two_way_ancestry correctly counts two-way admixed individuals", {
    # Row 1: KAFR=0.5, KEUR=0.5 -> two_way (both >= 0.10, sum >= 0.9)
    # Row 2: KAFR=0.95, KEUR=0.03 -> not two_way (x2 < min_prop)
    # Row 3: KAFR=0.03, KEUR=0.95 -> not two_way (x1 < min_prop)
    df <- make_adm(c(0.5, 0.95, 0.03), c(0.5, 0.03, 0.95))
    result <- get_two_way_ancestry(df, cohort_name = "cohort",
                                   threshold = 0.9, min_prop = 0.10)
    expect_equal(result$Count_two_way, 1L)
})

test_that("get_two_way_ancestry correctly identifies mono_ref1 (Ref2 gt90)", {
    # KAFR=0.02, KEUR=0.95 -> mono_ref1: x1 < min_prop, x2 >= threshold
    df <- make_adm(0.02, 0.95)
    result <- get_two_way_ancestry(df, cohort_name = "cohort")
    expect_equal(result$Excluded_Ref1_lt10_and_Ref2_gt90, 1L)
    expect_equal(result$Count_two_way, 0L)
})

test_that("get_two_way_ancestry correctly identifies mono_ref2 (Ref1 gt90)", {
    # KAFR=0.95, KEUR=0.02 -> mono_ref2: x1 >= threshold, x2 < min_prop
    df <- make_adm(0.95, 0.02)
    result <- get_two_way_ancestry(df, cohort_name = "cohort")
    expect_equal(result$Excluded_Ref2_lt10_and_Ref1_gt90, 1L)
    expect_equal(result$Count_two_way, 0L)
})

test_that("get_two_way_ancestry populates Cohort, Ref_Pop1, Ref_Pop2 columns", {
    df <- make_adm(c(0.5, 0.9), c(0.5, 0.05))
    result <- get_two_way_ancestry(df, cohort_name = "my_cohort")
    expect_equal(result$Cohort, "my_cohort")
    expect_equal(result$Ref_Pop1, "KAFR")
    expect_equal(result$Ref_Pop2, "KEUR")
})

test_that("get_two_way_ancestry n equals individuals with sum >= threshold", {
    # 3 individuals; sum >= 0.9 for rows 1 and 2 only
    df <- make_adm(c(0.5, 0.95, 0.2), c(0.45, 0.03, 0.3))
    result <- get_two_way_ancestry(df, cohort_name = "cohort",
                                   threshold = 0.9, min_prop = 0.10)
    # Row 1: 0.5 + 0.45 = 0.95 >= 0.9 yes
    # Row 2: 0.95 + 0.03 = 0.98 >= 0.9 yes
    # Row 3: 0.2 + 0.3 = 0.5 < 0.9 no
    expect_equal(result$n, 2L)
})

test_that("get_two_way_ancestry ignores non-K columns", {
    df <- data.frame(
        sample_id = 1:3,
        KAFR      = c(0.5, 0.9, 0.03),
        KEUR      = c(0.4, 0.03, 0.92),
        other_col = c(0.1, 0.07, 0.05)
    )
    result <- get_two_way_ancestry(df, cohort_name = "cohort")
    expect_equal(nrow(result), 1L)  # only one pair: KAFR-KEUR
    expect_equal(result$Ref_Pop1, "KAFR")
    expect_equal(result$Ref_Pop2, "KEUR")
})
