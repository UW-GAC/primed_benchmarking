context("PGS pipeline validation and helper functions")

# Tests for .validate_pgs_id -------------------------------------------------

test_that(".validate_pgs_id accepts valid IDs", {
    expect_silent(primedtools:::.validate_pgs_id("PGS000001"))
    expect_silent(primedtools:::.validate_pgs_id("PGS999999"))
    expect_silent(primedtools:::.validate_pgs_id("PGS000123"))
})

test_that(".validate_pgs_id rejects invalid IDs", {
    expect_error(primedtools:::.validate_pgs_id("pgs000001"),
                 "must match the format")
    expect_error(primedtools:::.validate_pgs_id("PGS0001"),
                 "must match the format")
    expect_error(primedtools:::.validate_pgs_id("PGS0000001"),
                 "must match the format")
    expect_error(primedtools:::.validate_pgs_id(""),
                 "non-empty character string")
    expect_error(primedtools:::.validate_pgs_id(123),
                 "non-empty character string")
    expect_error(primedtools:::.validate_pgs_id(character(0)),
                 "non-empty character string")
})

# Tests for .validate_genome_build -------------------------------------------

test_that(".validate_genome_build accepts valid builds", {
    expect_silent(primedtools:::.validate_genome_build("GRCh38"))
    expect_silent(primedtools:::.validate_genome_build("GRCh37"))
})

test_that(".validate_genome_build rejects invalid builds", {
    expect_error(primedtools:::.validate_genome_build("hg38"),
                 "must be one of")
    expect_error(primedtools:::.validate_genome_build("GRCh36"),
                 "must be one of")
    expect_error(primedtools:::.validate_genome_build("grch38"),
                 "must be one of")
})

# Tests for .get_workspace_attr -----------------------------------------------

test_that(".get_workspace_attr retrieves existing attribute", {
    data <- data.frame(
        key   = c("pgen", "psam", "pvar", "other"),
        value = c("gs://bucket/cohort.pgen",
                  "gs://bucket/cohort.psam",
                  "gs://bucket/cohort.pvar",
                  "something_else"),
        stringsAsFactors = FALSE
    )
    expect_equal(primedtools:::.get_workspace_attr(data, "pgen"),
                 "gs://bucket/cohort.pgen")
    expect_equal(primedtools:::.get_workspace_attr(data, "psam"),
                 "gs://bucket/cohort.psam")
    expect_equal(primedtools:::.get_workspace_attr(data, "pvar"),
                 "gs://bucket/cohort.pvar")
})

test_that(".get_workspace_attr errors when attribute missing", {
    data <- data.frame(
        key   = c("psam", "pvar"),
        value = c("gs://bucket/cohort.psam", "gs://bucket/cohort.pvar"),
        stringsAsFactors = FALSE
    )
    expect_error(primedtools:::.get_workspace_attr(data, "pgen"),
                 "workspace\\.pgen.*not found")
})

# Tests for .get_scorefile_path -----------------------------------------------

test_that(".get_scorefile_path returns correct path", {
    mock_tbl <- data.frame(
        pgs_scoring_file_id = c("a", "b"),
        pgs_model_id = c("PGS000001", "PGS000002"),
        file_path    = c("gs://bucket/PGS000001.txt.gz",
                         "gs://bucket/PGS000002.txt.gz"),
        file_type    = c("data", "data"),
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        avtable = function(...) mock_tbl,
        .package = "AnVILGCP",
        {
            result <- primedtools:::.get_scorefile_path(
                pgs_id = "PGS000001",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws"
            )
            expect_equal(result, "gs://bucket/PGS000001.txt.gz")
        }
    )
})

test_that(".get_scorefile_path errors when pgs_id not in table", {
    mock_tbl <- data.frame(
        pgs_scoring_file_id = "a",
        pgs_model_id = "PGS000002",
        file_path    = "gs://bucket/PGS000002.txt.gz",
        file_type    = "data",
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        avtable = function(...) mock_tbl,
        .package = "AnVILGCP",
        {
            expect_error(
                primedtools:::.get_scorefile_path(
                    pgs_id = "PGS000001",
                    workspace_namespace = "test-ns",
                    workspace_name = "test-ws"
                ),
                "No entry found"
            )
        }
    )
})

test_that(".get_scorefile_path errors when avtable call fails", {
    with_mocked_bindings(
        avtable = function(...) stop("table not found"),
        .package = "AnVILGCP",
        {
            expect_error(
                primedtools:::.get_scorefile_path(
                    pgs_id = "PGS000001",
                    workspace_namespace = "test-ns",
                    workspace_name = "test-ws"
                ),
                "Could not read 'pgs_scoring_file' table"
            )
        }
    )
})

# Tests for get_cohort_files --------------------------------------------------

test_that("get_cohort_files returns pgen/psam/pvar from workspace data", {
    mock_data <- data.frame(
        key   = c("pgen", "psam", "pvar"),
        value = c("gs://bucket/cohort.pgen",
                  "gs://bucket/cohort.psam",
                  "gs://bucket/cohort.pvar"),
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avdata = function(...) mock_data,
        .package = "AnVILGCP",
        {
            result <- get_cohort_files()
            expect_equal(result$pgen, "gs://bucket/cohort.pgen")
            expect_equal(result$psam, "gs://bucket/cohort.psam")
            expect_equal(result$pvar, "gs://bucket/cohort.pvar")
        }
    )
})

test_that("get_cohort_files errors when pvar is missing", {
    mock_data <- data.frame(
        key   = c("pgen", "psam"),
        value = c("gs://bucket/cohort.pgen", "gs://bucket/cohort.psam"),
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avdata = function(...) mock_data,
        .package = "AnVILGCP",
        {
            expect_error(get_cohort_files(), "workspace\\.pvar.*not found")
        }
    )
})

# Tests for run_pgs_pipeline input validation ---------------------------------

test_that("run_pgs_pipeline validates pgs_id format", {
    expect_error(
        run_pgs_pipeline(
            pgs_id = "bad_id",
            genome_build = "GRCh38",
            dest_bucket = "gs://bucket",
            sampleset_name = "cohort"
        ),
        "must match the format"
    )
})

test_that("run_pgs_pipeline validates genome_build", {
    expect_error(
        run_pgs_pipeline(
            pgs_id = "PGS000001",
            genome_build = "hg38",
            dest_bucket = "gs://bucket",
            sampleset_name = "cohort"
        ),
        "must be one of"
    )
})

test_that("run_pgs_pipeline errors when ancestry_adjust=TRUE but pcs=NULL", {
    expect_error(
        run_pgs_pipeline(
            pgs_id = "PGS000001",
            genome_build = "GRCh38",
            dest_bucket = "gs://bucket",
            sampleset_name = "cohort",
            ancestry_adjust = TRUE,
            pcs = NULL
        ),
        "'pcs' must be provided"
    )
})

# Tests for submit_fetch_pgs_workflow use_call_cache --------------------------

test_that("submit_fetch_pgs_workflow passes useCallCache=TRUE by default", {
    captured_args <- list()

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_configuration_get = function(...) list(),
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(config, ...) {
            captured_args <<- list(...)
            list(submissionId = "sub-fetch-default")
        },
        .package = "AnVILGCP",
        {
            result <- submit_fetch_pgs_workflow(
                pgs_id = "PGS000001",
                genome_build = "GRCh38",
                dest_bucket = "gs://bucket",
                model_url = "https://example.com/model.json",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws"
            )
        }
    )

    expect_equal(result, "sub-fetch-default")
    expect_true(captured_args$useCallCache)
})

test_that("submit_fetch_pgs_workflow passes useCallCache=FALSE when requested", {
    captured_args <- list()

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_configuration_get = function(...) list(),
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(config, ...) {
            captured_args <<- list(...)
            list(submissionId = "sub-fetch-nocache")
        },
        .package = "AnVILGCP",
        {
            result <- submit_fetch_pgs_workflow(
                pgs_id = "PGS000001",
                genome_build = "GRCh38",
                dest_bucket = "gs://bucket",
                model_url = "https://example.com/model.json",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws",
                use_call_cache = FALSE
            )
        }
    )

    expect_equal(result, "sub-fetch-nocache")
    expect_false(captured_args$useCallCache)
})

# Tests for submit_calc_pgs_workflow use_call_cache ---------------------------

test_that("submit_calc_pgs_workflow passes useCallCache=TRUE by default", {
    captured_args <- list()

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_configuration_get = function(...) list(),
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(config, ...) {
            captured_args <<- list(...)
            list(submissionId = "sub-calc-default")
        },
        .package = "AnVILGCP",
        {
            result <- submit_calc_pgs_workflow(
                pgs_model_id = "PGS000001",
                scorefile = "gs://bucket/score.txt.gz",
                genome_build = "GRCh38",
                pgen = "gs://bucket/cohort.pgen",
                psam = "gs://bucket/cohort.psam",
                pvar = "gs://bucket/cohort.pvar",
                min_overlap = 0.75,
                sampleset_name = "cohort",
                dest_bucket = "gs://bucket/results",
                model_url = "https://example.com/model.json",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws"
            )
        }
    )

    expect_equal(result, "sub-calc-default")
    expect_true(captured_args$useCallCache)
})

test_that("submit_calc_pgs_workflow passes useCallCache=FALSE when requested", {
    captured_args <- list()

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_configuration_get = function(...) list(),
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(config, ...) {
            captured_args <<- list(...)
            list(submissionId = "sub-calc-nocache")
        },
        .package = "AnVILGCP",
        {
            result <- submit_calc_pgs_workflow(
                pgs_model_id = "PGS000001",
                scorefile = "gs://bucket/score.txt.gz",
                genome_build = "GRCh38",
                pgen = "gs://bucket/cohort.pgen",
                psam = "gs://bucket/cohort.psam",
                pvar = "gs://bucket/cohort.pvar",
                min_overlap = 0.75,
                sampleset_name = "cohort",
                dest_bucket = "gs://bucket/results",
                model_url = "https://example.com/model.json",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws",
                use_call_cache = FALSE
            )
        }
    )

    expect_equal(result, "sub-calc-nocache")
    expect_false(captured_args$useCallCache)
})
