context("PGS pipeline validation and helper functions")

# Tests for .validate_pgs_id -------------------------------------------------

test_that(".validate_pgs_id accepts valid IDs", {
    expect_silent(primed_benchmarking:::.validate_pgs_id("PGS000001"))
    expect_silent(primed_benchmarking:::.validate_pgs_id("PGS999999"))
    expect_silent(primed_benchmarking:::.validate_pgs_id("PGS000123"))
})

test_that(".validate_pgs_id rejects invalid IDs", {
    expect_error(primed_benchmarking:::.validate_pgs_id("pgs000001"),
                 "must match the format")
    expect_error(primed_benchmarking:::.validate_pgs_id("PGS0001"),
                 "must match the format")
    expect_error(primed_benchmarking:::.validate_pgs_id("PGS0000001"),
                 "must match the format")
    expect_error(primed_benchmarking:::.validate_pgs_id(""),
                 "non-empty character string")
    expect_error(primed_benchmarking:::.validate_pgs_id(123),
                 "non-empty character string")
    expect_error(primed_benchmarking:::.validate_pgs_id(character(0)),
                 "non-empty character string")
})

# Tests for .validate_genome_build -------------------------------------------

test_that(".validate_genome_build accepts valid builds", {
    expect_silent(primed_benchmarking:::.validate_genome_build("GRCh38"))
    expect_silent(primed_benchmarking:::.validate_genome_build("GRCh37"))
})

test_that(".validate_genome_build rejects invalid builds", {
    expect_error(primed_benchmarking:::.validate_genome_build("hg38"),
                 "must be one of")
    expect_error(primed_benchmarking:::.validate_genome_build("GRCh36"),
                 "must be one of")
    expect_error(primed_benchmarking:::.validate_genome_build("grch38"),
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
    expect_equal(primed_benchmarking:::.get_workspace_attr(data, "pgen"),
                 "gs://bucket/cohort.pgen")
    expect_equal(primed_benchmarking:::.get_workspace_attr(data, "psam"),
                 "gs://bucket/cohort.psam")
    expect_equal(primed_benchmarking:::.get_workspace_attr(data, "pvar"),
                 "gs://bucket/cohort.pvar")
})

test_that(".get_workspace_attr errors when attribute missing", {
    data <- data.frame(
        key   = c("psam", "pvar"),
        value = c("gs://bucket/cohort.psam", "gs://bucket/cohort.pvar"),
        stringsAsFactors = FALSE
    )
    expect_error(primed_benchmarking:::.get_workspace_attr(data, "pgen"),
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
            result <- primed_benchmarking:::.get_scorefile_path(
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
                primed_benchmarking:::.get_scorefile_path(
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
                primed_benchmarking:::.get_scorefile_path(
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

# Tests for .find_successful_submission ---------------------------------------

test_that(".find_successful_submission returns NULL when no jobs exist", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = character(0),
                status         = character(0),
                succeeded      = integer(0),
                submissionRoot = character(0),
                stringsAsFactors = FALSE
            )
        },
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            expect_null(result)
        }
    )
})

test_that(".find_successful_submission returns NULL when no Done+succeeded jobs", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = c("sub-running", "sub-failed"),
                status         = c("Submitted", "Done"),
                succeeded      = c(0L, 0L),
                submissionRoot = c(
                    "gs://bucket/submissions/sub-running/primed_fetch_pgs_catalog/wf1",
                    "gs://bucket/submissions/sub-failed/primed_fetch_pgs_catalog/wf2"
                ),
                stringsAsFactors = FALSE
            )
        },
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            expect_null(result)
        }
    )
})

test_that(".find_successful_submission returns NULL when workflow name does not match", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = "sub-other",
                status         = "Done",
                succeeded      = 1L,
                submissionRoot = "gs://bucket/submissions/sub-other/primed_calc_pgs/wf1",
                stringsAsFactors = FALSE
            )
        },
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            expect_null(result)
        }
    )
})

test_that(".find_successful_submission returns NULL when submissionRoot path is unexpected", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = "sub-weird",
                status         = "Done",
                succeeded      = 1L,
                submissionRoot = "gs://bucket/other_prefix/sub-weird/primed_fetch_pgs_catalog/wf1",
                stringsAsFactors = FALSE
            )
        },
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            expect_null(result)
        }
    )
})

test_that(".find_successful_submission returns most recent matching submission", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = c("sub-newest", "sub-older"),
                status         = c("Done", "Done"),
                succeeded      = c(1L, 1L),
                submissionRoot = c(
                    "gs://bucket/submissions/sub-newest/primed_fetch_pgs_catalog/wf1",
                    "gs://bucket/submissions/sub-older/primed_fetch_pgs_catalog/wf2"
                ),
                stringsAsFactors = FALSE
            )
        },
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            # avworkflow_jobs returns sorted by submissionDate desc, so first = newest
            expect_equal(result, "sub-newest")
        }
    )
})

test_that(".find_successful_submission returns NULL when avworkflow_jobs errors", {
    with_mocked_bindings(
        avworkflow_jobs = function(...) stop("not in workspace"),
        .package = "AnVILGCP",
        {
            result <- primed_benchmarking:::.find_successful_submission(
                "primed_fetch_pgs_catalog",
                namespace = "test-ns",
                name = "test-ws"
            )
            expect_null(result)
        }
    )
})

# Tests for skip_if_complete in submit_fetch_pgs_workflow ---------------------

test_that("submit_fetch_pgs_workflow skips when prior success found", {
    avworkflow_run_called <- FALSE

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = "sub-prior",
                status         = "Done",
                succeeded      = 1L,
                submissionRoot = paste0(
                    "gs://bucket/submissions/sub-prior/",
                    "primed_fetch_pgs_catalog/wf1"
                ),
                stringsAsFactors = FALSE
            )
        },
        avworkflow_configuration_get    = function(...) { avworkflow_run_called <<- TRUE; list() },
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(...) { avworkflow_run_called <<- TRUE; list(submissionId = "new") },
        .package = "AnVILGCP",
        {
            result <- expect_message(
                submit_fetch_pgs_workflow(
                    pgs_id = "PGS000001",
                    genome_build = "GRCh38",
                    dest_bucket = "gs://bucket",
                    model_url = "https://example.com/model.json",
                    workspace_namespace = "test-ns",
                    workspace_name = "test-ws",
                    skip_if_complete = TRUE
                ),
                "prior successful primed_fetch_pgs_catalog"
            )
        }
    )

    expect_equal(result, "sub-prior")
    expect_false(avworkflow_run_called)
})

test_that("submit_fetch_pgs_workflow submits when skip_if_complete=TRUE but no prior", {
    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = character(0),
                status         = character(0),
                succeeded      = integer(0),
                submissionRoot = character(0),
                stringsAsFactors = FALSE
            )
        },
        avworkflow_configuration_get    = function(...) list(),
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(...) list(submissionId = "sub-new"),
        .package = "AnVILGCP",
        {
            result <- submit_fetch_pgs_workflow(
                pgs_id = "PGS000001",
                genome_build = "GRCh38",
                dest_bucket = "gs://bucket",
                model_url = "https://example.com/model.json",
                workspace_namespace = "test-ns",
                workspace_name = "test-ws",
                skip_if_complete = TRUE
            )
        }
    )

    expect_equal(result, "sub-new")
})

# Tests for skip_if_complete in submit_calc_pgs_workflow ----------------------

test_that("submit_calc_pgs_workflow skips when prior success found", {
    avworkflow_run_called <- FALSE

    with_mocked_bindings(
        avworkspace_namespace = function() "test-ns",
        avworkspace_name      = function() "test-ws",
        avworkflow_jobs = function(...) {
            data.frame(
                submissionId   = "sub-calc-prior",
                status         = "Done",
                succeeded      = 1L,
                submissionRoot = "gs://bucket/submissions/sub-calc-prior/primed_calc_pgs/wf1",
                stringsAsFactors = FALSE
            )
        },
        avworkflow_configuration_get    = function(...) { avworkflow_run_called <<- TRUE; list() },
        avworkflow_configuration_update = function(config, ...) config,
        avworkflow_run = function(...) { avworkflow_run_called <<- TRUE; list(submissionId = "new") },
        .package = "AnVILGCP",
        {
            result <- expect_message(
                submit_calc_pgs_workflow(
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
                    skip_if_complete = TRUE
                ),
                "prior successful primed_calc_pgs"
            )
        }
    )

    expect_equal(result, "sub-calc-prior")
    expect_false(avworkflow_run_called)
})
