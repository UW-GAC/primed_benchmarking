context("local ancestry workflow helpers")

# ---------------------------------------------------------------------------
# create_local_ancestry_summary_table
# ---------------------------------------------------------------------------

test_that("create_local_ancestry_summary_table errors on bad 'cohort'", {
    expect_error(
        create_local_ancestry_summary_table(
            cohort           = 123,
            cohort.namespace = "ns",
            cohort.name      = "ws"
        ),
        "non-NA scalar character"
    )
    expect_error(
        create_local_ancestry_summary_table(
            cohort           = NA_character_,
            cohort.namespace = "ns",
            cohort.name      = "ws"
        ),
        "non-NA scalar character"
    )
})

test_that("create_local_ancestry_summary_table errors on bad 'cohort.namespace'", {
    expect_error(
        create_local_ancestry_summary_table(
            cohort           = "C1",
            cohort.namespace = NULL,
            cohort.name      = "ws"
        ),
        "non-NA scalar character"
    )
})

test_that("create_local_ancestry_summary_table errors on bad 'cohort.name'", {
    expect_error(
        create_local_ancestry_summary_table(
            cohort           = "C1",
            cohort.namespace = "ns",
            cohort.name      = c("ws1", "ws2")
        ),
        "non-NA scalar character"
    )
})

test_that("create_local_ancestry_summary_table imports table when missing", {
    import_called <- FALSE

    with_mocked_bindings(
        avtables = function(namespace, name) {
            data.frame(table = character(0), stringsAsFactors = FALSE)
        },
        avtable_import = function(x, namespace, name) {
            import_called <<- TRUE
            invisible(NULL)
        },
        .package = "AnVILGCP",
        {
            result <- expect_message(
                create_local_ancestry_summary_table(
                    cohort           = "MY_COHORT",
                    cohort.namespace = "test-ns",
                    cohort.name      = "test-ws"
                ),
                "Creating table"
            )
        }
    )

    expect_identical(result, FALSE)
    expect_true(import_called)
})

test_that("create_local_ancestry_summary_table skips import when table exists", {
    import_called <- FALSE

    with_mocked_bindings(
        avtables = function(namespace, name) {
            data.frame(table = "local_ancestry_summary", stringsAsFactors = FALSE)
        },
        avtable_import = function(x, namespace, name) {
            import_called <<- TRUE
            invisible(NULL)
        },
        .package = "AnVILGCP",
        {
            result <- expect_message(
                create_local_ancestry_summary_table(
                    cohort           = "MY_COHORT",
                    cohort.namespace = "test-ns",
                    cohort.name      = "test-ws"
                ),
                "already exists"
            )
        }
    )

    expect_identical(result, TRUE)
    expect_false(import_called)
})

# ---------------------------------------------------------------------------
# set_up_step1c_summarize_local_ancestry_proportions
# ---------------------------------------------------------------------------

## Helper: build a minimal inputs/outputs data frame matching what
## avworkflow_configuration_inputs() / _outputs() would return.
.make_inputs <- function() {
    data.frame(
        name = c(
            "run_notebook_workflow.cohort",
            "run_notebook_workflow.data_file",
            "run_notebook_workflow.github_raw_url_notebook_file",
            "run_notebook_workflow.other_param"
        ),
        attribute = c("", "", "", "keep_me"),
        stringsAsFactors = FALSE
    )
}

.make_outputs <- function() {
    data.frame(
        name = c(
            "run_notebook_workflow.executed_notebook",
            "run_notebook_workflow.notebook_html",
            "run_notebook_workflow.summarize_local_ancestry_proportions",
            "run_notebook_workflow.two_way_proportion_IDs",
            "run_notebook_workflow.unrelated_output"
        ),
        attribute = c("", "", "", "", "unchanged"),
        stringsAsFactors = FALSE
    )
}

test_that("set_up_step1c mutates inputs correctly and returns config", {
    captured_set_args  <- list()
    captured_run_args  <- list()

    fake_config <- list(fake = TRUE)

    with_mocked_bindings(
        avworkflow_configuration_get = function(...) fake_config,
        avworkflow_configuration_inputs = function(config) .make_inputs(),
        avworkflow_configuration_outputs = function(config) .make_outputs(),
        avworkflow_configuration_update = function(config, inputs, outputs) {
            list(inputs = inputs, outputs = outputs)
        },
        avworkflow_configuration_set = function(x, namespace, name, dry) {
            captured_set_args <<- list(namespace = namespace, name = name, dry = dry)
            invisible(NULL)
        },
        avworkflow_run = function(config, entityName, entityType,
                                  deleteIntermediateOutputFiles,
                                  useCallCache, useReferenceDisks,
                                  namespace, name, dry) {
            captured_run_args <<- list(
                entityName = entityName, entityType = entityType,
                dry = dry, namespace = namespace, name = name
            )
            invisible(NULL)
        },
        .package = "AnVILGCP",
        {
            result <- set_up_step1c_summarize_local_ancestry_proportions(
                cohort                      = "MY_COHORT",
                cohort.namespace            = "test-ns",
                cohort.name                 = "test-ws",
                merged.6.ancestry_frac_path = "gs://bucket/merged.tsv"
            )
        }
    )

    ## inputs were mutated correctly
    inp <- result$inputs
    expect_equal(
        inp$attribute[inp$name == "run_notebook_workflow.cohort"],
        "\"MY_COHORT\""
    )
    expect_equal(
        inp$attribute[inp$name == "run_notebook_workflow.data_file"],
        "\"gs://bucket/merged.tsv\""
    )
    expect_true(grepl(
        "summarize_local_ancestry_proportions.ipynb",
        inp$attribute[inp$name == "run_notebook_workflow.github_raw_url_notebook_file"]
    ))
    ## untouched row preserved
    expect_equal(
        inp$attribute[inp$name == "run_notebook_workflow.other_param"],
        "keep_me"
    )

    ## outputs were mutated correctly
    out <- result$outputs
    expect_equal(
        out$attribute[out$name == "run_notebook_workflow.executed_notebook"],
        "this.notebook"
    )
    expect_equal(
        out$attribute[out$name == "run_notebook_workflow.notebook_html"],
        "this.notebook_html"
    )
    expect_equal(
        out$attribute[out$name == "run_notebook_workflow.summarize_local_ancestry_proportions"],
        "this.summarize_local_ancestry_proportions"
    )
    expect_equal(
        out$attribute[out$name == "run_notebook_workflow.two_way_proportion_IDs"],
        "this.two_way_proportion_IDs"
    )
    ## untouched output preserved
    expect_equal(
        out$attribute[out$name == "run_notebook_workflow.unrelated_output"],
        "unchanged"
    )

    ## dry-run set was called with dry = TRUE
    expect_true(captured_set_args$dry)
    expect_equal(captured_set_args$namespace, "test-ns")
    expect_equal(captured_set_args$name, "test-ws")

    ## dry-run submit was called with dry = TRUE
    expect_true(captured_run_args$dry)
    expect_equal(captured_run_args$entityName, "MY_COHORT")
    expect_equal(captured_run_args$entityType, "local_ancestry_summary")
})

# ---------------------------------------------------------------------------
# run_step1c_summarize_local_ancestry_proportions
# ---------------------------------------------------------------------------

test_that("run_step1c does nothing when run_now = FALSE", {
    set_called <- FALSE
    run_called <- FALSE

    with_mocked_bindings(
        avworkflow_configuration_set = function(...) { set_called <<- TRUE },
        avworkflow_run               = function(...) { run_called <<- TRUE },
        .package = "AnVILGCP",
        {
            result <- run_step1c_summarize_local_ancestry_proportions(
                cohort                      = "MY_COHORT",
                cohort.namespace            = "test-ns",
                cohort.name                 = "test-ws",
                merged.6.ancestry_frac_path = "gs://bucket/merged.tsv",
                run_now                     = FALSE,
                new_config                  = list()
            )
        }
    )

    expect_false(set_called)
    expect_false(run_called)
    expect_null(result)
})

test_that("run_step1c calls set and run when run_now = TRUE", {
    set_called <- FALSE
    run_called <- FALSE
    set_dry    <- NA
    run_dry    <- NA

    with_mocked_bindings(
        avworkflow_configuration_set = function(x, namespace, name, dry) {
            set_called <<- TRUE
            set_dry    <<- dry
            invisible(NULL)
        },
        avworkflow_run = function(config, entityName, entityType,
                                  deleteIntermediateOutputFiles,
                                  useCallCache, useReferenceDisks,
                                  namespace, name, dry) {
            run_called <<- TRUE
            run_dry    <<- dry
            invisible(NULL)
        },
        .package = "AnVILGCP",
        {
            run_step1c_summarize_local_ancestry_proportions(
                cohort                      = "MY_COHORT",
                cohort.namespace            = "test-ns",
                cohort.name                 = "test-ws",
                merged.6.ancestry_frac_path = "gs://bucket/merged.tsv",
                run_now                     = TRUE,
                new_config                  = list()
            )
        }
    )

    expect_true(set_called)
    expect_true(run_called)
    expect_false(set_dry)
    expect_false(run_dry)
})
