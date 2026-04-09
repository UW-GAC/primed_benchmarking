#' Create local ancestry summary table
#'
#' Checks whether the workspace data table \code{local_ancestry_summary} exists
#' in the specified AnVIL workspace. If the table is absent, a minimal table
#' with a single \code{local_ancestry_summary_id} column (populated by
#' \code{cohort}) is created and imported. If it already exists, a message is
#' emitted and nothing is changed.
#'
#' @param cohort Character scalar. The cohort/entity ID to use as the initial
#'   \code{local_ancestry_summary_id} value.
#' @param cohort.namespace Character scalar. AnVIL workspace namespace.
#' @param cohort.name Character scalar. AnVIL workspace name.
#'
#' @return Invisibly returns \code{FALSE} if the table was created, or
#'   \code{TRUE} if it already existed.
#'
#' @examples
#' \dontrun{
#' create_local_ancestry_summary_table(
#'   cohort           = "MY_COHORT",
#'   cohort.namespace = "my-namespace",
#'   cohort.name      = "my-workspace"
#' )
#' }
#'
#' @importFrom AnVILGCP avtables avtable_import
#' @importFrom dplyr filter
#' @importFrom tibble as_tibble
#' @export
create_local_ancestry_summary_table <- function(cohort, cohort.namespace, cohort.name) {
    if (!is.character(cohort) || length(cohort) != 1L || is.na(cohort)) {
        stop("'cohort' must be a non-NA scalar character string.")
    }
    if (!is.character(cohort.namespace) || length(cohort.namespace) != 1L ||
            is.na(cohort.namespace)) {
        stop("'cohort.namespace' must be a non-NA scalar character string.")
    }
    if (!is.character(cohort.name) || length(cohort.name) != 1L ||
            is.na(cohort.name)) {
        stop("'cohort.name' must be a non-NA scalar character string.")
    }

    info <- AnVILGCP::avtables(namespace = cohort.namespace, name = cohort.name)
    local_info <- dplyr::filter(info, .data$table == "local_ancestry_summary")

    if (nrow(local_info) == 0L) {
        message("Creating table in workspace: ", cohort.namespace, "/", cohort.name)
        summary_tbl <- tibble::as_tibble(
            data.frame(local_ancestry_summary_id = cohort,
                       stringsAsFactors = FALSE)
        )
        AnVILGCP::avtable_import(summary_tbl,
                                 namespace = cohort.namespace,
                                 name      = cohort.name)
        return(invisible(FALSE))
    }

    message("Table already exists in workspace")
    invisible(TRUE)
}


#' Set up Step 1C: summarize local ancestry proportions (dry run)
#'
#' Retrieves the \code{step1c_summarize_local_ancestry_proportions} workflow
#' configuration from the workspace, updates its inputs and outputs for the
#' specified cohort and data file, then performs a dry-run validation
#' (\code{\link[AnVILGCP]{avworkflow_configuration_set}} with \code{dry = TRUE})
#' and a dry-run submission
#' (\code{\link[AnVILGCP]{avworkflow_run}} with \code{dry = TRUE}).
#' Use \code{\link{run_step1c_summarize_local_ancestry_proportions}} to submit
#' the actual job.
#'
#' @param cohort Character scalar. The entity name to run on.
#' @param cohort.namespace Character scalar. AnVIL workspace namespace.
#' @param cohort.name Character scalar. AnVIL workspace name.
#' @param merged.6.ancestry_frac_path Character scalar. GCS path to the merged
#'   6-ancestry fraction file produced by an upstream workflow step.
#'
#' @return The updated workflow configuration object (as returned by
#'   \code{\link[AnVILGCP]{avworkflow_configuration_update}}).
#'
#' @examples
#' \dontrun{
#' new_config <- set_up_step1c_summarize_local_ancestry_proportions(
#'   cohort                    = "MY_COHORT",
#'   cohort.namespace          = "my-namespace",
#'   cohort.name               = "my-workspace",
#'   merged.6.ancestry_frac_path = "gs://my-bucket/merged_6anc_frac.tsv"
#' )
#' }
#'
#' @importFrom AnVILGCP avworkflow_configuration_get
#' @importFrom AnVILGCP avworkflow_configuration_inputs
#' @importFrom AnVILGCP avworkflow_configuration_outputs
#' @importFrom AnVILGCP avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_configuration_set
#' @importFrom AnVILGCP avworkflow_run
#' @importFrom dplyr mutate
#' @export
set_up_step1c_summarize_local_ancestry_proportions <- function(
    cohort,
    cohort.namespace,
    cohort.name,
    merged.6.ancestry_frac_path
) {
    config <- AnVILGCP::avworkflow_configuration_get(
        workflow_namespace = "primed_benchmarking",
        workflow_name      = "step1c_summarize_local_ancestry_proportions",
        namespace          = cohort.namespace,
        name               = cohort.name
    )

    ## --- inputs -------------------------------------------------------
    inputs <- AnVILGCP::avworkflow_configuration_inputs(config)

    inputs <- dplyr::mutate(
        inputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.cohort",
            paste0("\"", cohort, "\""),
            .data$attribute
        )
    )
    inputs <- dplyr::mutate(
        inputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.data_file",
            paste0("\"", merged.6.ancestry_frac_path, "\""),
            .data$attribute
        )
    )
    inputs <- dplyr::mutate(
        inputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.github_raw_url_notebook_file",
            paste0("\"https://raw.githubusercontent.com/UW-GAC/primed_benchmarking/",
                   "refs/heads/initial_workspace_setup_workflow_configuration/",
                   "notebook_templates/",
                   "summarize_local_ancestry_proportions.ipynb\""),
            .data$attribute
        )
    )

    ## --- outputs ------------------------------------------------------
    outputs <- AnVILGCP::avworkflow_configuration_outputs(config)

    outputs <- dplyr::mutate(
        outputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.executed_notebook",
            "this.notebook",
            .data$attribute
        )
    )
    outputs <- dplyr::mutate(
        outputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.notebook_html",
            "this.notebook_html",
            .data$attribute
        )
    )
    outputs <- dplyr::mutate(
        outputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.summarize_local_ancestry_proportions",
            "this.summarize_local_ancestry_proportions",
            .data$attribute
        )
    )
    outputs <- dplyr::mutate(
        outputs,
        attribute = ifelse(
            .data$name == "run_notebook_workflow.two_way_proportion_IDs",
            "this.two_way_proportion_IDs",
            .data$attribute
        )
    )

    ## --- update configuration -----------------------------------------
    new_config <- AnVILGCP::avworkflow_configuration_update(
        config  = config,
        inputs  = inputs,
        outputs = outputs
    )

    ## --- dry-run validate ---------------------------------------------
    AnVILGCP::avworkflow_configuration_set(
        new_config,
        namespace = cohort.namespace,
        name      = cohort.name,
        dry       = TRUE
    )

    ## --- dry-run submit -----------------------------------------------
    AnVILGCP::avworkflow_run(
        config                        = new_config,
        entityName                    = cohort,
        entityType                    = "local_ancestry_summary",
        deleteIntermediateOutputFiles = FALSE,
        useCallCache                  = TRUE,
        useReferenceDisks             = FALSE,
        namespace                     = cohort.namespace,
        name                          = cohort.name,
        dry                           = TRUE
    )

    new_config
}


#' Run Step 1C: summarize local ancestry proportions
#'
#' Applies the workflow configuration and submits the
#' \code{step1c_summarize_local_ancestry_proportions} workflow for real (i.e.,
#' not a dry run) when \code{run_now = TRUE}. When \code{run_now = FALSE}
#' (default) the function returns invisibly without doing anything, making it
#' safe to call in a script that is still being prepared.
#'
#' @param cohort Character scalar. The entity name to run on.
#' @param cohort.namespace Character scalar. AnVIL workspace namespace.
#' @param cohort.name Character scalar. AnVIL workspace name.
#' @param merged.6.ancestry_frac_path Character scalar. (Retained for API
#'   symmetry with
#'   \code{\link{set_up_step1c_summarize_local_ancestry_proportions}}; not used
#'   directly by this function.)
#' @param run_now Logical scalar. If \code{TRUE}, apply the configuration and
#'   submit the workflow. Default is \code{FALSE}.
#' @param new_config A workflow configuration object previously returned by
#'   \code{\link{set_up_step1c_summarize_local_ancestry_proportions}}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' new_config <- set_up_step1c_summarize_local_ancestry_proportions(
#'   cohort                    = "MY_COHORT",
#'   cohort.namespace          = "my-namespace",
#'   cohort.name               = "my-workspace",
#'   merged.6.ancestry_frac_path = "gs://my-bucket/merged_6anc_frac.tsv"
#' )
#' run_step1c_summarize_local_ancestry_proportions(
#'   cohort                    = "MY_COHORT",
#'   cohort.namespace          = "my-namespace",
#'   cohort.name               = "my-workspace",
#'   merged.6.ancestry_frac_path = "gs://my-bucket/merged_6anc_frac.tsv",
#'   run_now                   = TRUE,
#'   new_config                = new_config
#' )
#' }
#'
#' @importFrom AnVILGCP avworkflow_configuration_set avworkflow_run
#' @export
run_step1c_summarize_local_ancestry_proportions <- function(
    cohort,
    cohort.namespace,
    cohort.name,
    merged.6.ancestry_frac_path,
    run_now    = FALSE,
    new_config
) {
    if (isTRUE(run_now)) {
        AnVILGCP::avworkflow_configuration_set(
            new_config,
            namespace = cohort.namespace,
            name      = cohort.name,
            dry       = FALSE
        )

        AnVILGCP::avworkflow_run(
            config                        = new_config,
            entityName                    = cohort,
            entityType                    = "local_ancestry_summary",
            deleteIntermediateOutputFiles = FALSE,
            useCallCache                  = TRUE,
            useReferenceDisks             = FALSE,
            namespace                     = cohort.namespace,
            name                          = cohort.name,
            dry                           = FALSE
        )
    }

    invisible(NULL)
}
