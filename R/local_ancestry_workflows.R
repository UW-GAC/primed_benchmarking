#' Create Local Ancestry Summary Table
#'
#' This function checks for the workspace table "local_ancestry_summary" via AnVILGCP::avtables() and imports it if missing.
#'
#' @param cohort A character scalar for cohort ID.
#' @param cohort.namespace A character scalar for the cohort namespace.
#' @param cohort.name A character scalar for the cohort name.
#' @export
create_local_ancestry_summary_table <- function(cohort, cohort.namespace, cohort.name) {
    if (!AnVILGCP::avtables(namespace = cohort.namespace, name = cohort.name)) {
        message('Importing local ancestry summary table...')
        AnVILGCP::avtable_import(local_ancestry_summary_id = cohort)
    }
}

#' Set Up Step 1C: Summarize Local Ancestry Proportions
#'
#' This function modifies workflow configuration and performs a dry-run submission.
#'
#' @param cohort A character scalar for cohort ID.
#' @param cohort.namespace A character scalar for the cohort namespace.
#' @param cohort.name A character scalar for the cohort name.
#' @param merged.6.ancestry_frac_path A character scalar for the ancestry fraction path.
#' @export
set_up_step1c_summarize_local_ancestry_proportions <- function(cohort, cohort.namespace, cohort.name, merged.6.ancestry_frac_path) {
    config <- AnVILGCP::avworkflow_configuration_get(workflow_namespace = 'primed_benchmarking', workflow_name = 'step1c_summarize_local_ancestry_proportions', namespace = cohort.namespace, name = cohort.name)
    inputs <- AnVILGCP::avworkflow_configuration_inputs(config)
    inputs[['run_notebook_workflow.cohort']] <- cohort
    inputs[['run_notebook_workflow.data_file']] <- merged.6.ancestry_frac_path
    inputs[['run_notebook_workflow.github_raw_url_notebook_file']] <- 'https://raw.githubusercontent.com/UW-GAC/primed_benchmarking/refs/heads/initial_workspace_setup_workflow_configuration/notebook_templates/summarize_local_ancestry_proportions.ipynb'
    outputs <- AnVILGCP::avworkflow_configuration_outputs(config)
    outputs[['executed_notebook']] <- 'this.notebook'
    outputs[['notebook_html']] <- 'this.notebook.html'
    outputs[['summarize_local_ancestry_proportions']] <- 'this.proportions'
    outputs[['two_way_proportion_IDs']] <- 'this.id'
    AnVILGCP::avworkflow_configuration_update(config)
    AnVILGCP::avworkflow_configuration_set(dry = TRUE)
    new_config <- AnVILGCP::avworkflow_run(dry = TRUE, entityName = cohort, entityType = 'local_ancestry_summary', deleteIntermediateOutputFiles = FALSE, useCallCache = TRUE, useReferenceDisks = FALSE, namespace = cohort.namespace, name = cohort.name)
    return(new_config)
}

#' Run Step 1C: Summarize Local Ancestry Proportions
#'
#' This function runs the workflow immediately if run_now is TRUE.
#'
#' @param cohort A character scalar for cohort ID.
#' @param cohort.namespace A character scalar for the cohort namespace.
#' @param cohort.name A character scalar for the cohort name.
#' @param merged.6.ancestry_frac_path A character scalar for the ancestry fraction path.
#' @param run_now A logical indicating whether to run immediately.
#' @param new_config The new configuration to use.
#' @export
run_step1c_summarize_local_ancestry_proportions <- function(cohort, cohort.namespace, cohort.name, merged.6.ancestry_frac_path, run_now = FALSE, new_config) {
    if (run_now) {
        AnVILGCP::avworkflow_configuration_set(new_config, dry = FALSE, namespace = cohort.namespace, name = cohort.name)
        AnVILGCP::avworkflow_run(config = new_config, entityName = cohort, entityType = 'local_ancestry_summary', deleteIntermediateOutputFiles = FALSE, useCallCache = TRUE, useReferenceDisks = FALSE, namespace = cohort.namespace, name = cohort.name, dry = FALSE)
    }
}

#' @importFrom AnVILGCP avtables avtable_import avworkflow_configuration_get avworkflow_configuration_inputs avworkflow_configuration_outputs avworkflow_configuration_update avworkflow_configuration_set avworkflow_run
#' @importFrom dplyr filter mutate if_else
#' @importFrom tibble as_tibble
