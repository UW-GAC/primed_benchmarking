#' Run the PRIMED PGS pipeline for a given PGS Catalog score
#'
#' This function orchestrates the full PRIMED PGS pipeline for a given PGS
#' Catalog ID. It reads cohort genotype file paths from the workspace data
#' attributes (\code{workspace.pgen}, \code{workspace.psam},
#' \code{workspace.pvar}), then submits the \code{primed_fetch_pgs_catalog}
#' workflow to fetch the scoring file and import it into the workspace data
#' model, waits for that workflow to complete, retrieves the scoring file path
#' from the imported \code{pgs_scoring_file} data table, and finally submits
#' the \code{primed_calc_pgs} workflow to calculate individual-level scores
#' on the cohort.
#'
#' @name run_pgs_pipeline
#' @param pgs_id Character. ID of the score in the PGS Catalog,
#'   e.g. \code{"PGS000001"}.
#' @param genome_build Character. Genome build to use; either \code{"GRCh38"}
#'   or \code{"GRCh37"}.
#' @param dest_bucket Character. Google Cloud Storage bucket path (starting
#'   with \code{"gs://"}) where scoring and individual score files will be
#'   written.
#' @param sampleset_name Character. A name for the cohort sample set used to
#'   construct output file names.
#' @param model_url Character. URL to the PRIMED PGS data model JSON file.
#'   Defaults to the main branch of the PRIMED data models repository.
#' @param min_overlap Numeric. Minimum fraction of score variants that must
#'   be present in the genotype data (e.g. \code{0.75} for 75\% overlap).
#'   Default is \code{0.75}.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configurations in the workspace. Defaults to \code{workspace_namespace}.
#' @param overwrite Logical. Whether to overwrite existing rows in the data
#'   tables. Default is \code{FALSE}.
#' @param ancestry_adjust Logical. Whether to perform ancestry-based score
#'   adjustment. If \code{TRUE}, the \code{pcs} parameter must be provided.
#'   Default is \code{FALSE}.
#' @param pcs Character or \code{NULL}. Google Cloud Storage path to a file
#'   containing principal components for ancestry adjustment. Required when
#'   \code{ancestry_adjust = TRUE}.
#' @param primed_dataset_id Character or \code{NULL}. Optional PRIMED dataset
#'   identifier to include in the output data tables.
#' @param poll_interval Numeric. Number of seconds to wait between polling
#'   the workflow status. Default is \code{60}.
#' @param timeout Numeric. Maximum number of seconds to wait for the
#'   \code{primed_fetch_pgs_catalog} workflow to complete before timing out.
#'   Default is \code{3600} (1 hour).
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   workflow submissions. Default is \code{TRUE}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{fetch_submission}}{Submission ID of the
#'       \code{primed_fetch_pgs_catalog} workflow run.}
#'     \item{\code{calc_submission}}{Submission ID of the
#'       \code{primed_calc_pgs} workflow run.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- run_pgs_pipeline(
#'   pgs_id = "PGS000001",
#'   genome_build = "GRCh38",
#'   dest_bucket = "gs://my-bucket/pgs_results",
#'   sampleset_name = "my_cohort"
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name avtable
#' @export
run_pgs_pipeline <- function(
    pgs_id,
    genome_build,
    dest_bucket,
    sampleset_name,
    model_url = paste0(
        "https://raw.githubusercontent.com/UW-GAC/primed_data_models/",
        "refs/heads/main/PRIMED_PGS_data_model.json"
    ),
    min_overlap = 0.75,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    overwrite = FALSE,
    ancestry_adjust = FALSE,
    pcs = NULL,
    primed_dataset_id = NULL,
    poll_interval = 60,
    timeout = 3600,
    use_call_cache = TRUE
) {
    .validate_pgs_id(pgs_id)
    .validate_genome_build(genome_build)

    if (ancestry_adjust && is.null(pcs)) {
        stop("'pcs' must be provided when 'ancestry_adjust = TRUE'")
    }

    # Read cohort genotype files from workspace data attributes
    cohort <- get_cohort_files(
        workspace_namespace = workspace_namespace,
        workspace_name = workspace_name
    )

    # Submit the primed_fetch_pgs_catalog workflow
    message("Submitting primed_fetch_pgs_catalog workflow for ", pgs_id, "...")
    fetch_submission <- submit_fetch_pgs_workflow(
        pgs_id = pgs_id,
        genome_build = genome_build,
        dest_bucket = dest_bucket,
        model_url = model_url,
        workspace_namespace = workspace_namespace,
        workspace_name = workspace_name,
        workflow_namespace = workflow_namespace,
        overwrite = overwrite,
        use_call_cache = use_call_cache
    )
    message("Submitted primed_fetch_pgs_catalog workflow: ", fetch_submission)

    # Wait for the fetch workflow to complete
    message("Waiting for primed_fetch_pgs_catalog workflow to complete...")
    wait_for_workflow(
        submission_id = fetch_submission,
        workspace_namespace = workspace_namespace,
        workspace_name = workspace_name,
        poll_interval = poll_interval,
        timeout = timeout
    )
    message("primed_fetch_pgs_catalog workflow completed.")

    # Retrieve the scoring file path from the pgs_scoring_file data table
    scorefile <- .get_scorefile_path(
        pgs_id = pgs_id,
        workspace_namespace = workspace_namespace,
        workspace_name = workspace_name
    )
    message("Using scoring file: ", scorefile)

    # Submit the primed_calc_pgs workflow
    message("Submitting primed_calc_pgs workflow for ", pgs_id, "...")
    calc_submission <- submit_calc_pgs_workflow(
        pgs_model_id = pgs_id,
        scorefile = scorefile,
        genome_build = genome_build,
        pgen = cohort$pgen,
        psam = cohort$psam,
        pvar = cohort$pvar,
        min_overlap = min_overlap,
        sampleset_name = sampleset_name,
        dest_bucket = dest_bucket,
        model_url = model_url,
        workspace_namespace = workspace_namespace,
        workspace_name = workspace_name,
        workflow_namespace = workflow_namespace,
        overwrite = overwrite,
        ancestry_adjust = ancestry_adjust,
        pcs = pcs,
        primed_dataset_id = primed_dataset_id,
        use_call_cache = use_call_cache
    )
    message("Submitted primed_calc_pgs workflow: ", calc_submission)

    list(
        fetch_submission = fetch_submission,
        calc_submission = calc_submission
    )
}


#' Read cohort genotype file paths from workspace data attributes
#'
#' Retrieves the paths to the cohort genotype files stored as workspace-level
#' data attributes named \code{pgen}, \code{psam}, and \code{pvar}.
#' In the AnVIL workspace, these are exposed as \code{workspace.pgen},
#' \code{workspace.psam}, and \code{workspace.pvar}.
#'
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#'
#' @return A named list with three character elements:
#'   \describe{
#'     \item{\code{pgen}}{Google Cloud Storage path to the \code{.pgen} file.}
#'     \item{\code{psam}}{Google Cloud Storage path to the \code{.psam} file.}
#'     \item{\code{pvar}}{Google Cloud Storage path to the \code{.pvar} file.}
#'   }
#'
#' @examples
#' \dontrun{
#' cohort <- get_cohort_files()
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name avdata
#' @export
get_cohort_files <- function(
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name()
) {
    data <- avdata(namespace = workspace_namespace, name = workspace_name)

    pgen <- .get_workspace_attr(data, "pgen")
    psam <- .get_workspace_attr(data, "psam")
    pvar <- .get_workspace_attr(data, "pvar")

    list(pgen = pgen, psam = psam, pvar = pvar)
}


#' Submit the primed_fetch_pgs_catalog workflow
#'
#' Configures and submits the \code{primed_fetch_pgs_catalog} workflow in the
#' current AnVIL workspace. This workflow fetches a scoring file from the PGS
#' Catalog, copies it to the specified Google Cloud bucket, and imports
#' metadata into the workspace data tables in the PRIMED PGS data model.
#'
#' The workflow method configuration named \code{primed_fetch_pgs_catalog}
#' must already be imported into the workspace before calling this function.
#'
#' @param pgs_id Character. ID of the score in the PGS Catalog,
#'   e.g. \code{"PGS000001"}.
#' @param genome_build Character. Genome build to use; either \code{"GRCh38"}
#'   or \code{"GRCh37"}.
#' @param dest_bucket Character. Google Cloud Storage bucket path where scoring
#'   files will be written.
#' @param model_url Character. URL to the PRIMED PGS data model JSON file.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configuration. Defaults to \code{workspace_namespace}.
#' @param overwrite Logical. Whether to overwrite existing rows in data tables.
#'   Default is \code{FALSE}.
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   the workflow submission. Default is \code{TRUE}.
#'
#' @return Character. The submission ID of the workflow run.
#'
#' @examples
#' \dontrun{
#' submission_id <- submit_fetch_pgs_workflow(
#'   pgs_id = "PGS000001",
#'   genome_build = "GRCh38",
#'   dest_bucket = "gs://my-bucket/pgs_catalog",
#'   model_url = paste0(
#'     "https://raw.githubusercontent.com/UW-GAC/primed_data_models/",
#'     "refs/heads/main/PRIMED_PGS_data_model.json"
#'   )
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name
#' @importFrom AnVILGCP avworkflow_configuration_get avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_run
#' @export
submit_fetch_pgs_workflow <- function(
    pgs_id,
    genome_build,
    dest_bucket,
    model_url,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    overwrite = FALSE,
    use_call_cache = TRUE
) {
        namespace = workspace_namespace,
        name = workspace_name,
        workflow_namespace = workflow_namespace,
        workflow_name = "primed_fetch_pgs_catalog"
    )

    inputs <- list(
        "primed_fetch_pgs_catalog.pgs_id"           = as.list(pgs_id),
        "primed_fetch_pgs_catalog.genome_build"     = genome_build,
        "primed_fetch_pgs_catalog.dest_bucket"      = dest_bucket,
        "primed_fetch_pgs_catalog.model_url"        = model_url,
        "primed_fetch_pgs_catalog.workspace_name"   = workspace_name,
        "primed_fetch_pgs_catalog.workspace_namespace" = workspace_namespace,
        "primed_fetch_pgs_catalog.overwrite"        = overwrite,
        "primed_fetch_pgs_catalog.import_tables"    = TRUE
    )

    config <- avworkflow_configuration_update(config, inputs = inputs)
    result <- avworkflow_run(config,
                             namespace = workspace_namespace,
                             name = workspace_name,
                             submit = TRUE,
                             useCallCache = use_call_cache)
    result$submissionId
}


#' Submit the primed_calc_pgs workflow
#'
#' Configures and submits the \code{primed_calc_pgs} workflow in the current
#' AnVIL workspace. This workflow matches the provided scoring file to the
#' cohort genotype data, calculates individual-level polygenic scores using
#' PLINK2, optionally adjusts for ancestry using PCs, and imports results
#' into the workspace data tables.
#'
#' The workflow method configuration named \code{primed_calc_pgs} must
#' already be imported into the workspace before calling this function.
#'
#' @param pgs_model_id Character. Identifier for the PGS model, typically the
#'   PGS Catalog ID, e.g. \code{"PGS000001"}.
#' @param scorefile Character. Google Cloud Storage path to the scoring file
#'   fetched from the PGS Catalog.
#' @param genome_build Character. Genome build of the scoring file and genotype
#'   data; either \code{"GRCh38"} or \code{"GRCh37"}.
#' @param pgen Character. Google Cloud Storage path to the cohort \code{.pgen}
#'   file.
#' @param psam Character. Google Cloud Storage path to the cohort \code{.psam}
#'   file.
#' @param pvar Character. Google Cloud Storage path to the cohort \code{.pvar}
#'   file.
#' @param min_overlap Numeric. Minimum fraction of score variants that must be
#'   present in the genotype data.
#' @param sampleset_name Character. Name used to construct output file names.
#' @param dest_bucket Character. Google Cloud Storage bucket path where output
#'   score files will be written.
#' @param model_url Character. URL to the PRIMED PGS data model JSON file.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configuration. Defaults to \code{workspace_namespace}.
#' @param overwrite Logical. Whether to overwrite existing rows in data tables.
#'   Default is \code{FALSE}.
#' @param ancestry_adjust Logical. Whether to adjust scores for ancestry.
#'   Default is \code{FALSE}.
#' @param pcs Character or \code{NULL}. Google Cloud Storage path to a file
#'   with principal components for ancestry adjustment.
#' @param primed_dataset_id Character or \code{NULL}. Optional PRIMED dataset
#'   identifier.
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   the workflow submission. Default is \code{TRUE}.
#'
#' @return Character. The submission ID of the workflow run.
#'
#' @examples
#' \dontrun{
#' submission_id <- submit_calc_pgs_workflow(
#'   pgs_model_id = "PGS000001",
#'   scorefile = "gs://my-bucket/pgs_catalog/PGS000001_hmPOS_GRCh38.txt.gz",
#'   genome_build = "GRCh38",
#'   pgen = "gs://my-bucket/cohort/cohort.pgen",
#'   psam = "gs://my-bucket/cohort/cohort.psam",
#'   pvar = "gs://my-bucket/cohort/cohort.pvar",
#'   min_overlap = 0.75,
#'   sampleset_name = "my_cohort",
#'   dest_bucket = "gs://my-bucket/pgs_results",
#'   model_url = paste0(
#'     "https://raw.githubusercontent.com/UW-GAC/primed_data_models/",
#'     "refs/heads/main/PRIMED_PGS_data_model.json"
#'   )
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name
#' @importFrom AnVILGCP avworkflow_configuration_get avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_run
#' @export
submit_calc_pgs_workflow <- function(
    pgs_model_id,
    scorefile,
    genome_build,
    pgen,
    psam,
    pvar,
    min_overlap,
    sampleset_name,
    dest_bucket,
    model_url,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    overwrite = FALSE,
    ancestry_adjust = FALSE,
    pcs = NULL,
    primed_dataset_id = NULL,
    use_call_cache = TRUE
) {
    .validate_genome_build(genome_build)

    config <- avworkflow_configuration_get(
        config = "primed_calc_pgs",
        namespace = workspace_namespace,
        name = workspace_name,
        workflow_namespace = workflow_namespace,
        workflow_name = "primed_calc_pgs"
    )

    inputs <- list(
        "primed_calc_pgs.pgs_model_id"          = pgs_model_id,
        "primed_calc_pgs.scorefile"             = scorefile,
        "primed_calc_pgs.genome_build"          = genome_build,
        "primed_calc_pgs.pgen"                  = pgen,
        "primed_calc_pgs.psam"                  = psam,
        "primed_calc_pgs.pvar"                  = pvar,
        "primed_calc_pgs.min_overlap"           = min_overlap,
        "primed_calc_pgs.sampleset_name"        = sampleset_name,
        "primed_calc_pgs.dest_bucket"           = dest_bucket,
        "primed_calc_pgs.model_url"             = model_url,
        "primed_calc_pgs.workspace_name"        = workspace_name,
        "primed_calc_pgs.workspace_namespace"   = workspace_namespace,
        "primed_calc_pgs.overwrite"             = overwrite,
        "primed_calc_pgs.import_tables"         = TRUE,
        "primed_calc_pgs.ancestry_adjust"       = ancestry_adjust
    )

    if (!is.null(pcs)) {
        inputs[["primed_calc_pgs.pcs"]] <- pcs
    }

    if (!is.null(primed_dataset_id)) {
        inputs[["primed_calc_pgs.primed_dataset_id"]] <- primed_dataset_id
    }

    config <- avworkflow_configuration_update(config, inputs = inputs)
    result <- avworkflow_run(config,
                             namespace = workspace_namespace,
                             name = workspace_name,
                             submit = TRUE,
                             useCallCache = use_call_cache)
    result$submissionId
}


#' Wait for an AnVIL workflow submission to complete
#'
#' Polls the status of an AnVIL workflow submission at regular intervals until
#' all workflows in the submission have either succeeded or failed.
#'
#' @param submission_id Character. The submission ID returned by
#'   \code{\link{submit_fetch_pgs_workflow}} or
#'   \code{\link{submit_calc_pgs_workflow}}.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param poll_interval Numeric. Number of seconds to wait between status
#'   checks. Default is \code{60}.
#' @param timeout Numeric. Maximum number of seconds to wait before stopping
#'   with a timeout error. Default is \code{3600} (1 hour).
#'
#' @return Invisibly returns the final job status tibble.
#'
#' @examples
#' \dontrun{
#' wait_for_workflow("abc123-submission-id")
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name avworkflow_jobs
#' @export
wait_for_workflow <- function(
    submission_id,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    poll_interval = 60,
    timeout = 3600
) {
    terminal_states <- c("Succeeded", "Failed", "Aborted")
    start_time <- proc.time()[["elapsed"]]

    repeat {
        jobs <- avworkflow_jobs(
            submissionId = submission_id,
            namespace = workspace_namespace,
            name = workspace_name
        )

        statuses <- jobs$status
        message("  Workflow statuses: ",
                paste(table(statuses), names(table(statuses)), sep = " ", collapse = ", "))

        if (all(statuses %in% terminal_states)) {
            failed <- jobs[jobs$status %in% c("Failed", "Aborted"), ]
            if (nrow(failed) > 0) {
                stop("Workflow(s) in submission ", submission_id,
                     " did not succeed. Failed/Aborted workflows:\n",
                     paste(failed$workflowId, failed$status, sep = ": ",
                           collapse = "\n"))
            }
            return(invisible(jobs))
        }

        elapsed <- proc.time()[["elapsed"]] - start_time
        if (elapsed > timeout) {
            stop("Timeout: submission ", submission_id,
                 " did not complete within ", timeout, " seconds")
        }

        Sys.sleep(poll_interval)
    }
}


# ---- Internal helpers --------------------------------------------------------

#' @keywords internal
.validate_pgs_id <- function(pgs_id) {
    if (!is.character(pgs_id) || length(pgs_id) == 0) {
        stop("'pgs_id' must be a non-empty character string")
    }
    if (!all(grepl("^PGS[0-9]{6}$", pgs_id))) {
        stop("'pgs_id' must match the format 'PGS' followed by six digits ",
             "(e.g. 'PGS000001')")
    }
}


#' @keywords internal
.validate_genome_build <- function(genome_build) {
    valid <- c("GRCh38", "GRCh37")
    if (!genome_build %in% valid) {
        stop("'genome_build' must be one of: ",
             paste(valid, collapse = ", "))
    }
}


#' @keywords internal
.get_workspace_attr <- function(data, key) {
    val <- data$value[data$key == key]
    if (length(val) == 0 || is.na(val)) {
        stop("Workspace attribute '", key, "' (workspace.", key,
             ") not found in workspace data. ",
             "Please set this attribute in the workspace Data tab.")
    }
    val
}


#' @keywords internal
.get_scorefile_path <- function(pgs_id, workspace_namespace, workspace_name) {
    tbl <- tryCatch(
        avtable("pgs_scoring_file",
                namespace = workspace_namespace,
                name = workspace_name),
        error = function(e) {
            stop("Could not read 'pgs_scoring_file' table from workspace. ",
                 "Ensure primed_fetch_pgs_catalog completed successfully. ",
                 "Original error: ", conditionMessage(e))
        }
    )

    rows <- tbl[tbl$pgs_model_id == pgs_id, ]
    if (nrow(rows) == 0) {
        stop("No entry found in 'pgs_scoring_file' table for pgs_id '",
             pgs_id, "'. Ensure primed_fetch_pgs_catalog completed successfully.")
    }

    # Return the file path for the data file type
    data_rows <- rows[rows$file_type == "data", ]
    if (nrow(data_rows) == 0) {
        stop("No 'data' type file found in 'pgs_scoring_file' table for ",
             "pgs_id '", pgs_id, "'.")
    }

    data_rows$file_path[1]
}
