#' Submit the gaudi_prep workflow
#'
#' Configures and submits the \code{gaudi_prep} workflow in the current AnVIL
#' workspace. This workflow converts an array of per-chromosome VCF files to
#' PLINK2 format, infers local ancestry with FLARE, and converts the local
#' ancestry output to the \code{.lanc} format required by
#' \code{\link{submit_make_fbm_workflow}}.
#'
#' The workflow method configuration named \code{gaudi_prep} must already be
#' imported into the workspace from Dockstore
#' (\code{github.com/UW-GAC/gaudi_prep_wdl/gaudi_prep:gaudi_prep_wdl})
#' before calling this function.
#'
#' @param vcf_files Character vector. Google Cloud Storage paths to the
#'   per-chromosome VCF files. These are used both as input to PLINK2 conversion
#'   and as the target files for FLARE local ancestry inference.
#' @param ref_file_list Character vector. Google Cloud Storage paths to the
#'   per-chromosome reference VCF files for FLARE.
#' @param out_prefix_list Character vector. Output prefixes for FLARE, one per
#'   chromosome (e.g. \code{c("chr1", "chr2", ..., "chr22")}).
#' @param genetic_map_file Character. Google Cloud Storage path to the genetic
#'   map file required by FLARE.
#' @param reference_map_file Character. Google Cloud Storage path to the
#'   reference population map file for FLARE.
#' @param samples_keep Character or \code{NULL}. Optional Google Cloud Storage
#'   path to a file containing sample IDs to retain in both the PLINK2
#'   conversion and the FLARE run.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configuration. Defaults to \code{workspace_namespace}.
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   the workflow submission. Default is \code{TRUE}.
#' @param skip_if_complete Logical. If \code{TRUE}, the submission is skipped
#'   when a prior successful run for \code{gaudi_prep} already exists in the
#'   workspace, and the existing submission ID is returned instead. Default is
#'   \code{FALSE}.
#'
#' @return Character. The submission ID of the workflow run.
#'
#' @examples
#' \dontrun{
#' submission_id <- submit_gaudi_prep_workflow(
#'   vcf_files = paste0("gs://my-bucket/vcf/chr", 1:22, ".vcf.gz"),
#'   ref_file_list = paste0("gs://my-bucket/ref/chr", 1:22, "REF.vcf.gz"),
#'   out_prefix_list = paste0("chr", 1:22),
#'   genetic_map_file = "gs://my-bucket/ref/genetic_map.map",
#'   reference_map_file = "gs://my-bucket/ref/reference.pop"
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name
#' @importFrom AnVILGCP avworkflow_configuration_get avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_run avworkflow_jobs
#' @export
submit_gaudi_prep_workflow <- function(
    vcf_files,
    ref_file_list,
    out_prefix_list,
    genetic_map_file,
    reference_map_file,
    samples_keep = NULL,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    use_call_cache = TRUE,
    skip_if_complete = FALSE
) {
    if (!is.character(vcf_files) || length(vcf_files) == 0) {
        stop("'vcf_files' must be a non-empty character vector")
    }
    if (!is.character(ref_file_list) || length(ref_file_list) == 0) {
        stop("'ref_file_list' must be a non-empty character vector")
    }
    if (!is.character(out_prefix_list) || length(out_prefix_list) == 0) {
        stop("'out_prefix_list' must be a non-empty character vector")
    }
    if (length(vcf_files) != length(ref_file_list) ||
        length(vcf_files) != length(out_prefix_list)) {
        stop("'vcf_files', 'ref_file_list', and 'out_prefix_list' must all ",
             "have the same length (one entry per chromosome)")
    }

    if (skip_if_complete) {
        existing <- .find_successful_submission(
            "gaudi_prep",
            namespace = workspace_namespace,
            name = workspace_name
        )
        if (!is.null(existing)) {
            message("Found prior successful gaudi_prep submission: ",
                    existing, ". Skipping.")
            return(existing)
        }
    }

    config <- avworkflow_configuration_get(
        config = "gaudi_prep",
        namespace = workspace_namespace,
        name = workspace_name,
        workflow_namespace = workflow_namespace,
        workflow_name = "gaudi_prep"
    )

    inputs <- list(
        "gaudi_prep.vcf_files"         = as.list(vcf_files),
        "gaudi_prep.ref_file_list"     = as.list(ref_file_list),
        "gaudi_prep.out_prefix_list"   = as.list(out_prefix_list),
        "gaudi_prep.genetic_map_file"  = genetic_map_file,
        "gaudi_prep.reference_map_file" = reference_map_file
    )

    if (!is.null(samples_keep)) {
        inputs[["gaudi_prep.samples_keep"]] <- samples_keep
    }

    config <- avworkflow_configuration_update(config, inputs = inputs)
    result <- avworkflow_run(config,
                             namespace = workspace_namespace,
                             name = workspace_name,
                             submit = TRUE,
                             useCallCache = use_call_cache)
    result$submissionId
}


#' Submit the make_fbm workflow
#'
#' Configures and submits the \code{make_fbm} workflow in the current AnVIL
#' workspace. This workflow converts per-chromosome \code{.lanc} local ancestry
#' files and the corresponding PLINK2 files into a Filebacked Big Matrix (FBM)
#' compatible with the HAUDI and GAUDI methods.
#'
#' The workflow method configuration named \code{make_fbm} must already be
#' imported into the workspace from Dockstore
#' (\code{github.com/frankp-0/HAUDI_workflow/make_fbm:main})
#' before calling this function.
#'
#' @param lanc_files Character vector. Google Cloud Storage paths to the
#'   per-chromosome \code{.lanc} local ancestry files. These can be produced by
#'   the \code{gaudi_prep} or \code{convert_lanc} workflows.
#' @param pgen_files Character vector. Google Cloud Storage paths to the
#'   per-chromosome PLINK2 \code{.pgen} files.
#' @param pvar_files Character vector. Google Cloud Storage paths to the
#'   per-chromosome PLINK2 \code{.pvar} files.
#' @param psam_files Character vector. Google Cloud Storage paths to the
#'   per-chromosome PLINK2 \code{.psam} files.
#' @param fbm_prefix Character. Output prefix for the FBM files
#'   (e.g. \code{"cohort"}). The backing file will be named
#'   \code{<fbm_prefix>.bk}.
#' @param anc_names Character vector. Ancestry names in the same order as the
#'   integer codes (\code{0, 1, ...}) used in the \code{.lanc} files
#'   (e.g. \code{c("AFR", "EUR")}).
#' @param variants_file Character or \code{NULL}. Optional Google Cloud Storage
#'   path to a file with one variant ID per line used to subset the FBM.
#' @param min_ac Integer or \code{NULL}. Optional minimum allele count to
#'   retain a column (ancestry-specific or total genotype) in the FBM.
#' @param samples_file Character or \code{NULL}. Optional Google Cloud Storage
#'   path to a file with one sample ID per line used to subset the FBM.
#' @param chunk_size Integer. Maximum number of variants to read from the
#'   \code{.pgen} file at a time. Default is \code{400}.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configuration. Defaults to \code{workspace_namespace}.
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   the workflow submission. Default is \code{TRUE}.
#' @param skip_if_complete Logical. If \code{TRUE}, the submission is skipped
#'   when a prior successful run for \code{make_fbm} already exists in the
#'   workspace, and the existing submission ID is returned instead. Default is
#'   \code{FALSE}.
#'
#' @return Character. The submission ID of the workflow run.
#'
#' @examples
#' \dontrun{
#' submission_id <- submit_make_fbm_workflow(
#'   lanc_files  = paste0("gs://my-bucket/lanc/chr", 1:22, ".lanc"),
#'   pgen_files  = paste0("gs://my-bucket/plink/chr", 1:22, ".pgen"),
#'   pvar_files  = paste0("gs://my-bucket/plink/chr", 1:22, ".pvar"),
#'   psam_files  = paste0("gs://my-bucket/plink/chr", 1:22, ".psam"),
#'   fbm_prefix  = "cohort",
#'   anc_names   = c("AFR", "EUR")
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name
#' @importFrom AnVILGCP avworkflow_configuration_get avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_run avworkflow_jobs
#' @export
submit_make_fbm_workflow <- function(
    lanc_files,
    pgen_files,
    pvar_files,
    psam_files,
    fbm_prefix,
    anc_names,
    variants_file = NULL,
    min_ac = NULL,
    samples_file = NULL,
    chunk_size = 400L,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    use_call_cache = TRUE,
    skip_if_complete = FALSE
) {
    if (!is.character(lanc_files) || length(lanc_files) == 0) {
        stop("'lanc_files' must be a non-empty character vector")
    }
    if (!is.character(pgen_files) || length(pgen_files) == 0) {
        stop("'pgen_files' must be a non-empty character vector")
    }
    if (!is.character(pvar_files) || length(pvar_files) == 0) {
        stop("'pvar_files' must be a non-empty character vector")
    }
    if (!is.character(psam_files) || length(psam_files) == 0) {
        stop("'psam_files' must be a non-empty character vector")
    }
    n <- length(lanc_files)
    if (length(pgen_files) != n || length(pvar_files) != n ||
        length(psam_files) != n) {
        stop("'lanc_files', 'pgen_files', 'pvar_files', and 'psam_files' must ",
             "all have the same length (one entry per chromosome)")
    }
    if (!is.character(anc_names) || length(anc_names) < 2) {
        stop("'anc_names' must be a character vector with at least two ancestry names")
    }

    if (skip_if_complete) {
        existing <- .find_successful_submission(
            "make_fbm",
            namespace = workspace_namespace,
            name = workspace_name
        )
        if (!is.null(existing)) {
            message("Found prior successful make_fbm submission: ",
                    existing, ". Skipping.")
            return(existing)
        }
    }

    config <- avworkflow_configuration_get(
        config = "make_fbm",
        namespace = workspace_namespace,
        name = workspace_name,
        workflow_namespace = workflow_namespace,
        workflow_name = "make_fbm"
    )

    inputs <- list(
        "make_fbm.lanc_files"  = as.list(lanc_files),
        "make_fbm.pgen_files"  = as.list(pgen_files),
        "make_fbm.pvar_files"  = as.list(pvar_files),
        "make_fbm.psam_files"  = as.list(psam_files),
        "make_fbm.fbm_prefix"  = fbm_prefix,
        "make_fbm.anc_names"   = as.list(anc_names),
        "make_fbm.chunk_size"  = as.integer(chunk_size)
    )

    if (!is.null(variants_file)) {
        inputs[["make_fbm.variants_file"]] <- variants_file
    }
    if (!is.null(min_ac)) {
        inputs[["make_fbm.min_ac"]] <- as.integer(min_ac)
    }
    if (!is.null(samples_file)) {
        inputs[["make_fbm.samples_file"]] <- samples_file
    }

    config <- avworkflow_configuration_update(config, inputs = inputs)
    result <- avworkflow_run(config,
                             namespace = workspace_namespace,
                             name = workspace_name,
                             submit = TRUE,
                             useCallCache = use_call_cache)
    result$submissionId
}


#' Submit the fit_haudi workflow
#'
#' Configures and submits the \code{fit_haudi} workflow in the current AnVIL
#' workspace. This workflow fits a HAUDI or GAUDI polygenic score model using
#' a Filebacked Big Matrix (FBM) produced by
#' \code{\link{submit_make_fbm_workflow}} and a phenotype file.
#'
#' The workflow method configuration named \code{fit_haudi} must already be
#' imported into the workspace from Dockstore
#' (\code{github.com/frankp-0/HAUDI_workflow/fit_haudi:main})
#' before calling this function.
#'
#' @param method Character. PGS method to use: either \code{"HAUDI"} or
#'   \code{"GAUDI"}.
#' @param bk_file Character. Google Cloud Storage path to the FBM backing file
#'   (\code{.bk}) produced by \code{\link{submit_make_fbm_workflow}}.
#' @param info_file Character. Google Cloud Storage path to the FBM column info
#'   file produced by \code{\link{submit_make_fbm_workflow}}.
#' @param dims_file Character. Google Cloud Storage path to the FBM dimensions
#'   file produced by \code{\link{submit_make_fbm_workflow}}.
#' @param fbm_samples_file Character. Google Cloud Storage path to the FBM
#'   samples file produced by \code{\link{submit_make_fbm_workflow}}.
#' @param phenotype_file Character. Google Cloud Storage path to a
#'   tab/space/comma-separated phenotype file. Must contain a column
#'   \code{"#IID"} with sample IDs and at least one phenotype column.
#' @param phenotype Character. Name of the phenotype column in
#'   \code{phenotype_file} to use as the response variable.
#' @param output_prefix Character. Prefix for output files
#'   (model, effects, and PGS results).
#' @param family Character or \code{NULL}. Model family: \code{"gaussian"}
#'   (default) for continuous phenotypes or \code{"binomial"} for binary
#'   phenotypes (HAUDI only). Defaults to \code{NULL} (workflow default:
#'   \code{"gaussian"}).
#' @param training_samples_file Character or \code{NULL}. Optional Google Cloud
#'   Storage path to a file with training sample IDs (one per line).
#' @param gamma_min Numeric. Minimum value for the gamma tuning parameter.
#'   Default is \code{0.01}.
#' @param gamma_max Numeric. Maximum value for the gamma tuning parameter.
#'   Default is \code{5}.
#' @param n_gamma Numeric. Number of gamma values to test.
#'   Default is \code{5}.
#' @param variants_file Character or \code{NULL}. Optional Google Cloud Storage
#'   path to a file with one variant ID per line used to subset variants for
#'   model fitting.
#' @param n_folds Integer. Number of cross-validation folds. Default is
#'   \code{5}.
#' @param workspace_namespace Character. AnVIL workspace namespace.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_namespace}()}.
#' @param workspace_name Character. AnVIL workspace name.
#'   Defaults to \code{\link[AnVILGCP]{avworkspace_name}()}.
#' @param workflow_namespace Character. Namespace of the workflow method
#'   configuration. Defaults to \code{workspace_namespace}.
#' @param use_call_cache Logical. Whether to enable Cromwell call caching for
#'   the workflow submission. Default is \code{TRUE}.
#' @param skip_if_complete Logical. If \code{TRUE}, the submission is skipped
#'   when a prior successful run for \code{fit_haudi} already exists in the
#'   workspace, and the existing submission ID is returned instead. Default is
#'   \code{FALSE}.
#'
#' @return Character. The submission ID of the workflow run.
#'
#' @examples
#' \dontrun{
#' submission_id <- submit_fit_haudi_workflow(
#'   method         = "HAUDI",
#'   bk_file        = "gs://my-bucket/fbm/cohort.bk",
#'   info_file      = "gs://my-bucket/fbm/cohort_info.txt",
#'   dims_file      = "gs://my-bucket/fbm/cohort_dims.txt",
#'   fbm_samples_file = "gs://my-bucket/fbm/cohort_samples.txt",
#'   phenotype_file = "gs://my-bucket/pheno/cohort.pheno",
#'   phenotype      = "BMI",
#'   output_prefix  = "cohort_BMI"
#' )
#' }
#'
#' @importFrom AnVILGCP avworkspace_namespace avworkspace_name
#' @importFrom AnVILGCP avworkflow_configuration_get avworkflow_configuration_update
#' @importFrom AnVILGCP avworkflow_run avworkflow_jobs
#' @export
submit_fit_haudi_workflow <- function(
    method,
    bk_file,
    info_file,
    dims_file,
    fbm_samples_file,
    phenotype_file,
    phenotype,
    output_prefix,
    family = NULL,
    training_samples_file = NULL,
    gamma_min = 0.01,
    gamma_max = 5,
    n_gamma = 5,
    variants_file = NULL,
    n_folds = 5L,
    workspace_namespace = avworkspace_namespace(),
    workspace_name = avworkspace_name(),
    workflow_namespace = workspace_namespace,
    use_call_cache = TRUE,
    skip_if_complete = FALSE
) {
    .validate_haudi_method(method)

    if (!is.null(family)) {
        valid_families <- c("gaussian", "binomial")
        if (!family %in% valid_families) {
            stop("'family' must be one of: ",
                 paste(valid_families, collapse = ", "))
        }
        if (method == "GAUDI" && family == "binomial") {
            stop("'binomial' family is not supported by the GAUDI method")
        }
    }

    if (skip_if_complete) {
        existing <- .find_successful_submission(
            "fit_haudi",
            namespace = workspace_namespace,
            name = workspace_name
        )
        if (!is.null(existing)) {
            message("Found prior successful fit_haudi submission: ",
                    existing, ". Skipping.")
            return(existing)
        }
    }

    config <- avworkflow_configuration_get(
        config = "fit_haudi",
        namespace = workspace_namespace,
        name = workspace_name,
        workflow_namespace = workflow_namespace,
        workflow_name = "fit_haudi"
    )

    inputs <- list(
        "fit_haudi.method"           = method,
        "fit_haudi.bk_file"          = bk_file,
        "fit_haudi.info_file"        = info_file,
        "fit_haudi.dims_file"        = dims_file,
        "fit_haudi.fbm_samples_file" = fbm_samples_file,
        "fit_haudi.phenotype_file"   = phenotype_file,
        "fit_haudi.phenotype"        = phenotype,
        "fit_haudi.output_prefix"    = output_prefix,
        "fit_haudi.gamma_min"        = gamma_min,
        "fit_haudi.gamma_max"        = gamma_max,
        "fit_haudi.n_gamma"          = n_gamma,
        "fit_haudi.n_folds"          = as.integer(n_folds)
    )

    if (!is.null(family)) {
        inputs[["fit_haudi.family"]] <- family
    }
    if (!is.null(training_samples_file)) {
        inputs[["fit_haudi.training_samples_file"]] <- training_samples_file
    }
    if (!is.null(variants_file)) {
        inputs[["fit_haudi.variants_file"]] <- variants_file
    }

    config <- avworkflow_configuration_update(config, inputs = inputs)
    result <- avworkflow_run(config,
                             namespace = workspace_namespace,
                             name = workspace_name,
                             submit = TRUE,
                             useCallCache = use_call_cache)
    result$submissionId
}


# ---- Internal helpers --------------------------------------------------------

#' @keywords internal
.validate_haudi_method <- function(method) {
    valid <- c("HAUDI", "GAUDI")
    if (!is.character(method) || length(method) != 1 || !method %in% valid) {
        stop("'method' must be one of: ", paste(valid, collapse = ", "))
    }
}
