context("HAUDI/GAUDI workflow submission and validation")

# ---- .validate_haudi_method -------------------------------------------------

test_that(".validate_haudi_method accepts valid methods", {
    expect_silent(primedtools:::.validate_haudi_method("HAUDI"))
    expect_silent(primedtools:::.validate_haudi_method("GAUDI"))
})

test_that(".validate_haudi_method rejects invalid methods", {
    expect_error(primedtools:::.validate_haudi_method("haudi"),
                 "must be one of")
    expect_error(primedtools:::.validate_haudi_method("BOTH"),
                 "must be one of")
    expect_error(primedtools:::.validate_haudi_method(1),
                 "must be one of")
    expect_error(primedtools:::.validate_haudi_method(character(0)),
                 "must be one of")
})

# ---- submit_gaudi_prep_workflow input validation ----------------------------

test_that("submit_gaudi_prep_workflow errors on empty vcf_files", {
    expect_error(
        submit_gaudi_prep_workflow(
            vcf_files = character(0),
            ref_file_list = "gs://b/ref.vcf.gz",
            out_prefix_list = "chr1",
            genetic_map_file = "gs://b/map.map",
            reference_map_file = "gs://b/ref.pop"
        ),
        "non-empty character vector"
    )
})

test_that("submit_gaudi_prep_workflow errors on non-character vcf_files", {
    expect_error(
        submit_gaudi_prep_workflow(
            vcf_files = 1:3,
            ref_file_list = paste0("gs://b/chr", 1:3, "REF.vcf.gz"),
            out_prefix_list = paste0("chr", 1:3),
            genetic_map_file = "gs://b/map.map",
            reference_map_file = "gs://b/ref.pop"
        ),
        "non-empty character vector"
    )
})

test_that("submit_gaudi_prep_workflow errors when array lengths differ", {
    expect_error(
        submit_gaudi_prep_workflow(
            vcf_files       = paste0("gs://b/chr", 1:3, ".vcf.gz"),
            ref_file_list   = paste0("gs://b/chr", 1:2, "REF.vcf.gz"),
            out_prefix_list = paste0("chr", 1:3),
            genetic_map_file = "gs://b/map.map",
            reference_map_file = "gs://b/ref.pop"
        ),
        "same length"
    )
})

# ---- submit_make_fbm_workflow input validation ------------------------------

test_that("submit_make_fbm_workflow errors on empty lanc_files", {
    expect_error(
        submit_make_fbm_workflow(
            lanc_files = character(0),
            pgen_files = "gs://b/chr1.pgen",
            pvar_files = "gs://b/chr1.pvar",
            psam_files = "gs://b/chr1.psam",
            fbm_prefix = "cohort",
            anc_names  = c("AFR", "EUR")
        ),
        "non-empty character vector"
    )
})

test_that("submit_make_fbm_workflow errors when file arrays have different lengths", {
    expect_error(
        submit_make_fbm_workflow(
            lanc_files = paste0("gs://b/chr", 1:3, ".lanc"),
            pgen_files = paste0("gs://b/chr", 1:2, ".pgen"),
            pvar_files = paste0("gs://b/chr", 1:3, ".pvar"),
            psam_files = paste0("gs://b/chr", 1:3, ".psam"),
            fbm_prefix = "cohort",
            anc_names  = c("AFR", "EUR")
        ),
        "same length"
    )
})

test_that("submit_make_fbm_workflow errors with fewer than two ancestry names", {
    expect_error(
        submit_make_fbm_workflow(
            lanc_files = "gs://b/chr1.lanc",
            pgen_files = "gs://b/chr1.pgen",
            pvar_files = "gs://b/chr1.pvar",
            psam_files = "gs://b/chr1.psam",
            fbm_prefix = "cohort",
            anc_names  = "AFR"
        ),
        "at least two ancestry names"
    )
})

test_that("submit_make_fbm_workflow errors with non-character anc_names", {
    expect_error(
        submit_make_fbm_workflow(
            lanc_files = "gs://b/chr1.lanc",
            pgen_files = "gs://b/chr1.pgen",
            pvar_files = "gs://b/chr1.pvar",
            psam_files = "gs://b/chr1.psam",
            fbm_prefix = "cohort",
            anc_names  = c(0, 1)
        ),
        "at least two ancestry names"
    )
})

# ---- submit_fit_haudi_workflow input validation -----------------------------

test_that("submit_fit_haudi_workflow errors on invalid method", {
    expect_error(
        submit_fit_haudi_workflow(
            method           = "LMM",
            bk_file          = "gs://b/cohort.bk",
            info_file        = "gs://b/cohort_info.txt",
            dims_file        = "gs://b/cohort_dims.txt",
            fbm_samples_file = "gs://b/cohort_samples.txt",
            phenotype_file   = "gs://b/cohort.pheno",
            phenotype        = "BMI",
            output_prefix    = "out"
        ),
        "must be one of"
    )
})

test_that("submit_fit_haudi_workflow errors on invalid family", {
    expect_error(
        submit_fit_haudi_workflow(
            method           = "HAUDI",
            bk_file          = "gs://b/cohort.bk",
            info_file        = "gs://b/cohort_info.txt",
            dims_file        = "gs://b/cohort_dims.txt",
            fbm_samples_file = "gs://b/cohort_samples.txt",
            phenotype_file   = "gs://b/cohort.pheno",
            phenotype        = "BMI",
            output_prefix    = "out",
            family           = "poisson"
        ),
        "must be one of"
    )
})

test_that("submit_fit_haudi_workflow errors when GAUDI is used with binomial family", {
    expect_error(
        submit_fit_haudi_workflow(
            method           = "GAUDI",
            bk_file          = "gs://b/cohort.bk",
            info_file        = "gs://b/cohort_info.txt",
            dims_file        = "gs://b/cohort_dims.txt",
            fbm_samples_file = "gs://b/cohort_samples.txt",
            phenotype_file   = "gs://b/cohort.pheno",
            phenotype        = "T2D",
            output_prefix    = "out",
            family           = "binomial"
        ),
        "not supported by the GAUDI method"
    )
})

test_that("submit_fit_haudi_workflow accepts NULL family (uses workflow default)", {
    with_mocked_bindings(
        avworkflow_configuration_get   = function(...) list(),
        avworkflow_configuration_update = function(config, inputs) config,
        avworkflow_run = function(config, ...) list(submissionId = "sub-null-family"),
        .package = "AnVIL",
        {
            result <- submit_fit_haudi_workflow(
                method           = "HAUDI",
                bk_file          = "gs://b/cohort.bk",
                info_file        = "gs://b/cohort_info.txt",
                dims_file        = "gs://b/cohort_dims.txt",
                fbm_samples_file = "gs://b/cohort_samples.txt",
                phenotype_file   = "gs://b/cohort.pheno",
                phenotype        = "BMI",
                output_prefix    = "out",
                family           = NULL,
                workspace_namespace = "test-ns",
                workspace_name   = "test-ws"
            )
            expect_equal(result, "sub-null-family")
        }
    )
})
