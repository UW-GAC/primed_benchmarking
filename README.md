# primed_benchmarking

R package with tools for running PRIMED PGS (Polygenic Score) benchmarking
workflows in an AnVIL/Terra workspace.

## Overview

The `primedtools` package provides functions for orchestrating the PRIMED PGS
pipeline within an AnVIL workspace:

1. **Fetch a scoring file** from the [PGS Catalog](https://www.pgscatalog.org/)
   using the `primed_fetch_pgs_catalog` workflow.
2. **Calculate individual-level scores** on a cohort using the
   `primed_calc_pgs` workflow.

Both workflows validate their outputs against the
[PRIMED PGS data model](https://github.com/UW-GAC/primed_data_models) and
import results into the workspace data tables.

## Installation

### Stable release

Once the package is published, install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("UW-GAC/primed_benchmarking")
```

### Development version

While the package is under active development, you can install directly from
a specific branch, for example `initial_workspace_setup_workflow_configuration`:

```r
# install.packages("remotes")
remotes::install_github("UW-GAC/primed_benchmarking",
                         ref = "initial_workspace_setup_workflow_configuration")
```

### Installing from a local clone

If you have cloned the repository and want to use the latest uncommitted
changes (for example, while developing or testing), use `devtools`:

```r
# install.packages("devtools")
devtools::load_all("/path/to/primed_benchmarking")
```

`load_all()` sources all R files in `R/` into your session without fully
building the package, which makes the edit-test cycle faster. To do a full
local install instead, use:

```r
devtools::install("/path/to/primed_benchmarking")
```

### Dependencies

The package depends on the [AnVILGCP](https://bioconductor.org/packages/AnVILGCP/)
Bioconductor package. Install it with BiocManager
before installing `primedtools`:

```r
# install.packages("BiocManager")
BiocManager::install("AnVILGCP")
```

### Running the tests

Tests use the [testthat](https://testthat.r-lib.org/) framework. To run them
from a local clone:

```r
# install.packages("devtools")
devtools::test("/path/to/primed_benchmarking")
```

## Prerequisites

### Workspace setup

Before using these functions, the following must be configured in your AnVIL
workspace:

1. **Workflow method configurations** — both `primed_fetch_pgs_catalog` and
   `primed_calc_pgs` must be imported into the workspace from Dockstore or the
   Broad Methods Repository.

2. **Cohort genotype files** — the paths to the cohort PLINK2 files must be
   stored as workspace-level data attributes:
   - `workspace.pgen` — path to the `.pgen` file
   - `workspace.psam` — path to the `.psam` file
   - `workspace.pvar` — path to the `.pvar` file

   > **Note:** The `.pvar` file must have variant IDs in the form
   > `chr:pos:ref:alt` without the `chr` prefix.

   These can be set in the Terra UI under **Data → Workspace Data**, or
   programmatically with `AnVILGCP::avdata_import()`.

## Usage

### Run the full PGS pipeline

```r
library(primed_benchmarking)

result <- run_pgs_pipeline(
  pgs_id         = "PGS000001",
  genome_build   = "GRCh38",
  dest_bucket    = "gs://my-bucket/pgs_results",
  sampleset_name = "my_cohort"
)

# Returns submission IDs for both workflows
result$fetch_submission
result$calc_submission
```

### Run with ancestry adjustment

```r
result <- run_pgs_pipeline(
  pgs_id          = "PGS000001",
  genome_build    = "GRCh38",
  dest_bucket     = "gs://my-bucket/pgs_results",
  sampleset_name  = "my_cohort",
  ancestry_adjust = TRUE,
  pcs             = "gs://my-bucket/cohort/cohort.pcs"
)
```

### Step-by-step usage

```r
# 1. Read cohort file paths from workspace attributes
cohort <- get_cohort_files()

# 2. Submit the fetch workflow
fetch_id <- submit_fetch_pgs_workflow(
  pgs_id       = "PGS000001",
  genome_build = "GRCh38",
  dest_bucket  = "gs://my-bucket/pgs_catalog",
  model_url    = paste0(
    "https://raw.githubusercontent.com/UW-GAC/primed_data_models/",
    "refs/heads/main/PRIMED_PGS_data_model.json"
  )
)

# 3. Wait for the fetch workflow to complete
wait_for_workflow(fetch_id)

# 4. Submit the calc workflow
calc_id <- submit_calc_pgs_workflow(
  pgs_model_id   = "PGS000001",
  scorefile      = "gs://my-bucket/pgs_catalog/PGS000001_hmPOS_GRCh38.txt.gz",
  genome_build   = "GRCh38",
  pgen           = cohort$pgen,
  psam           = cohort$psam,
  pvar           = cohort$pvar,
  min_overlap    = 0.75,
  sampleset_name = "my_cohort",
  dest_bucket    = "gs://my-bucket/pgs_results",
  model_url      = paste0(
    "https://raw.githubusercontent.com/UW-GAC/primed_data_models/",
    "refs/heads/main/PRIMED_PGS_data_model.json"
  )
)
```

## Function reference

### PRIMED PGS pipeline

| Function | Description |
|----------|-------------|
| `run_pgs_pipeline()` | Run the full PGS pipeline for a given PGS Catalog ID |
| `get_cohort_files()` | Read cohort pgen/psam/pvar paths from workspace data attributes |
| `submit_fetch_pgs_workflow()` | Submit the `primed_fetch_pgs_catalog` workflow |
| `submit_calc_pgs_workflow()` | Submit the `primed_calc_pgs` workflow |
| `wait_for_workflow()` | Poll a workflow submission until it completes |

### HAUDI/GAUDI pipeline

| Function | Description |
|----------|-------------|
| `submit_gaudi_prep_workflow()` | Submit the `gaudi_prep` workflow (VCF → PLINK2 + FLARE → .lanc) |
| `submit_make_fbm_workflow()` | Submit the `make_fbm` workflow (.lanc + PLINK2 → FBM) |
| `submit_fit_haudi_workflow()` | Submit the `fit_haudi` workflow (FBM + phenotype → HAUDI/GAUDI PGS model) |

## Workflows

### PRIMED PGS catalog workflows

This package wraps two WDL workflows from the
[primed-pgs-catalog](https://github.com/UW-GAC/primed-pgs-catalog) repository:

- **`primed_fetch_pgs_catalog`** — fetches a scoring file from the PGS Catalog
  and imports metadata into the workspace `pgs_model` and `pgs_scoring_file`
  data tables.
- **`primed_calc_pgs`** — applies the scoring file to cohort genotype data
  (pgen/psam/pvar format) and imports individual-level scores into the
  workspace `pgs_individual_file` data table.

### HAUDI/GAUDI workflows

This package also wraps three WDL workflows for running the
[HAUDI](https://github.com/frankp-0/HAUDI) and GAUDI ancestry-aware PGS methods:

- **`gaudi_prep`**
  ([github.com/UW-GAC/gaudi_prep_wdl](https://github.com/UW-GAC/gaudi_prep_wdl),
  branch `gaudi_prep_wdl`) — converts per-chromosome VCF files to PLINK2
  format, runs FLARE local ancestry inference, and converts FLARE output to the
  `.lanc` format required by `make_fbm`.
- **`make_fbm`**
  ([github.com/frankp-0/HAUDI_workflow](https://github.com/frankp-0/HAUDI_workflow),
  branch `main`) — converts `.lanc` local ancestry files and the matching PLINK2
  files into a Filebacked Big Matrix (FBM) compatible with HAUDI and GAUDI.
- **`fit_haudi`**
  ([github.com/frankp-0/HAUDI_workflow](https://github.com/frankp-0/HAUDI_workflow),
  branch `main`) — fits a HAUDI or GAUDI PGS model using an FBM and a phenotype
  file; outputs ancestry-specific effect estimates and individual-level scores.

#### HAUDI/GAUDI usage example

```r
library(primedtools)

# Step 1 — prepare PLINK2 + .lanc files from VCF inputs
prep_id <- submit_gaudi_prep_workflow(
  vcf_files          = paste0("gs://my-bucket/vcf/chr", 1:22, ".vcf.gz"),
  ref_file_list      = paste0("gs://my-bucket/ref/chr", 1:22, "REF.vcf.gz"),
  out_prefix_list    = paste0("chr", 1:22),
  genetic_map_file   = "gs://my-bucket/ref/genetic_map.map",
  reference_map_file = "gs://my-bucket/ref/reference.pop"
)
wait_for_workflow(prep_id)

# Step 2 — build the Filebacked Big Matrix
#   (supply the .lanc and PLINK2 outputs from Step 1)
fbm_id <- submit_make_fbm_workflow(
  lanc_files = paste0("gs://my-bucket/lanc/chr", 1:22, ".lanc"),
  pgen_files = paste0("gs://my-bucket/plink/chr", 1:22, ".pgen"),
  pvar_files = paste0("gs://my-bucket/plink/chr", 1:22, ".pvar"),
  psam_files = paste0("gs://my-bucket/plink/chr", 1:22, ".psam"),
  fbm_prefix = "cohort",
  anc_names  = c("AFR", "EUR")
)
wait_for_workflow(fbm_id)

# Step 3 — fit the HAUDI/GAUDI model
#   (supply the FBM outputs from Step 2)
fit_id <- submit_fit_haudi_workflow(
  method           = "HAUDI",
  bk_file          = "gs://my-bucket/fbm/cohort.bk",
  info_file        = "gs://my-bucket/fbm/cohort_info.txt",
  dims_file        = "gs://my-bucket/fbm/cohort_dims.txt",
  fbm_samples_file = "gs://my-bucket/fbm/cohort_samples.txt",
  phenotype_file   = "gs://my-bucket/pheno/cohort.pheno",
  phenotype        = "BMI",
  output_prefix    = "cohort_BMI"
)
wait_for_workflow(fit_id)
```
