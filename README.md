# primed.benchmarking

An R package with tools for running PRIMED PGS (Polygenic Score) benchmarking
workflows in an AnVIL/Terra workspace. The package wraps and orchestrates WDL
workflows, and provides helper utilities for workspace inputs/outputs, local
ancestry pipelines (HAUDI/GAUDI), and PGS Catalog metadata ingestion.

## Overview

`primed.benchmarking` supports two main use cases:

1. **PRIMED PGS Catalog pipeline**
   - Fetch a scoring file from the [PGS Catalog](https://www.pgscatalog.org/)
     using the `primed_fetch_pgs_catalog` workflow.
   - Compute individual-level PGS on cohort genotypes (PLINK2 pgen/psam/pvar)
     using the `primed_calc_pgs` workflow.
   - Optionally perform ancestry-based score adjustment using PCs.
   - Both workflows validate outputs against the
     [PRIMED PGS data model](https://github.com/UW-GAC/primed_data_models) and
     import results into the workspace data tables.

2. **HAUDI/GAUDI local-ancestry-aware pipeline**
   - Prepare local ancestry inputs (VCF → PLINK2 + FLARE → `.lanc`).
   - Build a Filebacked Big Matrix (FBM) for HAUDI/GAUDI.
   - Fit HAUDI or GAUDI models and generate ancestry-aware PGS.

Many functions assume you are authenticated to AnVIL/Terra and working inside
(or targeting) a specific workspace.

## Installation

### Stable release

Once the package is published, install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("UW-GAC/primed.benchmarking")
```

### Development version

While the package is under active development, you can install directly from
a specific branch:

```r
# install.packages("remotes")
remotes::install_github("UW-GAC/primed.benchmarking", ref = "main")
```

### Installing from a local clone

If you have cloned the repository and want to use the latest uncommitted
changes (for example, while developing or testing), use `devtools`:

```r
# install.packages("devtools")
devtools::load_all("/path/to/primed.benchmarking")
```

`load_all()` sources all R files in `R/` into your session without fully
building the package, which makes the edit-test cycle faster. To do a full
local install instead, use:

```r
devtools::install("/path/to/primed.benchmarking")
```

### Dependencies

The package depends on the [AnVILGCP](https://bioconductor.org/packages/AnVILGCP/)
Bioconductor package. Install it with BiocManager
before installing `primed.benchmarking`:

```r
# install.packages("BiocManager")
BiocManager::install("AnVILGCP")
```

### Running the tests

Tests use the [testthat](https://testthat.r-lib.org/) framework. To run them
from a local clone:

```r
# install.packages("devtools")
devtools::test("/path/to/primed.benchmarking")
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
library(primed.benchmarking)

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

## Exported API reference

All exported (public) functions are listed below, grouped by category.

---

### PRIMED PGS pipeline

#### `run_pgs_pipeline(pgs_id, genome_build, dest_bucket, sampleset_name, ...)`

Orchestrates the full end-to-end PRIMED PGS workflow for a given PGS Catalog
score:

1. Reads cohort genotype paths from workspace attributes.
2. Submits `primed_fetch_pgs_catalog` to fetch and import the scoring file.
3. Waits for that workflow to complete.
4. Looks up the fetched scoring file path from the `pgs_scoring_file` table.
5. Submits `primed_calc_pgs` to calculate individual-level scores.

| Argument | Type | Default | Description |
|---|---|---|---|
| `pgs_id` | `character` | — | PGS Catalog ID, e.g. `"PGS000001"`. |
| `genome_build` | `character` | — | `"GRCh38"` or `"GRCh37"`. |
| `dest_bucket` | `character` | — | `gs://` path where output files are written. |
| `sampleset_name` | `character` | — | Cohort name used in output file naming. |
| `model_url` | `character` | PRIMED data model URL | URL to the PRIMED PGS data model JSON. |
| `min_overlap` | `numeric` | `0.75` | Minimum fraction of score variants present in the genotype data. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configurations. |
| `overwrite` | `logical` | `FALSE` | Overwrite existing rows in the data tables. |
| `ancestry_adjust` | `logical` | `FALSE` | Enable ancestry-based score adjustment. |
| `pcs` | `character` or `NULL` | `NULL` | `gs://` path to PC file (required when `ancestry_adjust = TRUE`). |
| `primed_dataset_id` | `character` or `NULL` | `NULL` | Optional PRIMED dataset identifier. |
| `poll_interval` | `numeric` | `60` | Seconds between workflow status polls. |
| `timeout` | `numeric` | `3600` | Maximum seconds to wait for the fetch workflow. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Reuse prior successful submissions when available. |

**Returns** a named list with:
- `fetch_submission` — submission ID of the `primed_fetch_pgs_catalog` run.
- `calc_submission` — submission ID of the `primed_calc_pgs` run.

---

#### `get_cohort_files(workspace_namespace, workspace_name)`

Reads cohort genotype file paths stored as workspace-level data attributes
(`workspace.pgen`, `workspace.psam`, `workspace.pvar`) and returns them as a
named list.

| Argument | Type | Default | Description |
|---|---|---|---|
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |

**Returns** a named list with character elements `pgen`, `psam`, and `pvar`
(Google Cloud Storage paths).

---

#### `submit_fetch_pgs_workflow(pgs_id, genome_build, dest_bucket, model_url, ...)`

Configures and submits the `primed_fetch_pgs_catalog` workflow in the current
AnVIL workspace. This workflow fetches a scoring file from the PGS Catalog,
copies it to the specified GCS bucket, and imports metadata into the workspace
`pgs_model` and `pgs_scoring_file` data tables.

The workflow method configuration named `primed_fetch_pgs_catalog` must already
be imported into the workspace before calling this function.

| Argument | Type | Default | Description |
|---|---|---|---|
| `pgs_id` | `character` | — | PGS Catalog ID, e.g. `"PGS000001"`. |
| `genome_build` | `character` | — | `"GRCh38"` or `"GRCh37"`. |
| `dest_bucket` | `character` | — | `gs://` path where scoring files are written. |
| `model_url` | `character` | — | URL to the PRIMED PGS data model JSON. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configuration. |
| `overwrite` | `logical` | `FALSE` | Overwrite existing rows in the data tables. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Skip submission if a prior successful run exists; return its ID. |

**Returns** a character string: the workflow submission ID.

---

#### `submit_calc_pgs_workflow(pgs_model_id, scorefile, genome_build, pgen, psam, pvar, min_overlap, sampleset_name, dest_bucket, model_url, ...)`

Configures and submits the `primed_calc_pgs` workflow. This workflow matches
the scoring file to cohort genotype data, calculates individual-level polygenic
scores with PLINK2, optionally adjusts for ancestry using PCs, and imports
results into the workspace `pgs_individual_file` data table.

The workflow method configuration named `primed_calc_pgs` must already be
imported into the workspace before calling this function.

| Argument | Type | Default | Description |
|---|---|---|---|
| `pgs_model_id` | `character` | — | PGS model identifier, e.g. `"PGS000001"`. |
| `scorefile` | `character` | — | `gs://` path to the scoring file fetched from the PGS Catalog. |
| `genome_build` | `character` | — | `"GRCh38"` or `"GRCh37"`. |
| `pgen` | `character` | — | `gs://` path to the cohort `.pgen` file. |
| `psam` | `character` | — | `gs://` path to the cohort `.psam` file. |
| `pvar` | `character` | — | `gs://` path to the cohort `.pvar` file. |
| `min_overlap` | `numeric` | — | Minimum fraction of score variants present in the genotype data. |
| `sampleset_name` | `character` | — | Name used to construct output file names. |
| `dest_bucket` | `character` | — | `gs://` path where score output files are written. |
| `model_url` | `character` | — | URL to the PRIMED PGS data model JSON. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configuration. |
| `overwrite` | `logical` | `FALSE` | Overwrite existing rows in the data tables. |
| `ancestry_adjust` | `logical` | `FALSE` | Enable ancestry-based score adjustment. |
| `pcs` | `character` or `NULL` | `NULL` | `gs://` path to a PC file for ancestry adjustment. |
| `primed_dataset_id` | `character` or `NULL` | `NULL` | Optional PRIMED dataset identifier. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Skip submission if a prior successful run exists; return its ID. |

**Returns** a character string: the workflow submission ID.

---

#### `wait_for_workflow(submission_id, workspace_namespace, workspace_name, poll_interval, timeout)`

Polls an AnVIL workflow submission at regular intervals until all workflows in
the submission reach a terminal state (`Succeeded`, `Failed`, or `Aborted`).
Raises an error if any workflows fail or are aborted.

| Argument | Type | Default | Description |
|---|---|---|---|
| `submission_id` | `character` | — | Submission ID returned by `submit_fetch_pgs_workflow()` or `submit_calc_pgs_workflow()`. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `poll_interval` | `numeric` | `60` | Seconds between status checks. |
| `timeout` | `numeric` | `3600` | Maximum seconds to wait before timing out. |

**Returns** invisibly: the final job status tibble.

---

### HAUDI/GAUDI workflows

This package wraps three WDL workflows for running the
[HAUDI](https://github.com/frankp-0/HAUDI) and GAUDI ancestry-aware PGS
methods.

#### `submit_gaudi_prep_workflow(vcf_files, ref_file_list, out_prefix_list, genetic_map_file, reference_map_file, ...)`

Configures and submits the `gaudi_prep` workflow
([github.com/UW-GAC/gaudi_prep_wdl](https://github.com/UW-GAC/gaudi_prep_wdl),
branch `gaudi_prep_wdl`). This workflow converts per-chromosome VCF files to
PLINK2 format, runs FLARE local ancestry inference, and converts FLARE output
to the `.lanc` format required by `submit_make_fbm_workflow()`.

The workflow method configuration named `gaudi_prep` must already be imported
into the workspace from Dockstore
(`github.com/UW-GAC/gaudi_prep_wdl/gaudi_prep:gaudi_prep_wdl`) before calling
this function.

| Argument | Type | Default | Description |
|---|---|---|---|
| `vcf_files` | `character` vector | — | `gs://` paths to per-chromosome VCF files. |
| `ref_file_list` | `character` vector | — | `gs://` paths to per-chromosome reference VCF files for FLARE. |
| `out_prefix_list` | `character` vector | — | Output prefixes for FLARE, one per chromosome (e.g. `c("chr1", ..., "chr22")`). |
| `genetic_map_file` | `character` | — | `gs://` path to the genetic map file for FLARE. |
| `reference_map_file` | `character` | — | `gs://` path to the FLARE reference population map file. |
| `samples_keep` | `character` or `NULL` | `NULL` | Optional `gs://` path to a file of sample IDs to retain. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configuration. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Skip submission if a prior successful run exists; return its ID. |

**Returns** a character string: the workflow submission ID.

---

#### `submit_make_fbm_workflow(lanc_files, pgen_files, pvar_files, psam_files, fbm_prefix, anc_names, ...)`

Configures and submits the `make_fbm` workflow
([github.com/frankp-0/HAUDI_workflow](https://github.com/frankp-0/HAUDI_workflow),
branch `main`). This workflow converts per-chromosome `.lanc` local ancestry
files and matching PLINK2 files into a Filebacked Big Matrix (FBM) compatible
with HAUDI and GAUDI.

The workflow method configuration named `make_fbm` must already be imported
into the workspace from Dockstore
(`github.com/frankp-0/HAUDI_workflow/make_fbm:main`) before calling this
function.

| Argument | Type | Default | Description |
|---|---|---|---|
| `lanc_files` | `character` vector | — | `gs://` paths to per-chromosome `.lanc` files. |
| `pgen_files` | `character` vector | — | `gs://` paths to per-chromosome `.pgen` files. |
| `pvar_files` | `character` vector | — | `gs://` paths to per-chromosome `.pvar` files. |
| `psam_files` | `character` vector | — | `gs://` paths to per-chromosome `.psam` files. |
| `fbm_prefix` | `character` | — | Output prefix for the FBM files (the backing file will be `<fbm_prefix>.bk`). |
| `anc_names` | `character` vector | — | Ancestry names in the same order as the integer codes used in the `.lanc` files (e.g. `c("AFR", "EUR")`). |
| `variants_file` | `character` or `NULL` | `NULL` | Optional `gs://` path to a file of variant IDs used to subset the FBM. |
| `min_ac` | `integer` or `NULL` | `NULL` | Optional minimum allele count to retain a column in the FBM. |
| `samples_file` | `character` or `NULL` | `NULL` | Optional `gs://` path to a file of sample IDs used to subset the FBM. |
| `chunk_size` | `integer` | `400` | Maximum number of variants to read from the `.pgen` file at a time. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configuration. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Skip submission if a prior successful run exists; return its ID. |

**Returns** a character string: the workflow submission ID.

---

#### `submit_fit_haudi_workflow(method, bk_file, info_file, dims_file, fbm_samples_file, phenotype_file, phenotype, output_prefix, ...)`

Configures and submits the `fit_haudi` workflow
([github.com/frankp-0/HAUDI_workflow](https://github.com/frankp-0/HAUDI_workflow),
branch `main`). This workflow fits a HAUDI or GAUDI polygenic score model using
the FBM produced by `submit_make_fbm_workflow()` and a phenotype file, and
outputs ancestry-specific effect estimates and individual-level PGS.

The workflow method configuration named `fit_haudi` must already be imported
into the workspace from Dockstore
(`github.com/frankp-0/HAUDI_workflow/fit_haudi:main`) before calling this
function.

| Argument | Type | Default | Description |
|---|---|---|---|
| `method` | `character` | — | `"HAUDI"` or `"GAUDI"`. |
| `bk_file` | `character` | — | `gs://` path to the FBM backing file (`.bk`) from `submit_make_fbm_workflow()`. |
| `info_file` | `character` | — | `gs://` path to the FBM column info file. |
| `dims_file` | `character` | — | `gs://` path to the FBM dimensions file. |
| `fbm_samples_file` | `character` | — | `gs://` path to the FBM samples file. |
| `phenotype_file` | `character` | — | `gs://` path to a phenotype file. Must contain a `"#IID"` column and at least one phenotype column. |
| `phenotype` | `character` | — | Name of the phenotype column to use as the response variable. |
| `output_prefix` | `character` | — | Prefix for output files (model, effects, PGS results). |
| `family` | `character` or `NULL` | `NULL` (→ `"gaussian"`) | Model family: `"gaussian"` or `"binomial"` (HAUDI only). |
| `training_samples_file` | `character` or `NULL` | `NULL` | Optional `gs://` path to a file with training sample IDs. |
| `gamma_min` | `numeric` | `0.01` | Minimum value of the gamma tuning parameter. |
| `gamma_max` | `numeric` | `5` | Maximum value of the gamma tuning parameter. |
| `n_gamma` | `numeric` | `5` | Number of gamma values to evaluate. |
| `variants_file` | `character` or `NULL` | `NULL` | Optional `gs://` path to a file of variant IDs to use for model fitting. |
| `n_folds` | `integer` | `5` | Number of cross-validation folds. |
| `workspace_namespace` | `character` | current workspace | AnVIL workspace namespace. |
| `workspace_name` | `character` | current workspace | AnVIL workspace name. |
| `workflow_namespace` | `character` | `workspace_namespace` | Namespace of the workflow method configuration. |
| `use_call_cache` | `logical` | `TRUE` | Enable Cromwell call caching. |
| `skip_if_complete` | `logical` | `FALSE` | Skip submission if a prior successful run exists; return its ID. |

**Returns** a character string: the workflow submission ID.

---

### Local ancestry summary workflow

#### `create_local_ancestry_summary_table(cohort, cohort.namespace, cohort.name)`

Checks whether the `local_ancestry_summary` workspace data table exists in the
specified AnVIL workspace. If absent, creates a minimal table with a single
`local_ancestry_summary_id` column populated by `cohort`. If it already exists,
emits a message and makes no changes.

| Argument | Type | Description |
|---|---|---|
| `cohort` | `character` | Cohort/entity ID for the initial `local_ancestry_summary_id` value. |
| `cohort.namespace` | `character` | AnVIL workspace namespace. |
| `cohort.name` | `character` | AnVIL workspace name. |

**Returns** invisibly: `FALSE` if the table was created, `TRUE` if it already
existed.

---

#### `set_up_step1c_summarize_local_ancestry_proportions(cohort, cohort.namespace, cohort.name, merged.6.ancestry_frac_path)`

Retrieves the `step1c_summarize_local_ancestry_proportions` workflow
configuration from the workspace, updates its inputs and outputs for the given
cohort and data file, then performs a **dry-run** validation and dry-run
submission. Use
`run_step1c_summarize_local_ancestry_proportions()` to actually submit the job.

| Argument | Type | Description |
|---|---|---|
| `cohort` | `character` | Entity name to run on. |
| `cohort.namespace` | `character` | AnVIL workspace namespace. |
| `cohort.name` | `character` | AnVIL workspace name. |
| `merged.6.ancestry_frac_path` | `character` | `gs://` path to the merged 6-ancestry fraction file from an upstream step. |

**Returns** the updated workflow configuration object.

---

#### `run_step1c_summarize_local_ancestry_proportions(cohort, cohort.namespace, cohort.name, merged.6.ancestry_frac_path, run_now, new_config)`

Applies the workflow configuration and submits the
`step1c_summarize_local_ancestry_proportions` workflow for real when
`run_now = TRUE`. When `run_now = FALSE` (default) returns invisibly without
doing anything, making it safe to call in scripts still being prepared.

| Argument | Type | Default | Description |
|---|---|---|---|
| `cohort` | `character` | — | Entity name to run on. |
| `cohort.namespace` | `character` | — | AnVIL workspace namespace. |
| `cohort.name` | `character` | — | AnVIL workspace name. |
| `merged.6.ancestry_frac_path` | `character` | — | GCS path to the merged 6-ancestry fraction file. Accepted for API symmetry with `set_up_step1c_summarize_local_ancestry_proportions()`; the configuration is already embedded in `new_config` so this argument is not read again here. |
| `run_now` | `logical` | `FALSE` | If `TRUE`, apply the configuration and submit the workflow. |
| `new_config` | workflow config | — | Configuration object returned by `set_up_step1c_summarize_local_ancestry_proportions()`. |

**Returns** invisibly `NULL`.

---

### Ancestry utilities

#### `get_two_way_ancestry(admixture_anc_prop_list, cohort_name, threshold, min_prop)`

Computes two-way ancestry counts from individual-level admixture proportion
data. For every pair of reference populations (columns whose names begin with
`"K"`, e.g. `KAFR`, `KEUR`), the function counts individuals meeting two-way
criteria and related exclusion categories.

| Argument | Type | Default | Description |
|---|---|---|---|
| `admixture_anc_prop_list` | `data.frame` / `tibble` | — | Data frame with at least two numeric columns whose names begin with `"K"`. Each row is one individual. |
| `cohort_name` | `character` | — | Cohort name; populates the `Cohort` column in the result. |
| `threshold` | `numeric` | `0.9` | Minimum combined ancestry proportion (`x1 + x2`) for an individual to be considered for any category. |
| `min_prop` | `numeric` | `0.10` | Minimum individual ancestry proportion for classification as admixed. |

**Returns** a tibble with one row per pair of reference populations and the
following columns:

| Column | Description |
|---|---|
| `Cohort` | Cohort name. |
| `Ref_Pop1` | Name of the first reference population column. |
| `Ref_Pop2` | Name of the second reference population column. |
| `Count_two_way` | Individuals with `x1 >= min_prop`, `x2 >= min_prop`, and `x1 + x2 >= threshold`. |
| `Excluded_Ref1_lt10_and_Ref2_lt90` | Individuals with `x1 < min_prop`, `x2 < threshold`, and `x1 + x2 >= threshold`. |
| `Excluded_Ref2_lt10_and_Ref1_lt90` | Individuals with `x2 < min_prop`, `x1 < threshold`, and `x1 + x2 >= threshold`. |
| `Excluded_Ref1_lt10_and_Ref2_gt90` | Individuals with `x1 < min_prop` and `x2 >= threshold`. |
| `Excluded_Ref2_lt10_and_Ref1_gt90` | Individuals with `x2 < min_prop` and `x1 >= threshold`. |
| `n` | Total individuals with `x1 + x2 >= threshold`. |

---

### PGS Catalog utilities

#### `read_pgs_all_metadata_scores(url, split_list_columns)`

Downloads and cleans the PGS Catalog "all metadata scores" CSV from the EBI
FTP server. Column names are standardized to `snake_case`, whitespace is
trimmed, obviously numeric or date-like columns are converted to their native
types, and optionally columns with delimited list values are split into
list-columns.

| Argument | Type | Default | Description |
|---|---|---|---|
| `url` | `character` | EBI FTP URL | URL of the `pgs_all_metadata_scores.csv` file. |
| `split_list_columns` | `logical` | `TRUE` | Split columns with pipe/semicolon/comma-space delimiters into list-columns. |

**Returns** a tibble with cleaned columns ready for analysis and joins. When
`split_list_columns = TRUE`, some columns will be list-columns (each element a
character vector).

---

### Google Cloud Storage utilities

#### `get_file_from_bucket(gspath, newfilename)`

Downloads a file from a GCS bucket to a local path using `gsutil cp`. If the
destination file already exists locally, the copy is skipped.

Requires `gsutil` (Google Cloud SDK) to be installed and available on `PATH`.

| Argument | Type | Description |
|---|---|---|
| `gspath` | `character` | GCS path of the file to download (e.g. `"gs://my-bucket/path/file.txt"`). |
| `newfilename` | `character` | Local destination path where the file should be saved. |

**Returns** invisibly `NULL`.

---

#### `copy_file_to_bucket(filename, gspath, newfilename)`

Uploads a local file to a GCS bucket using `gsutil cp`, then confirms the
transfer with `gsutil ls -l`.

Requires `gsutil` (Google Cloud SDK) to be installed and available on `PATH`.

| Argument | Type | Description |
|---|---|---|
| `filename` | `character` | Local path of the file to upload. |
| `gspath` | `character` | GCS bucket or prefix path to copy the file into (e.g. `"gs://my-bucket/output"`). |
| `newfilename` | `character` | Name to give the file in the bucket. The file is written to `<gspath>/<newfilename>`. |

**Returns** invisibly `NULL`.

---

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
library(primed.benchmarking)

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

