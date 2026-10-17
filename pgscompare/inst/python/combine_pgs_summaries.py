"""
Combine <DATASET>_summary.csv files across all PGS scores for a given cohort.

Run this INSIDE an AnVIL/Terra Jupyter environment (not locally), since the
match_summary files live in a secure GCS bucket that your notebook's pet
service account has access to.

Generalizable across workspaces/cohorts: just edit the CONFIG block below.
WORKSPACE_NAMESPACE / WORKSPACE_NAME / WORKSPACE_BUCKET are picked up
automatically from the environment of whatever workspace you're running in.

Requires: pandas, firecloud (fiss), and `gsutil` on PATH — all already
present in standard Terra/AnVIL Jupyter images.
"""

import io
import os
import subprocess
import pandas as pd
import firecloud.api as fapi

# =============================== CONFIG ======================================
DATASET_NAME = "GeneSTAR"          # cohort name, e.g. "GeneSTAR", "WHI", "MESA"...
ENTITY_TYPE = "pgs_scoring_file"   # data table name (assumed same across workspaces)
BUCKET_SUBFOLDER = "pgs_combined_summaries"  # subfolder under the workspace bucket
MATCH_STATUS_FILTER = "matched"    # status used for the flipped/not-flipped split
# ==============================================================================

SUMMARY_FILENAME = f"{DATASET_NAME}_summary.csv"
OUTPUT_COMBINED = f"{DATASET_NAME}_summary_combined.csv"
OUTPUT_FLIPPED = f"{DATASET_NAME}_summary_{MATCH_STATUS_FILTER}_flipped.csv"
OUTPUT_NOT_FLIPPED = f"{DATASET_NAME}_summary_{MATCH_STATUS_FILTER}_not_flipped.csv"

# ---- 1. Pull the pgs_scoring_file data table directly from the workspace --
# Terra/AnVIL notebook runtimes set these env vars automatically for whatever
# workspace the notebook is running in, so no hardcoding needed per-cohort.
WORKSPACE_NAMESPACE = os.environ["WORKSPACE_NAMESPACE"]
WORKSPACE_NAME = os.environ["WORKSPACE_NAME"]

response = fapi.get_entities_tsv(
    WORKSPACE_NAMESPACE, WORKSPACE_NAME, ENTITY_TYPE, model="flexible"
)
response.raise_for_status()

scoring_table_full = pd.read_csv(io.StringIO(response.text), sep="\t")
print(f"Pulled {len(scoring_table_full)} rows from data table '{ENTITY_TYPE}' "
      f"in {WORKSPACE_NAMESPACE}/{WORKSPACE_NAME}")

# Some workspaces' pgs_scoring_file tables may span multiple cohorts/datasets
# (e.g. a shared table with both GeneSTAR and WHI rows). Filter down to just
# the rows whose match_summary points at this dataset's summary file.
scoring_table = scoring_table_full[
    scoring_table_full["match_summary"].str.endswith(SUMMARY_FILENAME)
].copy()
print(f"Filtered to {len(scoring_table)} rows matching '{SUMMARY_FILENAME}'")

if scoring_table.empty:
    raise ValueError(
        f"No rows found with match_summary ending in '{SUMMARY_FILENAME}'. "
        f"Check DATASET_NAME, or inspect scoring_table_full['match_summary'] "
        f"to see what naming convention this workspace actually uses."
    )

# ---- 2. Helper to read a CSV directly from a gs:// path --------------------
def read_gs_csv(gs_path: str) -> pd.DataFrame:
    raw = subprocess.run(
        ["gsutil", "cat", gs_path],
        check=True,
        capture_output=True,
    ).stdout
    return pd.read_csv(io.BytesIO(raw))


# ---- 3. Loop over every row, fetch its summary file, and stack ------------
combined_frames = []
failures = []

for _, row in scoring_table.iterrows():
    pgs_id = row["pgs_model_id"]
    gs_path = row["match_summary"]

    try:
        df = read_gs_csv(gs_path)
    except subprocess.CalledProcessError as e:
        failures.append((pgs_id, gs_path, e.stderr.decode(errors="replace")))
        continue

    # Sanity check: the file's own 'accession' column should match pgs_model_id
    if "accession" in df.columns and not (df["accession"] == pgs_id).all():
        print(f"WARNING: accession mismatch in file for {pgs_id}")

    df["pgs_scoring_file_id"] = row["entity:pgs_scoring_file_id"]
    combined_frames.append(df)

if failures:
    print(f"\n{len(failures)} file(s) failed to download:")
    for pgs_id, path, err in failures:
        print(f"  {pgs_id}: {path}\n    {err.strip()}")

# ---- 4. Concatenate and save ------------------------------------------------
combined = pd.concat(combined_frames, ignore_index=True)
print(f"\nCombined {len(combined_frames)} files -> {combined.shape[0]} rows, "
      f"{combined['accession'].nunique()} unique PGS accessions.")

combined.to_csv(OUTPUT_COMBINED, index=False)
print(f"Saved combined table to {OUTPUT_COMBINED}")

# ---- 5. (Optional) quick sanity summary ------------------------------------
print("\nRow counts by match_status across all PGS scores:")
print(combined.groupby("match_status")["count"].sum())

# Total variant count per PGS score (sum of 'count' across all its
# match_status rows), added as a column so each row shows its share of the
# whole for that accession.
combined["total_count"] = combined.groupby("accession")["count"].transform("sum")

# ---- 6. Filter by match_status / match_flipped -----------------------------
# match_flipped is blank ("") for unmatched rows, so coerce to a real bool
# before filtering rather than comparing against the string "true"/"false".
combined["match_flipped_bool"] = (
    combined["match_flipped"].astype(str).str.lower().map({"true": True, "false": False})
)

def filter_summary(df, match_status=None, match_flipped=None):
    """
    match_status: str or list of str, e.g. "matched" or ["matched", "excluded"]
    match_flipped: True, False, or None (None = don't filter on this column)
    """
    out = df.copy()
    if match_status is not None:
        statuses = [match_status] if isinstance(match_status, str) else match_status
        out = out[out["match_status"].isin(statuses)]
    if match_flipped is not None:
        out = out[out["match_flipped_bool"] == match_flipped]
    return out

matched_flipped = filter_summary(combined, match_status=MATCH_STATUS_FILTER, match_flipped=True)
print(f"\n{MATCH_STATUS_FILTER} + flipped: {len(matched_flipped)} rows, "
      f"{matched_flipped['count'].sum()} total variants across "
      f"{matched_flipped['accession'].nunique()} PGS scores")

matched_not_flipped = filter_summary(combined, match_status=MATCH_STATUS_FILTER, match_flipped=False)
print(f"{MATCH_STATUS_FILTER} + not flipped: {len(matched_not_flipped)} rows, "
      f"{matched_not_flipped['count'].sum()} total variants")

matched_flipped.to_csv(OUTPUT_FLIPPED, index=False)
matched_not_flipped.to_csv(OUTPUT_NOT_FLIPPED, index=False)

# ---- 7. Copy outputs back to the Terra workspace bucket --------------------
# Terra notebook runtimes set WORKSPACE_BUCKET automatically (e.g.
# "gs://fc-secure-xxxxxxxx-xxxx-...").
WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]
DEST_PREFIX = f"{WORKSPACE_BUCKET}/{BUCKET_SUBFOLDER}/"

output_files = [OUTPUT_COMBINED, OUTPUT_FLIPPED, OUTPUT_NOT_FLIPPED]

for f in output_files:
    subprocess.run(["gsutil", "cp", f, DEST_PREFIX], check=True)
    print(f"Uploaded {f} -> {DEST_PREFIX}{f}")
