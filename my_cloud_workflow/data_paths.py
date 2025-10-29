import os

# Base directory on the VM or local machine where your data will be mounted
# You can override this by setting the BASE_DIR environment variable.
# Default points to the local workspace path on macOS used in this project.
GCS_MOUNT_DIR = "/mnt/gcs"
BASE_DIR = os.path.join(GCS_MOUNT_DIR, "my_project_data")

#use os

# If DEFAULT_BASE doesn't exist, try a path relative to this file (repo layout)
if not os.path.exists(BASE_DIR):
    alt = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "my_project_data"))
    if os.path.exists(alt):
        BASE_DIR = alt
    else:
        raise FileNotFoundError(
            f"Base data directory not found. Tried:\n - {BASE_DIR}\n - {alt}\n\n"
            "Set the BASE_DIR environment variable to the correct path or create the folder."
        )

print("Base path:", BASE_DIR)
print("Files:", os.listdir(BASE_DIR))

# Data subdirectories and files, relative to BASE_DIR
FCS_FOLDER = os.path.join(BASE_DIR, "fcs_files/")
BARCODES_PATH = os.path.join(BASE_DIR, "barcodes.csv")
MERGED_METADATA_PATH = os.path.join(BASE_DIR, "biorepo_2025_df_clean_REDCap_records_merged.csv")
DF_CLEAN_PATH = os.path.join(BASE_DIR, "df_clean_ready_for_web_app.csv") # Fixed path for the clean DFtput directory, also on the mounted GCS bucket

OUTPUT_DIR = os.path.join(GCS_MOUNT_DIR, "my_workflow_output")
COMBINED_ADATA_PATH = os.path.join(OUTPUT_DIR, "combined_adata_all_cells.h5ad")

# Sampling and other parameters
MAX_CELLS_TOTAL = 1000
CD8_MIN_CELLS = 100
SAMPLES_PER_TIMEPOINT = 25

# Valid timepoints
TIMEPOINTS = [
    'Pre-Treatment', 'On Treatment', 'Follow Up', 'Follow-Up1', 'Follow-Up2',
    'Irae', 'Irae 4 Week Follow-Up', 'Irae 6 Month Follow-Up', 'Irae 12 Month Follow-Up',
    'Eosi', 'End Of Study'
]