import pandas as pd
import anndata as ad
from data_paths import FCS_FOLDER, BARCODES_PATH, MERGED_METADATA_PATH, OUTPUT_DIR, SAMPLES_PER_TIMEPOINT
from utils_io import setup_logging, list_fcs_files
from metadata_processing import build_filename_dict, collect_filename_metadata
from data_processing import downsample_marker_aware
from analysis import preprocess_and_plot
import os

barcodes = pd.read_csv(BARCODES_PATH)
merged_df = pd.read_csv(MERGED_METADATA_PATH)
fcs_files = list_fcs_files(FCS_FOLDER)

def main():
    setup_logging()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    barcodes = pd.read_csv(BARCODES_PATH)
    merged_df = pd.read_csv(MERGED_METADATA_PATH)
    fcs_files = list_fcs_files(FCS_FOLDER)
    filename_dict = build_filename_dict(fcs_files)
    fcs_barcodes = pd.unique(barcodes['Client Sample ID'])

    metadata_df = collect_filename_metadata(merged_df, barcodes, fcs_barcodes, filename_dict)
    metadata_df.to_csv(os.path.join(OUTPUT_DIR, "metadata_rows.csv"), index=False)
    print("Collected metadata.")

    filtered_metadata = metadata_df.groupby("Timepoint").apply(
        lambda x: x.sample(n=min(SAMPLES_PER_TIMEPOINT, len(x)), random_state=42)
    ).reset_index(drop=True)
    print("Filtered metadata.")

    
    adata_list = []

    for _, row in filtered_metadata.iterrows():
        file_path = row['Filename']
        timepoint = row['Timepoint']
        patient_id = row['Patient ID']
        barcode = row['Barcode']
        pop = row['Barcode #']
        pid = row['PID']

        try:
            adata = downsample_marker_aware(file_path, timepoint, patient_id, barcode, pop, pid)
            adata.obs['sample_id'] = file_path  
            adata_list.append(adata)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    print("Downsampling done.")

    for adata in adata_list:
        for col in adata.obs.select_dtypes(['category']).columns:
            adata.obs[col] = adata.obs[col].astype(str)

    # Ensure keys are unique for ad.concat, skip duplicates
    keys = []
    unique_adata_list = []
    seen = set()
    for a in adata_list:
        key = a.obs['sample_id'][0]
        if key in seen:
            print(f"Warning: Duplicate sample_id found and skipped: {key}")
            continue
        seen.add(key)
        keys.append(key)
        unique_adata_list.append(a)

    merged_adata = ad.concat(
        unique_adata_list,
        join="outer",
        axis=0,
        label="sample_id",
        keys=keys
    )

    print("Starting clustering and UMAP visualization...")
    preprocess_and_plot(merged_adata, OUTPUT_DIR)
    print("Clustering and plotting completed successfully!")
    print(f"Plots saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()