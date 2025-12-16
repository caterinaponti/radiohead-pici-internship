import pandas as pd
import anndata as ad
import numpy as np
from data_paths import MAX_CELLS_TOTAL, CD8_MIN_CELLS
from utils_io import parse_fcs
import logging

def load_adata(file_path, timepoint, patient_id, barcode, pop, pid):
    meta, data = parse_fcs(file_path)
    if data is None:
        raise ValueError("FCS parsing failed.")
    adata = ad.AnnData(data)
    adata.obs['Filename'] = file_path
    adata.obs['Timepoint'] = timepoint
    adata.obs['Patient ID'] = patient_id
    adata.obs['Barcode'] = barcode
    adata.obs['Barcode #'] = pop
    adata.obs['PID'] = pid
    return adata

def downsample_marker_aware(file_path, timepoint, patient_id, barcode, pop, pid):
    meta, data = parse_fcs(file_path)
    if data is None:
        raise ValueError("FCS parsing failed.")

    required_markers = ['CD8A', 'CD3']
    if not all(marker in data.columns for marker in required_markers):
        raise ValueError(f"Missing markers in: {file_path}")

    cd8_mask = (data['CD8A'] > 0.5) & (data['CD3'] > 0.5)
    cd8_cells = data[cd8_mask]
    non_cd8_cells = data[~cd8_mask]

    if len(cd8_cells) > CD8_MIN_CELLS:
        cd8_cells = cd8_cells.sample(n=CD8_MIN_CELLS, random_state=42)

    remaining = MAX_CELLS_TOTAL - len(cd8_cells)
    if remaining > 0 and len(non_cd8_cells) > remaining:
        non_cd8_cells = non_cd8_cells.sample(n=remaining, random_state=42)

    final_data = pd.concat([cd8_cells, non_cd8_cells])
    adata = ad.AnnData(final_data)
    adata.obs['Filename'] = file_path
    adata.obs['Timepoint'] = timepoint
    adata.obs['Patient ID'] = patient_id
    adata.obs['Barcode'] = barcode
    adata.obs['Barcode #'] = pop
    adata.obs['PID'] = pid
    return adata
