import pandas as pd
import os
from data_paths import TIMEPOINTS
from utils_io import parse_fcs

def build_filename_dict(fcs_files):
    import re
    filename_dict = {}
    for file_path in fcs_files:
        filename = os.path.basename(file_path)
        parts = filename.split('_')
        pid = None
        pop = None
        if len(parts) >= 5:
            for part in parts:
                if len(part) == 3:
                    pid = part
                if part.startswith('Pop'):
                    pop = part.split('.')[0].lower()
            if pid and pop:
                _, data = parse_fcs(file_path)
                num_cells = data.shape[0] if data is not None else 0
                filename_dict[file_path] = (pid, pop, num_cells)
    return filename_dict

def collect_filename_metadata(merged_df, barcodes_df, fcs_barcodes, filename_dict):
    metadata_rows = []
    for tp in TIMEPOINTS:
        tp_df = merged_df[merged_df['Barcode Timepoint'] == tp]
        tp_barcodes = pd.unique(tp_df['Barcode'])
        tp_fcs = set(fcs_barcodes) & set(tp_barcodes)
        tp_dict = {}

        for barcode in tp_fcs:
            row = barcodes_df[barcodes_df['Client Sample ID'] == barcode]
            if not row.empty:
                pid = row['PID'].values[0]
                pop = row['Barcode #'].values[0]
                tp_dict[barcode] = (pid, pop)

        for filepath, (pid, pop, _) in filename_dict.items():
            for barcode, (b_pid, b_pop) in tp_dict.items():
                if (pid, pop) == (b_pid, b_pop):
                    match = merged_df[merged_df['Barcode'] == barcode]
                    participant_id = match['Participant ID'].values[0] if not match.empty else None
                    metadata_rows.append({
                        'Filename': filepath,
                        'Timepoint': tp,
                        'Patient ID': participant_id,
                        'Barcode': barcode,
                        'Barcode #': b_pop,
                        'PID': b_pid
                    })
                    break
    return pd.DataFrame(metadata_rows)
