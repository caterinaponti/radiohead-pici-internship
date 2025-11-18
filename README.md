# My Cloud Workflow

This folder contains Python scripts and modules for processing, analyzing, and managing flow cytometry data in a cloud workflow. For cloud environment setup, see `CLOUD_SETUP.md`.

## Folder Structure
- `analysis.py` — Main analysis functions and workflows.
- `data_paths.py` — Centralized paths and configuration for data locations.
- `data_processing.py` — Data cleaning, transformation, and preprocessing utilities.
- `metadata_processing.py` — Functions for handling and processing metadata.
- `utils_io.py` — Input/output helper functions for reading and writing data.
- `main.py` — Example entry point for running the workflow.
- `requirements.txt` — Python dependencies for the workflow.
- `CLOUD_SETUP.md` — Step-by-step guide for setting up a cloud environment, GUI, and storage.
- `__pycache__/` — Compiled Python files (auto-generated).

## Quick Start
1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
2. **Configure data paths**:
   - Edit `data_paths.py` to set the correct paths to your data files and output directories.

3. **Run the workflow**:
   - You can run the main workflow or individual scripts as needed. For example:
   ```bash
   python main.py
   ```

## Requirements
- Python 3.7+
- pandas
- numpy
- (other packages as listed in `requirements.txt`)

## Cloud Setup
For full instructions on setting up a cloud VM, GUI, Python environments, and Google Cloud Storage mounting, see [`CLOUD_SETUP.md`](CLOUD_SETUP.md).

