import os
import logging
import fcsparser

def setup_logging(log_file="processing.log"):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        filemode="w"
    )
    logging.getLogger().addHandler(logging.StreamHandler())

def list_fcs_files(folder):
    return [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".fcs")]
    
def parse_fcs(file_path):
    try:
        return fcsparser.parse(file_path)
    except Exception as e:
        logging.error(f"Failed to parse {file_path}: {e}")
        return None, None
