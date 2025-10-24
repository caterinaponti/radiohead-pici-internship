import scanpy as sc
from data_paths import COMBINED_ADATA_PATH

def preprocess_and_plot(cd8_adata, output_dir):
    markers_to_exclude = [
        'Residual', 'Time', 'beadDist', '131Xe_Environ', 'Center', 'Width',
        'Event_length', 'Offset', '138Ba_Environ', '120Sn_Environ', '133Cs_Environ'
    ]

    markers_to_keep = ~cd8_adata.var_names.isin(markers_to_exclude)
    cd8_adata = cd8_adata[:, markers_to_keep].copy()

    sc.pp.neighbors(cd8_adata, n_pcs=30, n_neighbors=15)
    sc.tl.leiden(cd8_adata, resolution=0.1)
    sc.tl.umap(cd8_adata, random_state=42)
    sc.pl.umap(cd8_adata, color=["leiden"], save="_leiden.png")

    # Marker visualization
    sc.pl.umap(cd8_adata, color=["CD3", "CD4", "CD8A"], cmap="viridis", vmin=0, vmax=2, save="_markers.png")


#vmc install to run jupyter notebook - GUI installer - display port VMC server 
# runnging GUI - forces to a particulr port andconnect the monitor to the cloud display
# need Wwindows App microsoft remote desktop - display port and location needs to talk to external IP for cloud - set up password 
