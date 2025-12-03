import scanpy as sc
import numpy as np
import os
from data_paths import COMBINED_ADATA_PATH
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples
import plotly.express as px 


def save_subsampled_adata(adata, output_dir, n_cells=1000):
    # ensure obs_names are unique
    if not adata.obs_names.is_unique:
        import pandas as pd
        adata.obs_names = (
            adata.obs_names.astype(str) + "_" +
            pd.Series(range(adata.shape[0]), index=adata.obs_names).astype(str)
        )
    if adata.n_obs > n_cells:
        idx = np.random.choice(adata.obs_names, n_cells, replace=False)
        adata_sub = adata[idx, :].copy()
    else:
        adata_sub = adata.copy()
    sub_path = os.path.join(output_dir, "cd8_adata_subsampled.h5ad")
    adata_sub.write(sub_path)
    return adata_sub, sub_path

def plot_heatmap(adata, marker_cols, output_dir):
    heatmap_dir = os.path.join(output_dir, "Heatmaps")
    try:
        os.makedirs(heatmap_dir, exist_ok=True)
    except Exception as e:
        print(f"Warning: Could not create directory {heatmap_dir}: {e}")
        return
    import scanpy as sc
    sc.settings.figdir = heatmap_dir

    markers_to_exclude = ['Residual', 'Time', 'beadDist', '131Xe_Environ', 'Center', 'Width', 'Event_length', 'Offset',  '138Ba_Environ', '120Sn_Environ', '133Cs_Environ', 'Event_length' ]
    marker_cols = [v for v in marker_cols if v not in markers_to_exclude]
    try:
        sc.pl.matrixplot(
            adata,
            var_names=marker_cols,
            groupby="leiden",
            cmap="viridis",
            save="_cd8catalystheatmaps.png",
            show=False
        )
    except Exception as e:
        print(f"Warning: Could not save heatmap: {e}")
    

def preprocess_and_plot(cd8_adata, output_dir):
    sc.settings.figdir = os.path.join(output_dir, "CD8_UMAPs")
    os.makedirs(sc.settings.figdir, exist_ok=True)
    markers_to_exclude = [
        'Residual', 'Time', 'beadDist', '131Xe_Environ', 'Center', 'Width',
        'Event_length', 'Offset', '138Ba_Environ', '120Sn_Environ', '133Cs_Environ'
    ]

    markers_to_keep = ~cd8_adata.var_names.isin(markers_to_exclude)
    cd8_adata = cd8_adata[:, markers_to_keep].copy()

    sc.pp.neighbors(cd8_adata, n_pcs=30, n_neighbors=15)
    sc.tl.leiden(cd8_adata, resolution=0.1)
    sc.tl.umap(cd8_adata, random_state=42)
    sc.pl.umap(cd8_adata, color=["leiden"], save="_leiden.png", show=False)

    # marker visualization
    sc.pl.umap(cd8_adata, color=["CD3", "CD4", "CD8A"], cmap="viridis", vmin=0, vmax=2, save="_markers.png", show=False)

    # save AnnData object
    cd8_adata.write(os.path.join(output_dir, "cd8_adata_full.h5ad"))
    # subsample and save
    marker_cols = [v for v in cd8_adata.var_names if v not in markers_to_exclude]
    adata_sub, sub_path = save_subsampled_adata(cd8_adata, output_dir)
    # heatmap
    plot_heatmap(cd8_adata, marker_cols, output_dir)

def umap_prep(adata, n_pcs=30, n_neighbors=15, leiden_res=0.2, random_state=42):
    import scanpy as sc
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.leiden(adata, resolution=leiden_res)
    sc.tl.umap(adata, random_state=random_state)
    return adata

def pre_and_on_treatment_umap_workflow(merged_adata, output_dir):
    import scanpy as sc
    markers = ["leiden", "CD3", "CD8A", "CD4", "FOXP3", "CD14", "CD19", "CD56"]
    # pre-Treatment
    pre = merged_adata[merged_adata.obs['Timepoint'] == 'Pre-Treatment'].copy()
    print(f"Pre-Treatment shape: {pre.shape}")
    pre = umap_prep(pre)
    try:
        sc.pl.umap(pre, color=markers, save="_pre_treatment_umap.png", show=False)
    except Exception as e:
        print(f"Failed to save pre-treatment UMAP plot: {e}")
    
    # on-Treatment
    on = merged_adata[merged_adata.obs['Timepoint'] == 'On Treatment'].copy()
    print(f"On-Treatment shape: {on.shape}")
    on = umap_prep(on)
    try:
        sc.pl.umap(on, color=markers, save="_on_treatment_umap.png", show=False)
    except Exception as e:
        print(f"Failed to save on-treatment UMAP plot: {e}")
   
    return pre, on

import scanpy as sc
import os

def plot_markers_by_timepoint(pre_adata, on_adata, output_dir, prefix=""):
    """
    Plot UMAPs colored by clustering markers for pre-treatment and on-treatment samples.
    Saves PNGs for each marker in a subdirectory.
    """
    clustering_markers = [
        'CD45RA', 'CCR7', 'CD27', 'CD28', 'CD127', 'CD95',
        'PD1', 'CTLA4', 'LAG3', 'TIM3', 'TIGIT', 'CD69', 'ICOS', 'CD25'
    ]
    umap_dir = os.path.join(output_dir, "MarkerUMAPs")
    os.makedirs(umap_dir, exist_ok=True)

    # Helper to plot for one AnnData object
    def plot_umaps(adata, group_label):
        for marker in clustering_markers:
            if marker in adata.var_names:
                try:
                    sc.pl.umap(
                        adata,
                        color=[marker],
                        cmap="viridis",
                        vmin=0,
                        vmax=5,
                        save=f"_{prefix}{group_label}_{marker}.png",
                        show=False
                    )
                    print(f"Saved UMAP for {marker} ({group_label})")
                except Exception as e:
                    print(f"Could not plot {marker} for {group_label}: {e}")
            else:
                print(f"Marker {marker} not found in {group_label} AnnData.")

    # Set Scanpy's figure directory
    sc.settings.figdir = umap_dir

    plot_umaps(pre_adata, "Pre-Treatment")
    plot_umaps(on_adata, "On-Treatment")

    print(f"All marker UMAPs saved to: {umap_dir}") 

def save_silhouette_plots(adata, output_dir, cluster_range=[5, 6, 7, 8, 9, 10], umap_key='X_umap', prefix=""):
    """
    Compute and save silhouette plots for KMeans clustering on UMAP coordinates.
    """
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_samples, silhouette_score
    import numpy as np
     
    X = adata.obsm[umap_key]
    sil_dir = os.path.join(output_dir, "SilhouettePlots")
    os.makedirs(sil_dir, exist_ok=True)
    for n_clusters in cluster_range:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)
        ax1.set_xlim([-0.1, 1])
        ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(X)
        silhouette_avg = silhouette_score(X, cluster_labels)
        sample_silhouette_values = silhouette_samples(X, cluster_labels)
        y_lower = 10
        for i in range(n_clusters):
            ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
            ith_cluster_silhouette_values.sort()
            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i
            color = cm.nipy_spectral(float(i) / n_clusters)
            ax1.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                ith_cluster_silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7,
            )
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
            y_lower = y_upper + 10
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
        ax1.set_yticks([])
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
        colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(X[:, 0], X[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k")
        centers = clusterer.cluster_centers_
        ax2.scatter(centers[:, 0], centers[:, 1], marker="o", c="white", alpha=1, s=200, edgecolor="k")
        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")
        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("UMAP1")
        ax2.set_ylabel("UMAP2")
        plt.suptitle(f"Silhouette analysis for KMeans clustering (n_clusters = {n_clusters})", fontsize=14, fontweight="bold")
        fname = f"{prefix}silhouette_kmeans_{n_clusters}.png"
        plt.savefig(os.path.join(sil_dir, fname))
        plt.close(fig)

def plot_silhouette_scores_as_boxplot(adata, output_dir, cluster_range, umap_key='X_umap', prefix=""):
    """
    Compute silhouette scores for all samples across a range of K values 
    and display them as an interactive boxplot grouped by K using Plotly.
    """
    X = adata.obsm.get(umap_key)
    if X is None:
        print(f"Warning: UMAP coordinates not found in adata.obsm['{umap_key}']")
        return

    all_scores = []
    print("Calculating individual silhouette scores for boxplot...")
    for n_clusters in cluster_range:
        clusterer = KMeans(n_clusters=n_clusters, random_state=10, n_init='auto')
        cluster_labels = clusterer.fit_predict(X)
        sample_silhouette_values = silhouette_samples(X, cluster_labels)
        for score in sample_silhouette_values:
            all_scores.append({'K': n_clusters, 'Silhouette Score': score})

    df_scores = pd.DataFrame(all_scores)

    fig = px.box(
        df_scores, 
        x="K", 
        y="Silhouette Score", 
        color="K",
        title="Distribution of Silhouette Scores by Number of Clusters (K)",
        points="outliers"
    )
    fig.update_traces(quartilemethod="exclusive")
    fig.update_layout(
        xaxis_title="Number of Clusters (K)",
        yaxis_title="Individual Sample Silhouette Score",
        yaxis_range=[-0.2, 1.0]
    )

    boxplot_dir = os.path.join(output_dir, "ClusterEvaluationPlots")
    os.makedirs(boxplot_dir, exist_ok=True)
    fname = f"{prefix}silhouette_boxplot_K_range.html"
    fig.write_html(os.path.join(boxplot_dir, fname))
    print(f"Interactive boxplot saved to: {os.path.join(boxplot_dir, fname)}")


def plot_silhouette_boxplot_per_cluster(adata, output_dir, k=5, umap_key='X_umap', prefix=""):
    """
    Generate a boxplot of silhouette scores per cluster for a given K.
    """
    X = adata.obsm.get(umap_key)
    if X is None:
        print(f"Warning: UMAP coordinates not found in adata.obsm['{umap_key}']")
        return

    # KMeans clustering
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    cluster_labels = kmeans.fit_predict(X)

    # Silhouette scores
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    # DataFrame for plotting
    silhouette_df = pd.DataFrame({
        'Silhouette Score': sample_silhouette_values,
        'Cluster': cluster_labels
    })

    # Plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Cluster', y='Silhouette Score', data=silhouette_df)
    plt.title(f'Silhouette Score Distribution per Cluster for K={k}')
    plt.xlabel('Cluster Label')
    plt.ylabel('Silhouette Score')
    plt.grid(axis='y', linestyle='--')

    # Save
    boxplot_dir = os.path.join(output_dir, "ClusterEvaluationPlots")
    os.makedirs(boxplot_dir, exist_ok=True)
    fname = f"{prefix}silhouette_boxplot_per_cluster_k{k}.png"
    plt.savefig(os.path.join(boxplot_dir, fname), bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Per-cluster silhouette boxplot saved to: {os.path.join(boxplot_dir, fname)}")