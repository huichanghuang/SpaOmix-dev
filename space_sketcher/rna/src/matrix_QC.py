#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
# from scipy.sparse import csr_matrix
from pathlib import Path
import anndata
from space_sketcher.tools.utils import add_log
import numpy as np

# Custom color palette
MY36COLORS = [
    '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
    '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
    '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
    '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
    '#968175'
]

@add_log
def read_mtx_matrix(matrix_dir: Path, var_names="gene_symbols"):
    """Read STARsolo output matrix with spatial coordinates."""
    matrix_dir = Path(matrix_dir)
    suffix = ".gz" if (matrix_dir / "matrix.mtx.gz").exists() else ""
    
    adata = sc.read(matrix_dir / f"matrix.mtx{suffix}").T
    adata.X = adata.X.tocsr()
    
    # Read genes/features
    genes_file = matrix_dir / f"features.tsv{suffix}" if (matrix_dir / f"features.tsv{suffix}").exists() else matrix_dir / f"genes.tsv{suffix}"
    genes = pd.read_csv(genes_file, header=None, sep="\t")
    
    if var_names == "gene_symbols":
        adata.var_names = anndata.utils.make_index_unique(pd.Index(genes[1].values))
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    
    # Read barcodes
    barcodes = pd.read_csv(matrix_dir / f"barcodes.tsv{suffix}", header=None)
    adata.obs_names = barcodes[0].values
    # Add spatial coordinates to AnnData object
    coord_df = pd.read_csv(matrix_dir / f"spatial_location_information.txt{suffix}", sep='\t', header=0, index_col=0)
    ###sort coord_df name by matrix names
    coord_df = coord_df.reindex(list(adata.obs_names))

    adata.obsm['spatial'] = coord_df[['xcoord', 'ycoord']].values    
    return adata

def perform_qc_filtering(adata, min_features=5, min_cells=3):
    """Perform quality control filtering."""
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-', 'GRCh38_MT-', 'GRCm39_mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)    
    # Plot QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, show=False)    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_features)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    if adata.n_obs < 50:
        raise ValueError("Filtered cell count <50. Adjust 'pct_mt' or 'min_features'.")
    
    return adata

def process_expression_data(adata, n_top_genes=2000, n_pcs=30, resolution=0.5):
    """Process expression data through normalization, clustering, etc."""
    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)    
    # Feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
    adata = adata[:, adata.var.highly_variable]    
    # Dimensionality reduction
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)    
    # Clustering
    sc.tl.louvain(adata, resolution=resolution, flavor="igraph")
    
    return adata

@add_log
def generate_plots(adata, outdir):
    """Generate spatial and UMAP plots."""
    plot_data = pd.DataFrame({
        'CB': adata.obs_names,
        'xcoord': adata.obsm['spatial'][:, 0],
        'ycoord': adata.obsm['spatial'][:, 1],
        'nCount_Spatial': adata.obs['total_counts'],
        'UMAP1': adata.obsm['X_umap'][:, 0],
        'UMAP2': adata.obsm['X_umap'][:, 1],
        'Cluster': adata.obs['louvain'],
        'nFeature_Spatial': adata.obs['n_genes_by_counts'],
        'percent.mt': adata.obs['pct_counts_mt']
    })
    
    # UMI Counts Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    scatter1 = ax1.scatter(plot_data['xcoord'], plot_data['ycoord'],
                          c=plot_data['nCount_Spatial'], cmap='Spectral_r',
                          s=10, alpha=1, edgecolors='none')
    ax1.set_title('Spots colored by UMI counts', pad=20)
    ax1.axis('off')
    
    scatter2 = ax2.scatter(plot_data['UMAP1'], plot_data['UMAP2'],
                          c=plot_data['nCount_Spatial'], cmap='Spectral_r', s=5)
    ax2.set_title('UMAP projection', pad=20)
    ax2.axis('off')
    plt.colorbar(scatter2, ax=ax2, shrink=0.5).set_label("nCount_Spatial")
    plt.savefig(f"{outdir}/UMAP_umicounts.png", dpi=100, bbox_inches='tight')
    plt.close()
    
    # Cluster Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    for i, cluster in enumerate(sorted(plot_data['Cluster'].unique())):
        cluster_data = plot_data[plot_data['Cluster'] == cluster]
        ax1.scatter(cluster_data['xcoord'], cluster_data['ycoord'],
                   color=MY36COLORS[i % len(MY36COLORS)], label=cluster,
                   s=10, alpha=1, edgecolors='none')
        ax2.scatter(cluster_data['UMAP1'], cluster_data['UMAP2'],
                   color=MY36COLORS[i % len(MY36COLORS)], label=cluster, s=5)
    
    ax1.set_title('Spots Colored by Cluster', pad=20)
    ax1.axis('off')
    ax2.set_title('UMAP projection', pad=20)
    ax2.axis('off')
    ax2.legend(bbox_to_anchor=(1, 1), title="Cluster")
    plt.savefig(f"{outdir}/UMAP_cluster.png", dpi=100, bbox_inches='tight')
    plt.close()
    
    return plot_data

@add_log
def find_markers(adata, outdir):
    """Find differentially expressed genes."""
    sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', pts=True, n_genes=30)
    # 提取所有相关信息
    marker_genes = sc.get.rank_genes_groups_df(
        adata, group=None, key="rank_genes_groups"
    )
    marker_genes.columns = ["cluster", "gene", "scores", "log2FC", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"]
    if not marker_genes.empty:
        marker_genes = marker_genes.groupby('cluster').head(30)
        marker_genes.to_csv(f"{outdir}/top30DEG.txt", sep='\t', index=False)
    else:
        print("No DEG detected!")


def perform_matrix_qc(matrixdir, outdir, minfeatures, mincells,
                    nvariables, ndims, resolution):
    
    os.makedirs(outdir, exist_ok=True)
    # 1. Load data
    adata = read_mtx_matrix(matrixdir)
    # 2. Quality control
    adata = perform_qc_filtering(adata, minfeatures, mincells)    
    # 3. Process data
    adata = process_expression_data(adata, nvariables, ndims, resolution)        
    # 4. Find markers
    find_markers(adata, outdir)
    # 5. Generate plots and outputs
    plot_data = generate_plots(adata, outdir)
    plot_data.to_csv(f"{outdir}/UMAPpos.txt", sep='\t', index=False)
    # 6. Save results
    adata.write(f"{matrixdir}/cluster.h5ad")
    del adata
    # Print parameters
    print(f"Parameters used:\nmincells: {mincells}\nminfeatures: {minfeatures}")
    print(f"nvariables: {nvariables}\nndims: {ndims}\nresolution: {resolution}")

def parse_args():
    parser = argparse.ArgumentParser(description='Process spatial transcriptomics matrix data.')
    parser.add_argument(
        '-m', '--matrixdir', 
        required=True, 
        help='Directory containing the matrix (mtx/files) and coordfile')
    parser.add_argument(
        '-o', '--outdir', 
        required=True,
        help='Output directory')
    parser.add_argument(
        '-c','--mincells', 
        type=int, 
        default=3,
        help='Minimum cells per gene')
    parser.add_argument(
        '-f','--minfeatures', 
        type=int, 
        default=5,
        help='Minimum features per cell')
    parser.add_argument(
        '-v', '--nvariables', 
        type=int, 
        default=2000,
        help='Number of variable genes')
    parser.add_argument(
        '-d','--ndims', 
        type=int, 
        default=30,
        help='PCA dimensions')
    parser.add_argument(
        '-r','--resolution', 
        type=float, 
        default=0.5,
        help='Clustering resolution')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    perform_matrix_qc(args.matrixdir, args.outdir, args.minfeatures, args.mincells,
                    args.nvariables, args.ndims, args.resolution)
