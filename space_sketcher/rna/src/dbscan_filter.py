#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from joblib import Parallel, delayed
from tqdm import tqdm  # 进度条支持
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
import scipy.io
import random
import matplotlib.colors as mcolors
norm = mcolors.Normalize(vmin=0, vmax=3)


def read_data(infile):
    """Read input file and return DataFrame"""
    return pd.read_csv(infile, sep='\t')

def write_10x_format(outdir, adata, gene_symbols=True):
    """
    输出10X Genomics标准格式（Seurat兼容）
    
    参数:
        outdir: 输出目录路径
        adata: AnnData对象
        gene_symbols: 是否使用基因符号作为feature名称
    """
    os.makedirs(outdir, exist_ok=True)
    
    # 1. 写入barcodes.tsv
    barcodes_path = os.path.join(outdir, "barcodes.tsv")
    pd.DataFrame(adata.obs_names).to_csv(
        barcodes_path, 
        sep='\t', 
        header=False, 
        index=False
    )
    
    # 2. 写入features.tsv
    features_path = os.path.join(outdir, "genes.tsv")
    if gene_symbols and 'gene_symbols' in adata.var.columns:
        features_df = pd.DataFrame({
            'gene_ids': adata.var.index,
            'gene_symbols': adata.var['gene_symbols'],
            'feature_type': 'Gene Expression'
        })
    else:
        features_df = pd.DataFrame({
            'gene_ids': adata.var.index,
            'gene_symbols': adata.var.index,  # 如果没有gene_symbols则用gene_id
            'feature_type': 'Gene Expression'
        })
    features_df.to_csv(
        features_path, 
        sep='\t', 
        header=False, 
        index=False
    )
    
    # 3. 写入matrix.mtx (稀疏矩阵格式)
    matrix_path = os.path.join(outdir, "matrix.mtx")
    scipy.io.mmwrite(
        matrix_path,
        csr_matrix(adata.X.T)  # 注意需要转置
    )
    
    print(f"10X format files saved to: {outdir}")

def filter_by_umi(df, maxumi=5000, minumi=2):
    """Filter spatial barcodes by total UMI counts"""
    subsb_umi_summary = df.groupby('subsb')['umi_count'].sum().reset_index(name='umi_count_sum')
    subsb_umi_summary = subsb_umi_summary.sort_values('umi_count_sum', ascending=False)
    subsb_umi_summary['rank'] = range(1, len(subsb_umi_summary)+1)
    
    subsb_to_filter = subsb_umi_summary[
        (subsb_umi_summary['umi_count_sum'] > maxumi) | 
        (subsb_umi_summary['umi_count_sum'] < minumi)
    ]['subsb']
    
    return df[~df['subsb'].isin(subsb_to_filter)], subsb_umi_summary

def process_single_cell(barcode, cell_data, eps=150, min_samples=6):
    """
    处理单个细胞barcode的DBSCAN聚类
    参数:
        barcode: 当前处理的细胞barcode
        cell_data: 该细胞的所有空间坐标和UMI数据
        eps: DBSCAN半径参数
        min_samples: DBSCAN最小样本数
    返回:
        (barcode, centroid_coords, cluster_type)
    """
    # 数据准备
    coords = cell_data[['xcoord', 'ycoord']].values
    weights = np.sqrt(cell_data['umi_count'].values)
    
    # 跳过太小的簇
    if len(coords) < min_samples:
        return (barcode, None, "no-cluster")
    
    # DBSCAN聚类
    db = DBSCAN(
        eps=eps,
        min_samples=min_samples,
        metric='euclidean',
        algorithm='kd_tree',  # 更快的算法
        n_jobs=1  # 每个任务单线程
    ).fit(coords, sample_weight=weights)
    
    labels = db.labels_
    valid_labels = labels[labels != -1]  # 去除噪声点
    
    # 无有效簇的情况
    if len(valid_labels) == 0:
        return (barcode, None, "no-cluster")
    
    # 统计簇信息
    unique_labels, counts = np.unique(valid_labels, return_counts=True)
    main_cluster = unique_labels[np.argmax(counts)]
    main_mask = (labels == main_cluster)
    
    # 计算加权质心
    x_coords = cell_data['xcoord'].values[main_mask]
    y_coords = cell_data['ycoord'].values[main_mask]
    cluster_weights = cell_data['umi_count'].values[main_mask]
    
    centroid_x = round(np.sum(x_coords * cluster_weights) / np.sum(cluster_weights))
    centroid_y = round(np.sum(y_coords * cluster_weights) / np.sum(cluster_weights))
    
    # 判断簇类型
    if len(counts) > 1 and np.sort(counts)[-2] > max(counts)/2:
        cluster_type = "multi-cluster"
    else:
        cluster_type = "single-cluster"
    
    return (barcode, (centroid_x, centroid_y), cluster_type)

def perform_dbscan_parallel(df, eps=150, min_samples=6, n_jobs=4, batch_size=100):
    """
    并行DBSCAN处理
    参数:
        df: 包含所有细胞数据的DataFrame
        eps: DBSCAN半径参数
        min_samples: DBSCAN最小样本数 
        n_jobs: 并行任务数(-1表示使用所有核心)
        batch_size: 每个任务处理的细胞数
    返回:
        filtered_df, cb_cluster, cluster_stats
    """
    # 预处理：按细胞barcode分组
    grouped = df.groupby('cb')
    barcodes = list(grouped.groups.keys())
    
    # 并行处理
    results = Parallel(n_jobs=n_jobs, batch_size=batch_size)(
        delayed(process_single_cell)(barcode, grouped.get_group(barcode), eps, min_samples)
        for barcode in tqdm(barcodes, desc="Processing cells")
    )
    
    ## 有效细胞(有簇的)
    valid_results = [r for r in results if r[1] is not None]
    coord_data = {
        'cb': [r[0] for r in valid_results],
        'xcoord': [r[1][0] for r in valid_results],
        'ycoord': [r[1][1] for r in valid_results]
    }
    coord_df = pd.DataFrame(coord_data)
    
    ## 所有细胞的聚类类型
    cluster_data = {
        'cb': [r[0] for r in results],
        'cluster': [r[2] for r in results]
    }
    cb_cluster = pd.DataFrame(cluster_data)
    
    ## 聚类统计
    cluster_stats = cb_cluster['cluster'].value_counts().reset_index()
    cluster_stats.columns = ['cluster_type', 'counts']
    cluster_stats['ratio'] = np.round(
        cluster_stats['counts'] / cluster_stats['counts'].sum(), 
        4
    )
    
    return coord_df, cb_cluster, cluster_stats

def generate_plots(df, cb_cluster, subsb_umi_summary, outdir):
    """Generate diagnostic plots"""
    # Spatial barcode UMI knee plot
    plt.figure(figsize=(8, 6))
    plt.plot(subsb_umi_summary['rank'], subsb_umi_summary['umi_count_sum'], 'o-', markersize=3)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("UMI Count Sum of Each Spatial Barcode")
    plt.xlabel("Rank")
    plt.ylabel("UMI Count Sum")
    plt.savefig(os.path.join(outdir, "spatial_umi_knee_plot.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Cell barcode mean UMI plot
    df = df.merge(cb_cluster, on='cb', how='left')
    
    cb_umi_mean_top100 = df.groupby('cb').apply(
        lambda x: x.nlargest(100, 'umi_count')['umi_count'].mean()
    ).reset_index(name='umi_count_mean')
    cb_umi_mean_top100 = cb_umi_mean_top100.merge(cb_cluster, on='cb', how='left')
    cb_umi_mean_top100 = cb_umi_mean_top100.sort_values('umi_count_mean', ascending=False)
    cb_umi_mean_top100['rank'] = range(1, len(cb_umi_mean_top100)+1)
    

    colordict = {"single-cluster":"lightskyblue",
                 "multi-cluster": "lightsalmon",
                 "no-cluster": "limegreen"}
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=cb_umi_mean_top100, 
                    x='rank', y='umi_count_mean', 
                    hue='cluster', linewidth = 0,
                    palette=colordict,
                    legend=False)
    sns.lineplot(data=cb_umi_mean_top100, 
                    x='rank', y='umi_count_mean', 
                    hue='cluster', markers=True,
                    linewidth=2,
                    palette=colordict)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("Mean of UMI Count of Top 100 Spatial Barcodes in Each Cell")
    plt.xlabel("Rank")
    plt.ylabel("Mean of UMI count")
    plt.savefig(os.path.join(outdir, "cb_umi_mean_knee_plot.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Cell barcode vs spatial barcodes count plot
    sb_count = df.groupby(['cb', 'cluster'])['sb'].nunique().reset_index(name='sbcounts')
    sb_count = sb_count.sort_values('sbcounts', ascending=False)
    sb_count['rank'] = range(1, len(sb_count)+1)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=sb_count, 
                    x='rank', y='sbcounts', 
                    hue='cluster', linewidth = 0,
                    palette=colordict,
                    legend=False)
    sns.lineplot(data=sb_count, 
                    x='rank', y='sbcounts', 
                    hue='cluster', markers=True,
                    linewidth=2,
                    palette=colordict)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("Spatial Barcode Counts of Each Cell")
    plt.xlabel("Rank")
    plt.ylabel("sb counts")
    plt.savefig(os.path.join(outdir, "cb_sb_counts.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    return cb_umi_mean_top100

def filter_matrix(matrixdir, coord_df, outdir):
    """Filter matrix based on DBSCAN results"""
    adata = sc.read_10x_mtx(matrixdir, var_names='gene_symbols')

    adata_barcodes = adata.obs_names.values
    cell_barcodes = coord_df['cb'].tolist()
    # 验证barcode交集
    valid_barcodes = set(cell_barcodes) & set(adata_barcodes)
    print(f"Found {len(valid_barcodes)} matching barcodes (originally {len(cell_barcodes)})")

    # 创建安全的布尔索引
    mask = np.isin(adata_barcodes, list(valid_barcodes))
    # 安全过滤
    filtered_adata = adata[mask, :].copy()
    # 输出前检查
    assert len(filtered_adata) == len(valid_barcodes), "过滤后维度不匹配"    

    outmtxdir = os.path.join(outdir, "filtered_matrix")
    if os.path.exists(outmtxdir):
        import shutil
        shutil.rmtree(outmtxdir)
        
    write_10x_format(outmtxdir, filtered_adata)
    coord_df.to_csv(os.path.join(outmtxdir, "spatial_location_information.txt"), sep='\t', index=False)

    
def generate_summary(cellreads, coord_df, cb_umi_mean_top100, cluster_stats, outdir):
    """Generate summary statistics after dbscan filter"""
    cellreadsdata = pd.read_csv(cellreads, sep='\t')
    cellreadsdata = cellreadsdata[cellreadsdata['CB'] != "CBnotInPasslist"]

    cell_barcodes = coord_df['cb'].tolist()
    truecellreads = cellreadsdata[cellreadsdata['CB'].isin(cell_barcodes)]
    
    stats = {
        'Cells with Certain Location': len(truecellreads),
        'Unique Reads in Cells Mapped to Gene': truecellreads['featureU'].sum(),
        'Fraction of Unique Reads in Cells': round(truecellreads['featureU'].sum()/cellreadsdata['featureU'].sum(), 4),
        'Mean Reads per Cell': round(truecellreads['featureU'].mean(), 0),
        'Median Reads per Cell': round(truecellreads['featureU'].median(), 0),
        'UMIs in Cells': truecellreads['nUMIunique'].sum(),
        'Fraction of UMIs in Cells': round(truecellreads['nUMIunique'].sum()/cellreadsdata['nUMIunique'].sum(), 4),
        'Mean UMI per Cell': round(truecellreads['nUMIunique'].mean(), 0),
        'Median UMI per Cell': round(truecellreads['nUMIunique'].median(), 0),
        'Mean Gene per Cell': round(truecellreads['nGenesUnique'].mean(), 0),
        'Median Gene per Cell': round(truecellreads['nGenesUnique'].median(), 0),
        'Median top100 Spatial UMI Mean per Cell': round(cb_umi_mean_top100['umi_count_mean'].median(), 3),
        'Mean top100 Spatial UMI Mean per Cell': round(cb_umi_mean_top100['umi_count_mean'].mean(), 3)
    }
    
    stats.update(dict(zip(cluster_stats['cluster_type'], cluster_stats['ratio'])))
    
    with open(os.path.join(outdir, "dbscan_filtered_cells.summary.csv"), 'w') as f:
        for k, v in stats.items():
            f.write(f"{k},{v}\n")


def sb_umi_distribution(df, cbs, cluster_type, outdir):
    fig = plt.figure(figsize=(24,30))
    for i, cb in enumerate(cbs):
        df_cb = df[df["cb"] == cb]
        ax = fig.add_subplot(4, 5, i+1, projection='3d')
        ax.scatter(df_cb["xcoord"], df_cb["ycoord"], df_cb["umi_count"], c=df_cb["umi_count"], cmap='Blues', norm=norm)
        # ax.set_xlim(0, 10000)
        # ax.set_ylim(0, 10000)
        ax.set_title(cb, fontsize=10)
        ax.set_xlabel('X', fontsize=10)
        ax.set_ylabel('Y', fontsize=10)
        ax.set_zlabel('UMI_count', fontsize=10)
    plt.suptitle(cluster_type, fontsize=16)
    plt.tight_layout()
    outpic = os.path.join(outdir, cluster_type+".png")
    plt.savefig(outpic)
    print("Finish ploting "+outpic)
    return


def dbscan_filter(infile, outdir, maxumi, minumi, matrixdir, cellreads, eps, min_samples, n_jobs):
    os.makedirs(outdir, exist_ok=True)
    # Process data
    df = read_data(infile)
    filtered_df, subsb_umi_summary = filter_by_umi(df, maxumi, minumi)
    coord_df, cb_cluster, cluster_stats = perform_dbscan_parallel(filtered_df, eps, min_samples, n_jobs)
    
    # Save intermediate results
    cluster_stats.to_csv(os.path.join(outdir, "dbscan_cluster_distribution.csv"), index=False)
    cb_cluster.to_csv(os.path.join(outdir, "CB_cluster.txt"), sep='\t', index=False)
    subsb_umi_summary.to_csv(os.path.join(outdir, "spatial_umi_knee_plot.temp.txt"), sep='\t', index=False)
    
    # Generate plots and additional data
    cb_umi_mean_top100 = generate_plots(filtered_df, cb_cluster, subsb_umi_summary, outdir)
    cb_umi_mean_top100.to_csv(os.path.join(outdir, "cb_umi_mean_knee_plot.temp.txt"), sep='\t', index=False)

    # Filter matrix and generate final outputs
    filter_matrix(matrixdir, coord_df, outdir)
    generate_summary(cellreads, coord_df, cb_umi_mean_top100, cluster_stats, outdir)

    ###Plot sb-umi distribution of each cluster type
    singlecluster_cb, multicluster_cb, nocluster_cb = [], [], []
    if "single-cluster" in cb_cluster["cluster"].unique():
        singlecluster_cb = cb_cluster[cb_cluster["cluster"] == "single-cluster"]["cb"].tolist()
    if "multi-cluster" in cb_cluster["cluster"].unique():
        multicluster_cb = cb_cluster[cb_cluster["cluster"] == "multi-cluster"]["cb"].tolist()
    if "no-cluster" in cb_cluster["cluster"].unique():
        nocluster_cb = cb_cluster[cb_cluster["cluster"] == "no-cluster"]["cb"].tolist()
    
    if len(singlecluster_cb) > 20:
        singlecluster_cb = random.sample(singlecluster_cb, 20)
    if len(multicluster_cb) > 20:
        multicluster_cb = random.sample(multicluster_cb, 20)
    if len(nocluster_cb) > 20:
        nocluster_cb = random.sample(nocluster_cb, 20)

    cluster_cbs = [singlecluster_cb, multicluster_cb, nocluster_cb]
    cluster_types = ["single-cluster", "multi-cluster", "no-cluster"]
    for cbs, cluster_type in zip(cluster_cbs, cluster_types):
        if len(cbs) == 0:
            continue
        sb_umi_distribution(df, cbs, cluster_type, outdir)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Perform DBSCAN clustering for cell barcodes and filter matrix.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--infile',
        metavar='FILE',
        type=str,
        required=True,
        help='Input coordinate file containing cb, sb, xcoord, ycoord, umi_count'
        )
    parser.add_argument('-o', '--outdir',
        metavar='PATH',
        type=str,
        required=True,
        help='Output directory for results'
        )
    parser.add_argument('-m', '--matrixdir',
        metavar='PATH',
        type=str,
        required=True,
        help='Directory containing 10X format matrix'
        )
    parser.add_argument('-cr', '--cellreads',
        metavar='FILE',
        type=str,
        required=True,
        help='Cell reads file from STARsolo output'
        )
    parser.add_argument('-c', '--chemistry',
        metavar='STR',
        type=str,
        choices=['10X', 'leader_v1'],
        default='leader_v1',
        help='Chemistry version: 10X or leader_v1, [default: leader_v1]'
        )
    parser.add_argument('--maxumi',
        metavar='INT',
        type=int,
        default=5000,
        help='Maximum UMI count threshold for spatial barcodes filtering'
        )
    parser.add_argument('--minumi',
        metavar='INT',
        type=int,
        default=2,
        help='Minimum UMI count threshold for spatial barcodes filtering'
        )
    parser.add_argument('--eps',
        metavar='FLOAT',
        type=float,
        default=150.0,
        help='DBSCAN epsilon parameter (maximum distance between points)'
        )
    parser.add_argument('--min_samples',
        metavar='INT',
        type=int,
        default=6,
        help='DBSCAN min_samples parameter (minimum points to form cluster)'
        )
    parser.add_argument('--n_jobs',
        metavar='INT',
        type=int,
        default=-1,
        help='Number of CPU cores to use (-1 for all available)'
        )
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    dbscan_filter(args.infile, args.outdir, args.maxumi, 
                  args.minumi, args.matrixdir, args.cellreads, 
                  args.eps, args.min_samples, args.n_jobs)
