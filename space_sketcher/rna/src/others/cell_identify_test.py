#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmwrite
import gzip
import shutil
import scanpy as sc
from pathlib import Path
import anndata
from loguru import logger

def read_mtx_matrix(matrix_dir: Path, var_names="gene_ids"):
    """Read STARsolo raw matrix"""
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
    return adata

def write_10x_matrix(output_dir, adata):
    """Write AnnData object to 10x Genomics format - 优化版本"""    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 直接写入压缩的矩阵文件（避免中间文件）
    matrix_file = os.path.join(output_dir, 'matrix.mtx.gz')
    with gzip.open(matrix_file, 'wb') as f:
        # 使用内存缓冲避免磁盘IO
        import io
        buffer = io.BytesIO()
        mmwrite(buffer, adata.X.T)
        f.write(buffer.getvalue())
    
    # 2. 预准备特征数据（避免多次条件判断）
    if 'gene_ids' in adata.var.columns and 'gene_symbols' in adata.var.columns:
        gene_ids = adata.var['gene_ids'].values
        gene_symbols = adata.var['gene_symbols'].values
    elif 'gene_symbols' in adata.var.columns:
        gene_ids = adata.var_names.values
        gene_symbols = adata.var['gene_symbols'].values
    elif 'gene_ids' in adata.var.columns:
        gene_ids = adata.var['gene_ids'].values
        gene_symbols = adata.var_names.values
    else:
        gene_ids = adata.var_names.values
        gene_symbols = adata.var_names.values
    
    # 3. 使用numpy数组直接构建特征数据（避免DataFrame创建开销）
    features_file = os.path.join(output_dir, 'features.tsv.gz')
    with gzip.open(features_file, 'wt') as f:
        for i in range(adata.n_vars):
            f.write(f"{gene_ids[i]}\t{gene_symbols[i]}\tGene Expression\n")
    
    # 4. 直接写入barcodes（避免创建DataFrame）
    barcodes_file = os.path.join(output_dir, 'barcodes.tsv.gz')
    with gzip.open(barcodes_file, 'wt') as f:
        for barcode in adata.obs_names:
            f.write(f"{barcode}\n")


def umi_filter(adata, min_umi=200):
    """Simple binary threshold filter to replace emptyDrops"""
    # umi_counts = np.array(adata.X.sum(axis=1)).flatten()
    umi_counts = adata.X.sum(axis=1).A1
    # Binary threshold: keep barcodes with UMI >= min_umi
    selected_indices = np.where(umi_counts >= min_umi)[0]
    return selected_indices, umi_counts


def batch_replace(summary_data, mapping):
    """批量替换Summary文件中的值"""
    updated_data = []
    for line in summary_data:
        found = False
        for field, var_name in mapping.items():
            if line.startswith(field + ','):
                if var_name in globals():
                    new_value = globals()[var_name]
                    updated_data.append(f"{field},{new_value}")
                    found = True
                    break
        if not found:
            updated_data.append(line)
    return updated_data


def cell_identify(inputdir, outputdir, minumi=200, force_cells=0):
    """Main function to identify cells using binary threshold"""
    
    inputdir = "/data03/lead/userdata/huanghuichang/Work/Pipeline/Space-sketcher-dev/test/L3-7-10X-GM/01.count/Solo.out/GeneFull_Ex50pAS/raw"
    outputdir = "/data03/lead/userdata/huanghuichang/Work/Pipeline/Space-sketcher-dev/space_sketcher/rna/src/test"
    minumi=200
    force_cells=0
    # Remove existing output directory
    if os.path.exists(outputdir):
        shutil.rmtree(outputdir)
        print(f"Removed existing directory: {outputdir}")
    
    print(f"inputdir: {inputdir}\noutputdir: {outputdir}")
    print(f"min_umi: {minumi}\nforce_cells: {force_cells}")
    
    logger.info("Start reading matrix")
    adata = read_mtx_matrix(Path(inputdir))
    
    logger.info("Start cell filtering")
    if force_cells == 0:
        # Use binary threshold instead of emptyDrops
        selected_indices, umi_counts = umi_filter(adata, minumi)
        selected_barcodes = adata.obs_names[selected_indices]
        sorted_indices = np.argsort(umi_counts)[::-1]
    else:
        # Force cells mode
        umi_counts = adata.X.sum(axis=1).A1
        sorted_indices = np.argsort(umi_counts)[::-1]  # descending order
        selected_indices = sorted_indices[:min(force_cells, len(sorted_indices))]
        selected_barcodes = adata.obs_names.values[selected_indices]
    
    print(f"Number of filtered cell barcodes: {len(selected_indices)}")
    # Filter AnnData object and write output
    filtered_adata = adata[selected_indices, :]
    
    logger.info("Start writing matrix")
    write_10x_matrix(outputdir, filtered_adata)

    ###prepare plot data
    sorted_umi = umi_counts[sorted_indices]
    sorted_barcodes = adata.obs_names.values[sorted_indices]
    ranks = np.arange(1, len(sorted_barcodes) + 1)
    selected_dict = {barcode: True for barcode in selected_barcodes}
    is_cell_barcode = np.array([int(barcode in selected_dict) for barcode in sorted_barcodes], dtype=int)

    plotdata = pd.DataFrame({
        'barcode': sorted_barcodes,
        'UMI': sorted_umi,
        'is_cell_barcode': is_cell_barcode,
        'rank': ranks
    })
    # Output plot data
    outrank = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(inputdir))), 
                            'cell_rna_umi.rank.txt')
    plotdata.to_csv(outrank, sep='\t', index=False, float_format='%.0f')

    ### Process summary
    logger.info("Processing summary")
    summaryfile = os.path.join(os.path.dirname(inputdir), 'Summary.csv')
    cellreadsstat = os.path.join(os.path.dirname(inputdir), 'CellReads.stats')
    
    if os.path.exists(cellreadsstat):
        try:
            # Read cell reads stats
            cellreadsdata = pd.read_csv(cellreadsstat, sep='\t', header=0)
            
            # Check if required columns exist
            required_columns = ['CB', 'featureU', 'nUMIunique', 'nGenesUnique']
            if not all(col in cellreadsdata.columns for col in required_columns):
                print(f"Warning: CellReads.stats missing required columns. Expected: {required_columns}")
                return
            
            # Remove row that CB = CBnotInPasslist
            cellreadsdata = cellreadsdata[cellreadsdata['CB'] != 'CBnotInPasslist']
            truecellreads = cellreadsdata[cellreadsdata['CB'].isin(selected_barcodes)]
            
            # Calculate summary statistics
            estimated_number_of_cells = len(selected_barcodes)
            unique_reads_in_cells_mapped_to_gene = truecellreads['featureU'].sum()
            fraction_of_unique_reads_in_cells = round(unique_reads_in_cells_mapped_to_gene / cellreadsdata['featureU'].sum(), 4)
            mean_reads_per_cell = round(truecellreads['featureU'].mean(), 0)
            median_reads_per_cell = round(truecellreads['featureU'].median(), 0)
            umi_in_cells = truecellreads['nUMIunique'].sum()
            mean_umi_per_cell = round(truecellreads['nUMIunique'].mean(), 0)
            median_umi_per_cell = round(truecellreads['nUMIunique'].median(), 0)
            mean_gene_per_cell = round(truecellreads['nGenesUnique'].mean(), 0)
            median_gene_per_cell = round(truecellreads['nGenesUnique'].median(), 0)
            
            # Define field mapping
            field_var_mapping = {
                "Estimated Number of Cells": "estimated_number_of_cells",
                "Unique Reads in Cells Mapped to GeneFull_Ex50pAS": "unique_reads_in_cells_mapped_to_gene",
                "Fraction of Unique Reads in Cells": "fraction_of_unique_reads_in_cells",
                "Mean Reads per Cell": "mean_reads_per_cell",
                "Median Reads per Cell": "median_reads_per_cell",
                "UMIs in Cells": "umi_in_cells",
                "Mean UMI per Cell": "mean_umi_per_cell",
                "Median UMI per Cell": "median_umi_per_cell",
                "Mean GeneFull_Ex50pAS per Cell": "mean_gene_per_cell",
                "Median GeneFull_Ex50pAS per Cell": "median_gene_per_cell"
            }
            
            if os.path.exists(summaryfile):
                # Read and update summary file
                with open(summaryfile, 'r') as f:
                    summary_data = f.read().splitlines()
                
                updated_data = batch_replace(summary_data, field_var_mapping)
                
                outsummary = os.path.join(os.path.dirname(inputdir), 'Summary.callcell.csv')
                with open(outsummary, 'w') as f:
                    f.write('\n'.join(updated_data))
                
                print(f"Summary file updated successfully: {outsummary}")
            else:
                print(f"Summary file not found: {summaryfile}")
                
        except Exception as e:
            print(f"Error processing summary files: {e}")
            
    else:
        print(f"CellReads.stats file not found: {cellreadsstat}")
        print("Skipping summary processing")

def parse_args():
    parser = argparse.ArgumentParser(description='Filter cell barcodes using binary threshold')
    parser.add_argument('-i', '--inputdir', 
        metavar='PATH',
        type=str,
        required=True,
        help='Path of input (full path containing raw feature bc matrix)'
    )
    parser.add_argument('-o', '--outputdir', 
        metavar='PATH',
        type=str,
        help='Path of output, default filtered under inputdir'
    )
    parser.add_argument('-u', '--minumi', 
        metavar='INT',
        type=int,
        default=200,
        help='The minimum number of UMI for binary threshold [default: 200]'
    )
    parser.add_argument('-f', '--force_cells', 
        metavar='INT',
        type=int,
        default=0,
        help='Force cells [default: 0]'
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    
    # Set output directory
    if args.outputdir is None:
        args.outputdir = os.path.join(os.path.dirname(args.inputdir), 'callcell')
    
    # Call cell_identify function
    cell_identify(args.inputdir, args.outputdir, args.minumi, args.force_cells)