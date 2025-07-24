import os
import numpy as np
import time
import gzip
import polars as pl
import argparse
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline
from space_sketcher.tools.utils import judgeFilexits, csv2dict
from space_sketcher.__init__ import __root_dir__
import subprocess

"""
This function calculates various saturation statistics for scRNA-seq data,
such as mean reads per cell, median genes per cell, sequencing saturation,
and UMI saturation, at different sampling fractions.
"""

def prepare_readinfo(_indir, _threads, _lines:int):
    """
    Extract readinfo from bamfile
    """
    inbam = os.path.join(_indir,"Aligned.sortedByCoord.out.bam")
    judgeFilexits(inbam)

    software_version = subprocess.check_output(f"{__root_dir__}/software/sambamba --version", shell=True)
    print(f"sambamba 版本号：{software_version.decode('utf8')}")

    outtext = os.path.join(_indir, "readinfo.txt")
    extract_cmd = (
        f"{__root_dir__}/software/sambamba view {inbam} -t {_threads} |"
        f"head -n {_lines} | "
        f"""awk '{{print $1"\t"$27"\t"$20"\t"$28}}' > {outtext}"""
    )
    subprocess.check_call(extract_cmd, shell=True)


def load_filter_barcodes(filterbarcodefile):
    """
    Add prefix to filtered cell barcode
    """
    filterbarcodes = set()
    if filterbarcodefile.endswith("gz"):
        with gzip.open(filterbarcodefile, 'r') as f:
            for line in f:
                modifiedcb = "CB:Z:" + line.strip()
                filterbarcodes.add(modifiedcb)    
    else:
        with open(filterbarcodefile, 'r') as f:
            for line in f:
                modifiedcb = "CB:Z:" + line.strip()
                filterbarcodes.add(modifiedcb)
    return filterbarcodes


def downsample_and_calculate(infile, passcb):

    total_reads = int(os.popen(f"wc -l {infile}").read().split()[0])
    sampling_fractions = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    sampling_fractions_length = len(sampling_fractions)
    stats_df = pd.DataFrame(
        {"Sampling Fraction": np.zeros(sampling_fractions_length, np.float64),
        "Mean Reads per Cell": np.zeros(sampling_fractions_length, np.uint32),         
        "Median Genes per Cell": np.zeros(sampling_fractions_length, np.uint32),
        "Sequencing Saturation": np.zeros(sampling_fractions_length, np.float64),
        "Total mapped reads": np.zeros(sampling_fractions_length, np.uint32),},
        index=pd.Index(data=np.array(sampling_fractions), name="sampling_fraction"),)
    
    for fraction in sampling_fractions:
        if fraction == 0.0: continue
        sample_size = int(total_reads * fraction)
        sampled_indices = np.random.choice(total_reads, sample_size, replace=False)
        
        ###Load data with selected indices
        df_lazy = (
            pl.scan_csv(infile, has_header=False, separator="\t")
            .with_row_index("row_number")
            .filter(pl.col("row_number")
                    .is_in(sampled_indices))
                   )
        df = df_lazy.collect()
        del df_lazy
        df.columns = ["row_number", "readname", "Cell", "uGeneID", "UMI"]
        df = df.drop("row_number")

        ###Filter out rows with invalid values
        df = df.filter(pl.col("Cell")!="CB:Z:-").filter(pl.col("uGeneID")!="GX:Z:-")
        ###Calculate saturation
        yessubWLmatch_UniqueFeature = df.group_by("Cell").agg([pl.n_unique(["readname"]).alias("countedU")]).select([pl.col("countedU").sum()])["countedU"][0]
        yesUMIs = df.filter(pl.col("UMI")!="UB:Z:-").group_by("Cell").agg([pl.n_unique(['UMI']).alias("UMIcount")]).select([pl.col("UMIcount").sum()])["UMIcount"][0]
        Saturation = round((1- yesUMIs/yessubWLmatch_UniqueFeature)*100,3) if yessubWLmatch_UniqueFeature > 0 else 0

        ###Filter by passcb(Only keep cell barcodes that are in passcb)
        passdf = df.filter(pl.col("Cell").is_in(passcb))
        sampled_passcb = set(passdf["Cell"].to_list())
        mean_reads_per_cell = int(passdf.group_by("Cell").agg([pl.n_unique(["readname"]).alias("countedU")]).select([pl.col("countedU").sum()])["countedU"][0]/len(sampled_passcb))
        median_genes_per_cell = math.ceil(passdf.filter(pl.col("UMI")!="UB:Z:-").group_by("Cell").agg([pl.n_unique(["uGeneID"]).alias("uGenecount")]).select([pl.col("uGenecount").median()])["uGenecount"][0])
        del df
    
        stats_df.loc[fraction, "Sampling Fraction"] = fraction
        stats_df.loc[fraction, "Sequencing Saturation"] = Saturation
        stats_df.loc[fraction, "Mean Reads per Cell"] = mean_reads_per_cell
        stats_df.loc[fraction, "Median Genes per Cell"] = median_genes_per_cell
        stats_df.loc[fraction, "Total mapped reads"] = sample_size

    return stats_df

####Plot saturation
def to_percent(temp, position):
    return '%d'%(temp/1000) + 'k'

def umi_saturation(ax,table):

    ###check if duplicate number in Mean Reads per Cell, and replace with +1
    exits = []
    for i in range(0,len(table['Mean Reads per Cell'])):
        if table['Mean Reads per Cell'][i] in exits:
            table['Mean Reads per Cell'][i] = table['Mean Reads per Cell'][i] + 1
            exits.append(table['Mean Reads per Cell'][i])
        else:
            exits.append(table['Mean Reads per Cell'][i])
            
    xnew = np.linspace(table['Mean Reads per Cell'].min(),table['Mean Reads per Cell'].max(),300)
    smooth = make_interp_spline(table['Mean Reads per Cell'],table['Sequencing Saturation'])(xnew)
    ax.set_xlim([0, table['Mean Reads per Cell'].max()])
    ax.set_ylim([0, 0.9999])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(to_percent))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False)
    ax.grid(linestyle='--')
    ax.plot(xnew,smooth,linewidth=3.0)
    ax.axhline(y=0.9,ls="--",c="black",linewidth=2.0)
    ax.set(xlabel='Mean Reads per Cell', ylabel='Sequencing Saturation',title='Sequencing Saturation')
    return

def gene_saturation(ax,table):

    ###check if duplicate number in Mean Reads per Cell, and replace with +1
    exits = []
    for i in range(0,len(table['Mean Reads per Cell'])):
        if table['Mean Reads per Cell'][i] in exits:
            table['Mean Reads per Cell'][i] = table['Mean Reads per Cell'][i] + 1
            exits.append(table['Mean Reads per Cell'][i])
        else:
            exits.append(table['Mean Reads per Cell'][i])

    xnew = np.linspace(table['Mean Reads per Cell'].min(),table['Mean Reads per Cell'].max(),300)
    smooth = make_interp_spline(table['Mean Reads per Cell'],table['Median Genes per Cell'])(xnew)
    ax.set_xlim([0, table['Mean Reads per Cell'].max()])
    ax.set_ylim([0, table['Median Genes per Cell'].max()])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(to_percent))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False)
    ax.grid(linestyle='--')
    ax.plot(xnew,smooth,linewidth=3.0)
    ax.set(xlabel='Mean Reads per Cell', ylabel='Median Gene per Cell',title='Median Gene per Cell')

    return

def plot_saturation(inputfile, outdir):
    inputtable = pd.read_table(inputfile, sep="\t")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), tight_layout=True)
    umi_saturation(ax1,inputtable)
    gene_saturation(ax2,inputtable)
    fig.savefig(os.path.join(outdir, 'saturation.png'),facecolor='white',transparent=False,dpi=400)
    plt.close(fig)
    return


def count_saturation(indir, threads=4, lines=50000000):

    prepare_readinfo(indir, threads, lines)

    judgeFilexits(
        os.path.join(indir, "Solo.out/GeneFull_Ex50pAS/callcell/barcodes.tsv.gz"),
        os.path.join(indir, "Solo.out/GeneFull_Ex50pAS/Summary.callcell.csv"),
        os.path.join(indir, "readinfo.txt")
    )

    summary = csv2dict(
        os.path.join(indir, "Solo.out/GeneFull_Ex50pAS/Summary.callcell.csv")
    )

    filterbarcodes = load_filter_barcodes(os.path.join(indir, "Solo.out/GeneFull_Ex50pAS/callcell/barcodes.tsv.gz"))
    results = downsample_and_calculate(os.path.join(indir, "readinfo.txt"), filterbarcodes)
    one_ratio = results[results["Sampling Fraction"]==1]

    ###recalculate the sampling 
    final_res = pd.DataFrame()
    final_res["Sampling Fraction"] = results["Sampling Fraction"]
    final_res["Mean Reads per Cell"] = (int(summary["Mean Reads per Cell"])*(results["Mean Reads per Cell"]/int(one_ratio["Mean Reads per Cell"]))).astype(int)
    final_res["Median Genes per Cell"] = (int(summary["Median GeneFull_Ex50pAS per Cell"])*(results["Median Genes per Cell"]/int(one_ratio["Median Genes per Cell"]))).astype(int)
    final_res["Sequencing Saturation"] = round(float(summary["Sequencing Saturation"])*(results["Sequencing Saturation"]/int(one_ratio["Sequencing Saturation"])),3)

    resultsfile = os.path.join(indir,'saturation.xls')
    final_res.to_csv(resultsfile, sep="\t", index=False)
    plot_saturation(resultsfile, indir)


def parse_args():
    parser = argparse.ArgumentParser(description='sequencing saturation')
    parser.add_argument('--indir', 
        metavar='PATH', 
        type=str,
        help='input directory'
        )
    parser.add_argument(
        '--threads',
        metavar='INT',
        help='Analysis threads. [default: 4].',
        type=int,default=4
        )
    parser.add_argument(
        '--lines',
        metavar='INT',
        help='How many lines to calculate saturation',
        type=int,default=50000000
        )
    args = parser.parse_args()
    return args


if __name__=='__main__':
    args = parse_args()
    indir = args.indir
    threads = args.threads
    lines = args.lines
    
    count_saturation(indir, threads, lines)
