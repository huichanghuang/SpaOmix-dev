import pandas as pd 
import matplotlib.pyplot as plt
from Bio import Seq
from matplotlib.pyplot import MultipleLocator
import sys
import os

coordfile = sys.argv[1]
indir = os.path.abspath(sys.argv[2])

samplename = os.path.basename(indir)
cell_coord = os.path.join(indir, "02.oligo/cb_sb_coord.txt")
fastq_coord = os.path.join(indir, "02.oligo/spatial_umis.csv.gz")

inputdf = pd.read_csv(coordfile, header=None, names=["barcode", "xcoord", "ycoord"])
xmax = max(inputdf["xcoord"])
ymax = max(inputdf["ycoord"])
xmin = min(inputdf["xcoord"])
ymin = min(inputdf["ycoord"])

print(xmin, xmax, ymin, ymax)

celldf = pd.read_csv(cell_coord, header=0, sep="\t")
celldf = celldf[["sb", "xcoord", "ycoord"]]
celldf = celldf.drop_duplicates()
plt.figure(figsize=(int(xmax/1000),int(ymax/1000)))
plt.scatter(x=celldf["xcoord"], y = celldf["ycoord"],marker='o', s=0.05, alpha=1, rasterized=True)
x_major_locator = MultipleLocator(576.64)
y_major_locator = MultipleLocator(1089.68)
plt.xlim(0, xmax+100)
plt.ylim(0, ymax+100)
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.title(f"{samplename}.spatial_spot_in_cell")
plt.tight_layout()
plt.savefig(f"{indir}/{samplename}.spatial_spot_in_cell.png")


rawsampledf = pd.read_csv(fastq_coord, header=0,compression="gzip")
samplesbs = list(set(rawsampledf["Spatial_Barcode"]))
rc_samplesbs = [Seq.reverse_complement(i) for i in samplesbs if isinstance(i, str)]
sampledf = inputdf[inputdf["barcode"].isin(rc_samplesbs)]
sampledf = sampledf.drop_duplicates()

plt.figure(figsize=(int(xmax/1000),int(ymax/1000)))
plt.scatter(x=sampledf["xcoord"], y = sampledf["ycoord"],marker='o', s=0.05, alpha=1, rasterized=True)
x_major_locator = MultipleLocator(576.64)
y_major_locator = MultipleLocator(1089.68)
plt.xlim(0, xmax+100)
plt.ylim(0, ymax+100)
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.title(f"{samplename}.spatial_spot")
plt.tight_layout()
plt.savefig(f"{indir}/{samplename}.spatial_spot.png")
