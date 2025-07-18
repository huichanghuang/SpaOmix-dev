library(decontX)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(argparse)

rm(list = ls())
# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("--matrixdir", type = "character", help = "The directory of the matrix")
parser$add_argument("--samplename", type = "character", help = "The name of the sample")
parser$add_argument("--outdir", type = "character", help = "The output directory")
args <- parser$parse_args()



matrixdir <- args$matrixdir
samplename <- args$samplename
outdir <- args$outdir
# matrixdir <- "/data01/huanghuichang/scRNA_anlysis/20250109/02.count/0106-1229-WT-sampled/Solo.out/GeneFull_Ex50pAS/filtered"
# samplename<-"0106-1229-WT-sampled"
# outdir <- "/data01/home/huanghuichang/Test/09.decontX"
mincells <- 3
minfeatures <- 200
nvariables <- 5000
n_dims <- 30
reso <- 0.5

set.seed(1)  #设置随机数种子，使结果可重复
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


#####>>>>>Load Data And Quality Control<<<<<<#####
pbmc.data <- Read10X(data.dir = matrixdir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = samplename, min.cells = mincells, min.features = minfeatures)
pbmc$orig.ident <- samplename

#####>>>>>Normalization and Feature Selection<<<<<<#####
pbmc <- NormalizeData(pbmc)
# Identify the most highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nvariables)
# Scale the data
pbmc <- ScaleData(pbmc)

###Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
###Run non-linear dimensional reduction (UMAP)
pbmc <- FindNeighbors(pbmc, dims = 1:n_dims)
pbmc <- FindClusters(pbmc, resolution = reso)
pbmc <- RunUMAP(pbmc, dims = 1:n_dims, verbose = FALSE)

###decontX decontamination 
counts <- pbmc@assays$RNA$count
decontX_results <- decontX(counts)
#add contamination results to  metadata
pbmc$Contamination =decontX_results$contamination
head(pbmc@meta.data)

#plot contamination results
metadata <- pbmc@meta.data
outpic <- paste(outdir, paste(samplename, "RNA-Contamination.png", sep = "."), sep = "/")
png(file = outpic, width = 1600, height = 800, res=100)
p1 <- ggplot(metadata, aes(x = orig.ident, y = Contamination, fill = orig.ident))+
    geom_violin(scale = "width", trim = FALSE) +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplot inside
  theme_minimal() +  # Use a minimal theme
  labs(
    x = "Sample",  # X-axis label
    y = "Contamination",  # Y-axis label
    title = "Contamination Level by Sample"  # Plot title
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.line = element_line(color = "black", size = 0.5),  # Add axis lines
  )
p2 <- FeaturePlot(pbmc, features = c("Contamination"), cols = c("blue", "orange"))+
    scale_color_viridis_c()+
    ggtitle("RNA Contamination")
p1 + p2
dev.off()

# Calculate proportions
contamination_summary <- metadata %>%
  group_by(orig.ident) %>%  # Group by sample
  summarise(
    perc_gt_0.2 = mean(Contamination > 0.2) * 100,  # Percentage of cells with Contamination > 0.2
    perc_gt_0.5 = mean(Contamination > 0.5) * 100   # Percentage of cells with Contamination > 0.5
  )
# Convert to a matrix
contamination_matrix <- as.matrix(contamination_summary)  
colnames(contamination_matrix) <- c("sample", "Contamination>0.2", "Contamination>0.5")  

###output data
outfile1 <- paste(outdir, paste(samplename, "decontX-stat.txt", sep = "."), sep = "/")
outfile2 <- paste(outdir, paste(samplename, "decontX-metadata.txt", sep = "."), sep = "/")
write.table(contamination_matrix, file = outfile1, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(metadata, file = outfile2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
