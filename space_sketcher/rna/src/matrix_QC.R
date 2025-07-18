library(dplyr)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(argparse)
library(glmGamPoi)
library(cowplot)


# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("--matrixdir", type = "character", help = "The directory of the matrix")
parser$add_argument("--outdir", type = "character", help = "The output directory")
parser$add_argument("--mincells", type = "numeric", help = "The minimum number of cells", default = 3)
parser$add_argument("--minfeatures", type = "numeric", help = "The minimum number of features(for CreateSeuratObject filter)", default = 5)
parser$add_argument("--nvariables", type = "numeric", help = "The number of features for variable selection", default = 2000)
parser$add_argument("--ndims", type = "numeric", help = "The number of dimensions for PCA", default = 30)
parser$add_argument("--resolution", type = "numeric", help = "The resolution for clustering", default = 0.5)
args <- parser$parse_args()

matrixdir <- args$matrixdir
outdir <- args$outdir
mincells <- args$mincells
minfeatures <- args$minfeatures
nvariables <- args$nvariables
n_dims <- args$ndims
reso <- args$resolution


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


#####>>>>>Load Data And Quality Control as Spatial<<<<<<#####
###load the matrix
COUNTS_MTX <- Read10X(data.dir = matrixdir)
slide.seq = CreateSeuratObject(counts = COUNTS_MTX, assay="Spatial", min.cells = mincells, min.features = minfeatures)
##four patterns: ^MT-, ^mt-, ^GRCh38_MT-, ^GRCm39_mt-
slide.seq[["percent.mt"]] <- PercentageFeatureSet(slide.seq, pattern = "^MT-|^mt-|^GRCh38_MT-|^GRCm39_mt-|^GRCh38-MT-|^GRCm39-mt-")


###cellbarcode after filtering
filtered_spots <- colnames(slide.seq)
####>>>>>Add the spatial information<<<<<<#####
coord.file <- paste0(matrixdir, "/spatial_location_information.txt")
coord.df = read.table(coord.file, header = T, sep = "\t")
###deduplicate the barcode, select the first one
#check if the cb has duplicated
print(any(duplicated(coord.df$cb)))
coord.df <- coord.df[!duplicated(coord.df$cb),]
dim(coord.df)
###only keep xcoord, ycoord
coord.df <- coord.df[,c("cb", "xcoord", "ycoord")]
rownames(coord.df) <- coord.df$cb
select.df <- coord.df[coord.df$cb %in% filtered_spots,]
###remove the cb column
select.df <- select.df[,-1]
print("The number of spots finally selected: ")
print(dim(select.df)[1])

####>>>>>Add the spatial information <<<<<<#####
slide.seq@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = select.df
  )


####>>>>>normalize the data using sctransform and perform a standard scRNA-seq dimensionality reduction and clustering workflow.<<<<<<#####
# slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 2000, verbose = FALSE)
slide.seq <- NormalizeData(slide.seq, normalization.method = "LogNormalize", scale.factor = 10000)
slide.seq <- FindVariableFeatures(slide.seq, selection.method = "vst", nfeatures = nvariables)
all.genes <- rownames(slide.seq)
slide.seq <- ScaleData(slide.seq, features = all.genes)

slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:n_dims)
slide.seq <- FindNeighbors(slide.seq, reduction = "pca", dims = 1:n_dims)
slide.seq <- FindClusters(slide.seq, resolution = reso, verbose = FALSE)


#####>>>>>prepare plot data<<<<<<#####
coordplot <- SpatialFeaturePlot(slide.seq, features = "nCount_Spatial")
coorddata <- coordplot$data
umapplot <- FeaturePlot(slide.seq, features = c("nCount_Spatial"))
umapdata <- umapplot$data
umapdata$cells <- rownames(umapdata)
mergedf <- merge(coorddata, umapdata, by = "cells")
mergedf <- mergedf[,c("cells", "x", "y", "nCount_Spatial.x", "umap_1", "umap_2", "ident")]
colnames(mergedf) <- c("CB", "xcoord", "ycoord", "nCount_Spatial", "UMAP1", "UMAP2", "Cluster")
head(mergedf)
#                 CB xcoord ycoord nCount_Spatial     UMAP1      UMAP2 Cluster
# 1 AACTAGACAGCGACTG   2735   2249            730  2.343697  -3.703772       0
# 2 AACTAGACATCAGTAG   3872   2092             13 -4.195808 -12.848030       3
# 3 AACTAGACGCGCGACA   3193    997            370  2.409972  -4.225308       0
metadata <- slide.seq@meta.data
metadata <- metadata[c("nFeature_Spatial", "percent.mt")]
metadata$CB <- rownames(metadata)
mergedf <- merge(mergedf, metadata, by = "CB")
mergedf <- mergedf[,c("CB",  "nCount_Spatial","nFeature_Spatial", "percent.mt", "xcoord", "ycoord","UMAP1", "UMAP2", "Cluster")]
outumap <- paste(outdir, "UMAPpos.txt", sep = "/")
write.table(mergedf, outumap, sep = "\t", quote = F, row.names = F)


#####>>>>>plot by UMI counts <<<<<<#####
p1 <- ggplot(mergedf, aes(x = xcoord, y = ycoord, fill = nCount_Spatial)) +
    geom_point(shape = 21, size = 2, alpha = 1, stroke=0) +
    scale_fill_distiller(palette = "Spectral") +
    theme_void()+
    ggtitle("Spots colored by UMI counts")

p2 <- ggplot(mergedf, aes(x = UMAP1, y = UMAP2, color = nCount_Spatial)) +
    geom_point(size = 2) +
    scale_color_distiller(palette = "Spectral") +
    theme_minimal()+
    theme(plot.margin = margin(10, 10, 10, 10))+
    ggtitle("UMAP projection of spots colored by UMI counts")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")


outpic1 <- paste(outdir, "UMAP_umicounts.png", sep = "/")
png(file = outpic1, width = 1600, height = 800, res=100)
plot_grid(p1, p2, rel_widths = c(1.5, 1.2))
dev.off()

#####>>>>>plot by cluster <<<<<<#####
p3 <- ggplot(mergedf, aes(x = xcoord, y = ycoord, fill = Cluster)) +
    geom_point(shape = 21, size = 2, alpha = 1, stroke=0) +
    scale_fill_manual(values = my36colors) +
    theme_void()+
    theme(plot.margin = margin(10, 10, 10, 10))+
    ggtitle("Spots Colored by Cluster")

p4 <- ggplot(mergedf, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 2) +
    scale_color_manual(values = my36colors) +
    theme_minimal()+
    theme(plot.margin = margin(10, 10, 10, 10))+
    ggtitle("UMAP projection of spots colored by Cluster")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")

outpic2 <- paste(outdir, "UMAP_cluster.png", sep = "/")
png(file = outpic2, width = 1600, height = 800, res=100)
plot_grid(p3, p4, rel_widths = c(1.5, 1.2))
dev.off()


#####>>>>>Extract diff genes<<<<<<#####
marker_genes <- FindAllMarkers(slide.seq, only.pos = TRUE, min.pct = 0.25, log2FC.threshold = 0.25, test.use = "wilcox")

dim(marker_genes)
if (dim(marker_genes)[1] == 0){
    print("No DEG detected!")
}else{
    marker_genes %>%
        group_by(cluster) %>%
        top_n(30, avg_log2FC) -> topmarkers
    head(topmarkers)
    ###rearrange the data, move gene names to the first column, and cluster names to the second column
    topmarkers <- topmarkers %>%
        select(gene, cluster, everything())
    head(topmarkers)
    outdiff <- paste(outdir, "top30DEG.txt", sep = "/")
    write.table(topmarkers, outdiff, sep = "\t", quote = F, row.names = F)
}

###print all the parameter values used in the analysis
print(paste("mincells:", mincells))
print(paste("nvariables:", nvariables))
print(paste("n_dims:", n_dims))
print(paste("resolution:", reso))

###save RDS file
LayerData(slide.seq, "scale.data") <- NULL
saveRDS(slide.seq, paste0(matrixdir, "/Cluster.rds"))