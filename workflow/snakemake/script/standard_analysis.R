#!/usr/bin/env Rscript

" seurat 标准分析
Usage:
  analysis.R -i <file> -o <dir> [-s <species>] [--minFeature <int>] [--maxFeature <int> ] [--minCell <int>] [--pctMT <int>] [--npcs <int>] [--resolution <num>] [--method <str>] [-g <file>]
Options:
  -h --help                   Show this screen.
  -i,--input=file             配置文件
  -o,--outdir=dir             结果输出目录
  -s,--species=str            物种信息，可选 <human|mouse>  [default: human]
  --minFeature=int            每个细胞最小基因数 [default: 200]
  --maxFeature=int            每个细胞最大基因数 [default: Inf]
  --minCell=int               每个基因最小细胞数 [default: 3]
  --pctMT=int                 最大线粒体占比 [default: 20]
  --npcs=int                  使用的降维个数 [default: 30]
  --resolution=num            分辨率 [default: 0.8]
  --method=str                整合方法 [default: harmony]
  -g,--group=file             差异基因分组信息
" -> doc


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scCustomize)
library(scRNAtoolVis)
set.seed(322)


custom_colors <- c(ggsci::pal_d3("category20")(20), ggsci::pal_flatui()(10))

`%not_in%` <- function(x, table) {
  !(x %in% table)
}

validate_args <- function(args) {
  for (k in c("minFeature", "maxFeature", "minCell", "pctMT", "npcs", "resolution")) {
    args[[k]] <- suppressWarnings(as.numeric(args[[k]]))
    if (is.na(args[[k]]) || args[[k]] < 0) {
      stop(k, "值不应该为 NA 小于0")
    }
  }
  methods <- c("harmony", "scvi", "cca", "rpca")
  if (args$method %not_in% methods) {
    stop("method 必须为以下选项：", paste0(methods, collapse = ","))
  }
  args$input <- normalizePath(args$input, mustWork = TRUE)
  if (!is.null(args$group)) {
    args$group <- normalizePath(args$group, mustWork = TRUE)
  }
  args
}



#' 此函数用于 ide 交互式调试
get_args <- function() {
  args <- docopt::docopt(doc, args = c(
    "-i", "./config.csv",
    "-o", "./outdir",
    "-g", "./group.txt"
  ))
  args <- validate_args(args)
}


args <- docopt::docopt(doc)
args <- validate_args(args)

check_cols <- function(cols, config) {
  for (col in cols) {
    if (col %not_in% colnames(config)) {
      stop(sprintf("配置文件 %s 需要列 %s", args$input, col))
    }
  }
}

# 获取当前脚本目录
get_abs_dir <- function() {
  cmd <- commandArgs(F)
  file <- grep("--file=", cmd, value = TRUE)[1]
  if (is.na(file)) {
    path <- dirname(rstudioapi::getSourceEditorContext()$path)
  } else {
    file <- sub("--file=", "", file)
    path <- normalizePath(dirname(file))
  }
  path
}

source(file.path(get_abs_dir(), "utils.R"))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
setwd(args$outdir)

run_seurat <- function(args) {
  config <- data.table::fread(args$input, header = TRUE) |>
    unique()
  check_cols(c("sample", "path"), config)
  # to do: 是否需要先获取总细胞数，判断是否使用 BPcells

  obj_list <- apply(config, 1, function(x) {
    sample <- x["sample"]
    path <- x["path"]
    # 链接矩阵到当前交付目录
    upstream <- "01.upstream"
    dir <- file.path(upstream,sample)
    dir.create(dir,recursive = T,showWarnings = F)
    file.symlink(path,dir)

    counts <- read_sc_auto(path)

    obj <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
    obj$orig.ident <- as.character(sample)
    cols <- names(x[names(x) %not_in% c("sample", "path")])
    for (col in cols) {
      obj[[col]] <- as.character(x[col])
    }
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-|^GRCh38_MT-|^GRCm39_mt-")
    obj
  })

  # 是否需要整合
  is_integrated <- nrow(config) > 1


  if (is_integrated) {
    obj <- merge(obj_list[[1]], y = obj_list[-1])
  } else {
    obj <- obj_list[[1]]
  }
  rm(obj_list)
  invisible(gc())

  obj$orig.ident <- factor(obj$orig.ident, levels = unique(obj$orig.ident))
  qc_dir <- "02.qc"
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  plot_qc <- function(obj, prefix = paste0(qc_dir, "/raw_")) {
    for (feature in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
      p <- VlnPlot(obj, features = feature, group.by = "orig.ident", pt.size = 0.1)
      name <- gsub("_RNA|percent.", "", feature)
      ggsave(paste0(prefix, name, ".png"), p, width = 12, height = 8)
    }
    # 过滤前细胞数统计
    obj@meta.data |>
      group_by(orig.ident) |>
      summarise(
        cells = n(),
        median_gene = as.integer(median(nFeature_RNA)),
        median_umi = as.integer(median(nCount_RNA))
      ) |>
      data.table::fwrite(paste0(prefix, "qc_summary.csv"))
  }

  # 自动识别参数选项中的过滤选项
  data.frame(
    "过滤项" = c("基因过滤", "细胞过滤", "死细胞或破损细胞"),
    "阈值" = c(
      sprintf("过滤掉只在小于或等于 %s 个细胞中表达的基因", args$minCell),
      sprintf("细胞的表达基因数低于 %s、或高于 %s 时过滤", args$minFeature, args$maxFeature),
      sprintf("过滤掉线粒体基因比例大于 %s 的细胞", args$pctMT)
    )
  ) |>
    data.table::fwrite(file.path(qc_dir, "qc_options.csv"))

  plot_qc(obj)
  # 过滤
  obj <- filter_gene_merge(obj, min.cells = args$minCell)
  obj <- subset(obj, subset = nFeature_RNA > args$minFeature & nFeature_RNA < args$maxFeature & percent.mt < args$pctMT)
  plot_qc(obj, prefix = paste0(qc_dir, "/filtered_"))


  # 标准化
  cluster_dir <- "03.cluster"
  dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 50)

  p <- ElbowPlot(obj, ndims = 50) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    )
  ggsave(file.path(cluster_dir, "elbow_plot.png"), p, width = 8, height = 6)
  p <- DimPlot(obj, reduction = "pca", group.by = "orig.ident")
  ggsave(file.path(cluster_dir, "pca_plot.png"), p, width = 8, height = 6)


  # 整合
  # 聚类使用的分辨率
  res_used <- unique(c(seq(0.1, 2, by = 0.3), args$resolution), fromLast = TRUE)
  if (is_integrated) {
    obj <- RunUMAP(obj, dims = 1:args$npcs, reduction = "pca", reduction.name = "unintegrated")
    # methord
    mlist <- switch(args$method,
      "harmony" = list(method = HarmonyIntegration, new.reduction = "harmony"),
      "cca" = list(method = CCAIntegration, new.reduction = "integrated.cca"),
      "rpca" = list(method = RPCAIntegration, new.reduction = "integrated.rpca"),
      "scvi" = list(method = scVIIntegration, new.reduction = "scvi")
    )
    obj <- IntegrateLayers(
      object = obj, method = mlist$method,
      orig.reduction = "pca", new.reduction = mlist$new.reduction,
      verbose = TRUE
    )
    obj <- RunUMAP(obj, reduction = mlist$new.reduction, dims = 1:args$npcs, reduction.name = "umap")
    p1 <- DimPlot(obj, group.by = "orig.ident", reduction = "unintegrated") +
      labs(title = "unintegrated") +
      theme(plot.title = element_text(hjust = 0.5))
    p2 <- DimPlot(obj, group.by = "orig.ident", reduction = "umap") +
      labs(title = "integrated") +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(cluster_dir, "integrated.png"), p1 + p2, width = 12, height = 6)
    # 聚类
    obj <- FindNeighbors(obj, reduction = mlist$new.reduction, dims = 1:args$npcs)
    obj <- FindClusters(obj, resolution = res_used)
  } else {
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:args$npcs, reduction.name = "umap")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:args$npcs)
    obj <- FindClusters(obj, resolution = res_used)
  }
  obj <- JoinLayers(obj)

  plot_clustree(obj, filename = file.path(cluster_dir, "clustree.png"))

  use_colors <- rep(custom_colors, len = length(unique(obj$seurat_clusters)))
  # 聚类图
  p1 <- DimPlot(obj, group.by = "orig.ident", reduction = "umap") +
    labs(title = "samples") +
    theme(plot.title = element_text(hjust = 0.5))

  p2 <- DimPlot(obj, group.by = "seurat_clusters", reduction = "umap") +
    labs(title = "seurat_clusters") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = use_colors)

  ggsave(file.path(cluster_dir, "cluster.png"), p1 + p2, width = 12, height = 6)

  # 样本分开
  p1 <- DimPlot(obj, split.by = "orig.ident", reduction = "umap", group.by = "seurat_clusters", ncol = 6) +
    labs(title = "seurat_clusters") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = use_colors)

  ggsave(file.path(cluster_dir, "cluster_by_sample.png"), p1, width = 12, height = 6)
  # 类型分布图
  p <- cellRatioPlot(
    object = obj,
    sample.name = "orig.ident",
    celltype.name = "seurat_clusters",
    col.width = 0.7,
    fill.col = use_colors
  )

  ggsave(file.path(cluster_dir, "cell_ratio.png"), p, width = 12, height = 6)

  # marker 基因
  marker_dir <- "04.marker_gene"
  dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

  marker_genes <- FindAllMarkers(obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25) |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC)) |>
    ungroup() |>
    dplyr::relocate(cluster, gene)

  data.table::fwrite(marker_genes, file.path(marker_dir, "all_markers.csv"))

  n_genes <- 2
  top_markers <- marker_genes |>
    group_by(cluster) |>
    dplyr::filter(avg_log2FC > 0.25 & p_val_adj < 0.01) |>
    top_n(n = n_genes,wt = avg_log2FC) |>
    ungroup()

  # 气泡图
  p <- DotPlot(obj, features = unique(top_markers$gene)) +
    Seurat::RotatedAxis() +
    scale_colour_gradient2(low = "#2480b8", mid = "#E1E1DF", high = "#d5626a") +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white")
    )
  ggsave(file.path(marker_dir, "DotPlot.png"), width = 10, height = 6)


  p <- DoHeatmap(obj, features = top_markers$gene, size = 3, angle = 90, group.colors = use_colors) + NoLegend()
  ggsave(file.path(marker_dir, "heatmap.png"), plot = p, width = 10, height = 12)

  for (cluster_name in unique(top_markers$cluster)) {
    df <- top_markers |>
      dplyr::filter(cluster == cluster_name)
    genes <- df$gene
    p <- FeaturePlot(obj, features = genes, ncol = 3, reduction = "umap", combine = FALSE, order = TRUE)
    p <- wrap_plots(p, ncol = 3) + plot_layout(axis_titles = "collect") +
      plot_annotation(title = paste0("cluster_", cluster_name))
    filename <- paste0(marker_dir, "/cluster_", cluster_name, "_feature.png")
    ggsave(filename, p, width = n_genes * 5, height = 5)
  }

  p <- Stacked_VlnPlot(obj,
    features = unique(top_markers$gene), x_lab_rotate = TRUE,
    colors_use = use_colors, split.by = "seurat_clusters", plot_spacing = 0.3
  )

  ggsave(file.path(marker_dir, "Stacked_VlnPlot.png"), width = 10, height = 20)

  LayerData(obj, layer = "scale.data") <- NULL
  saveRDS(obj, "cluster.rds")
  obj
}

obj <- run_seurat(args)

cat("单细胞 meta data 如下：\n")
dplyr::glimpse(obj@meta.data)

# 差异基因
if (!is.null(args$group)) {
  rlang::check_installed("presto", "希望安装 presto 包加速")
  deg_dir <- "05.deg"
  dir.create(deg_dir, showWarnings = F, recursive = T)
  group <- data.table::fread(args$group, header = TRUE)
  if (any(duplicated(group$name))) {
    warning(args$group, "文件 name 列包含相同的名称，去除重复项")
    group <- group[!duplicated(group$name), ]
  }

  colname <- c("name", "colname", "group1", "group2")
  if (!all(colname %in% colnames(group))) {
    stop(args$group, "文件必须包含列名：", paste0(colname, collapse = ","))
  }

  apply(group, 1, function(row) {
    group_1 <- strsplit(row["group1"],";")[[1]]
    group_2 <- strsplit(row["group2"],";")[[1]]
    deg <- FindMarkers(obj, ident.1 = group_1, ident.2 = group_2, group.by = as.character(row["colname"]))
    deg$gene <- rownames(deg)
    if (nrow(deg) == 0) {
      warning("定义的分组文件中，行", row["name"], "未找到差异基因")
      return(NULL)
    } else {
      p <- EnhancedVolcano::EnhancedVolcano(deg, x = "avg_log2FC", y = "p_val", lab = rownames(deg), subtitle = row["name"], FCcutoff = 1, pCutoff = 0.05, pointSize = 3, labSize = 3)
      ggsave(file.path(deg_dir, paste0(row["name"], "_volcano.png")), p, width = 10, height = 12)
      data.table::fwrite(deg, file.path(deg_dir, paste0(row["name"], ".csv")))
    }
  })
}
