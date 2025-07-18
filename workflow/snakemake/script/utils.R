#' h5ad, rds 文件读取
read_h5ad_rds <- function(path) {
  if (grepl(".h5ad$", path)) {
    obj <- schard::h5ad2seurat(path, use.raw = TRUE, load.obsm = FALSE)
  } else if (grepl(".rds$", path)) {
    obj <- readRDS(obj)
  }
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  if (ncol(counts) == 0) {
    stop("counts 为空")
  }
}

#' csv 等文件读取为 counts
read_sc_txt <- function(file) {
  suppressMessages(require(Matrix))
  counts <- data.table::fread(file, header = "auto")
  # 基因名去重
  counts[[1]] <- make.unique(counts[[1]])
  counts <- as.matrix(counts, rownames = 1)
  counts <- as(counts, "dgCMatrix")
  counts
}

#' 简单的自动识别文件格式读取单细胞数据
read_sc_auto <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path)) {
    features <- file.path(path, "features.tsv")
    if (file.exists(features)) {
      counts <- ReadSTARsolo(path)
    } else {
      features_df <- data.table::fread(file.path(path, "features.tsv.gz"),sep = "\t",nrows = 10)
      gene.column <- if (ncol(features_df) == 1) 1 else 2
      counts <- Read10X(path,gene.column = gene.column)
    }
  } else {
    if (grepl("(csv|tsv|txt)(.gz)?$", path)) {
      counts <- read_sc_txt(path)
    } else if (grepl(".h5$", path)) {
      counts <- Read10X_h5(path)
    } else if (grepl("(.h5ad|.rds)$", path)) {
      counts <- read_h5ad_rds(path)
    } else {
      stop(path, "未知的扩展名，扩展名需为 rds,h5ad,csv.txt,tsv 等")
    }
  }
  counts
}


#' 绘制小提琴图
getvlnplot <- function(obj) {
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  vln <- VlnPlot(obj, features = feats, ncol = length(feats), cols = "#8491B4FF", combine = F)
  p1 <- vln[[1]] + ylab("Number")
  p1[["layers"]][[2]][["aes_params"]][["colour"]] <- "LightSkyBlue4"
  p1[["layers"]][[2]][["aes_params"]][["alpha"]] <- 0.5
  p2 <- vln[[2]] + ylab("Number")
  p2[["layers"]][[1]][["aes_params"]][["fill"]] <- "LightGoldenrod1"
  p2[["layers"]][[2]][["aes_params"]][["colour"]] <- "LightGoldenrod3"
  p2[["layers"]][[2]][["aes_params"]][["alpha"]] <- 0.5
  p3 <- vln[[3]] + ylab("Percentage(%)")
  p3[["layers"]][[2]][["aes_params"]][["colour"]] <- "DarkOliveGreen4"
  p3[["layers"]][[2]][["aes_params"]][["alpha"]] <- 0.5
  p <- p1 + p2 + p3 + plot_layout(ncol = length(feats)) & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 15)) & guides(fill = "none")
  p
}


filter_gene <- function(obj, min.cells = 0) {
  counts <- LayerData(obj, layer = "counts")
  if (inherits(x = counts, what = "IterableMatrix")) {
    rlang::check_installed(pkg = "BPCells", reason = "for working with BPCells")
    row_stat <- BPCells::matrix_stats(
      matrix = counts,
      row_stats = "nonzero"
    )$row_stats
    features.use <- which(row_stat >= min.cells)
  } else {
    features.use <- which(Matrix::rowSums(counts > 0) >= min.cells)
  }
  obj <- obj[features.use, ]
  obj
}

#' 合并后过滤基因
filter_gene_merge <- function(obj, min.cells = 0) {
  layers <- Layers(obj)
  keep_genes <- lapply(layers, function(layer) {
    counts <- LayerData(obj, layer = layer)
    if (inherits(x = counts, what = "IterableMatrix")) {
      rlang::check_installed(pkg = "BPCells", reason = "for working with BPCells")
      row_stat <- BPCells::matrix_stats(
        matrix = counts,
        row_stats = "nonzero"
      )$row_stats
      features.use <- which(row_stat >= min.cells)
    } else {
      features.use <- which(Matrix::rowSums(counts > 0) >= min.cells)
    }
    features.use
  }) |> unlist()
  keep_genes <- unique(names(keep_genes))
  obj <- obj[keep_genes, ]
}


#' clustree 绘制
#'
plot_clustree <- function(obj, filename = NULL) {
  # 分辨率聚类树
  # https://cloud.tencent.com/developer/article/1707335
  require(clustree)
  require(scales)
  # to do: 待解决：https://github.com/r-lib/scales/issues/413
  # number_si <- label_number(scale_cut = scales::cut_short_scale())
  res_tree <- clustree(obj) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set1") +
    scale_edge_color_continuous(low = "grey80", high = "#ee253a") +
    scale_size(range = c(4, 12)) +
    guides(edge_alpha = "none")
  if (!is.null(filename)) {
    ggsave(filename, res_tree, width = 10, height = 8)
  }
  res_tree
}
