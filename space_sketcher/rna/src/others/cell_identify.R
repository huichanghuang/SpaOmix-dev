#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-i", "--inputdir"), type="character",
                help="path of input(full path contain raw feature bc matrix)"),
    make_option(c("-o", "--outputdir"), type="character",
                help="path of output, default filtered under inputdir"),
    make_option(c("-u", "--minumi"), type = "integer", default = 200,
                help = "The minimum number of UMI"),
    make_option(c("-p", "--pvalue"), type="double", default=0.01,
                help="pvalue cutoff"),
    make_option(c("-f", "--force_cells"), type="integer", default=0,
                help="force cells")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$inputdir)){
  print_help(opt_parser)
  stop("No inputfile set", call.=FALSE)
}

if (!dir.exists(opt$inputdir)) {
    message("input dir does not exist")
    stop
}
inputdir = opt$inputdir

outputdir = ''
if(!is.null(opt$outputdir)){
    outputdir = opt$outputdir
}else{
    outputdir = paste0(dirname(opt$inputdir), '/callcell')
}
# 检查output目录是否已存在
if (dir.exists(outputdir)) {
    # 删除原有目录及其内容
    unlink(outputdir, recursive = TRUE)
    message(paste("Removed existing directory:", outputdir))
}

force_cells = opt$force_cells
minumi = opt$minumi
pvalue = opt$pvalue
print(paste0("inputdir:", opt$inputdir, "    outputdir:", outputdir,
             "    min_umi: ", minumi, "    pvalue:", pvalue, "    force_cells:", force_cells))


suppressMessages({
    library(DropletUtils)
    library(Matrix)
    library(dplyr)
    library(Seurat)
})

                             
output_10x_format_results <-function(outputdir,  matrix_read, matrix_counts, cell_barcode){
    # get cell only data
    gi <- rowData(matrix_read)[,1:1]
    gs <- rowData(matrix_read)[,2:2]
    bc <- colData(matrix_read)[,2][cell_barcode]
    umi.counts <- matrix_counts[,cell_barcode]     
    
    # output 10X format filtered_feature_bc_matrix
    write10xCounts(outputdir, umi.counts, gene.id=gi, gene.symbol=gs, barcodes=bc, version='3', overwrite=TRUE)
}


set.seed(0)
print("read and count raw matrix")
matrix_read <- read10xCounts(inputdir, col.names=TRUE, row.names = c("id", "symbol"))
matrix_counts <- counts(matrix_read)
barcode_counts <- colSums(matrix_counts)
if (force_cells == 0) {
    ###filter the cell barcodes by min umi
    out <- emptyDrops(counts(matrix_read), lower = minumi)
    nonan <- out[!is.na(out$FDR),]
    is.cell <- nonan[nonan$FDR < pvalue,]
    selected_barcodes <- row.names(is.cell)
    selected_indices <- which(names(barcode_counts) %in% selected_barcodes)
    print(paste("Number of filtered cell barcodes:", dim(is.cell)[1]))
    #output_10x_format_results
    output_10x_format_results(outputdir,  matrix_read, matrix_counts, selected_indices)

    ###extract the plot data
    br.out <- barcodeRanks(matrix_read, lower = minumi)
    plotdata <- data.frame(rank = br.out$rank, total = br.out$total, barcode = rownames(br.out))
    plotdata$is_cell_barcode <- ifelse(plotdata$barcode %in% row.names(is.cell), 1, 0)
    ###sort by rank
    plotdata <- plotdata[order(plotdata$rank),]
    plotdata <- plotdata[c("barcode", "total", "is_cell_barcode", "rank")]
    colnames(plotdata) <- c("barcode", "UMI", "is_cell_barcode", "rank")
    ###output the plot data
    outrank <- paste0(dirname(dirname(dirname(inputdir))), "/cell_rna_umi.rank.txt")
    write.table(plotdata, file = outrank, sep = "\t", quote = FALSE, row.names = FALSE)

} else {
    # force cell
    # Sort ALL barcodes by UMI (descending) and assign ranks
    all_barcodes <- names(barcode_counts)
    sorted_order <- order(barcode_counts, decreasing = TRUE)
    sorted_barcodes <- all_barcodes[sorted_order]

    selected_barcodes <- sorted_barcodes[1:min(force_cells, length(sorted_barcodes))]
    selected_indices <- which(names(barcode_counts) %in% selected_barcodes)
    #output_10x_format_results
    output_10x_format_results(outputdir,  matrix_read, matrix_counts, selected_indices)

    sorted_umi <- barcode_counts[sorted_order]
    ranks <- seq_along(sorted_barcodes)  # rank=1 is highest UMI
    # Generate plotdata (same structure as opt$force_cells==0)
    plotdata <- data.frame(
        barcode = sorted_barcodes,
        UMI = sorted_umi,
        rank = ranks,
        is_cell_barcode = ifelse(sorted_barcodes %in% selected_barcodes, 1, 0)
    )
    # Output plot data
    outrank <- paste0(dirname(dirname(dirname(inputdir))), "/cell_rna_umi.rank.txt")
    write.table(plotdata, file = outrank, sep = "\t", quote = FALSE, row.names = FALSE)
}

###Summary
summaryfile = paste0(dirname(inputdir), '/Summary.csv')
cellreadsstat = paste0(dirname(inputdir), '/CellReads.stats')

####cell reads, umi, genes summary
cellreadsdata <- read.table(cellreadsstat, header = TRUE, sep = "\t")
###remove row that CB = CBnotInPasslist
cellreadsdata <- cellreadsdata[cellreadsdata$CB != "CBnotInPasslist",]
truecellreads <- cellreadsdata[cellreadsdata$CB %in% selected_barcodes,]

estimated_number_of_cells <- dim(truecellreads)[1]
unique_reads_in_cells_mapped_to_gene <- sum(truecellreads$featureU)
fraction_of_unique_reads_in_cells <- round(unique_reads_in_cells_mapped_to_gene/sum(cellreadsdata$featureU), 4)
mean_reads_per_cell <- round(mean(truecellreads$featureU),0)
median_reads_per_cell <- round(median(truecellreads$featureU),0)
umi_in_cells <- sum(truecellreads$nUMIunique)
mean_umi_per_cell <- round(mean(truecellreads$nUMIunique),0)
median_umi_per_cell <- round(median(truecellreads$nUMIunique),0)
mean_gene_per_cell <- round(mean(truecellreads$nGenesUnique),0)
median_gene_per_cell <- round(median(truecellreads$nGenesUnique),0)

# 定义字段名与变量名的映射关系
field_var_mapping <- list(
  "Estimated Number of Cells" = "estimated_number_of_cells",
  "Unique Reads in Cells Mapped to GeneFull_Ex50pAS" = "unique_reads_in_cells_mapped_to_gene",
  "Fraction of Unique Reads in Cells" = "fraction_of_unique_reads_in_cells",
  "Mean Reads per Cell" = "mean_reads_per_cell",
  "Median Reads per Cell" = "median_reads_per_cell",
  "UMIs in Cells" = "umi_in_cells",
  "Mean UMI per Cell" = "mean_umi_per_cell",
  "Median UMI per Cell" = "median_umi_per_cell",
  "Mean GeneFull_Ex50pAS per Cell" = "median_gene_per_cell",
  "Median GeneFull_Ex50pAS per Cell" = "median_gene_per_cell"
)
# 批量替换函数
batch_replace <- function(data, mapping) {
  sapply(data, function(line) {
    for (field in names(mapping)) {
      if (grepl(paste0("^", field, ","), line)) {
        var_name <- mapping[[field]]
        if (exists(var_name, envir = .GlobalEnv)) {
          new_value <- get(var_name, envir = .GlobalEnv)
          return(paste0(field, ",", new_value))
        }
      }
    }
    return(line)
  }, USE.NAMES = FALSE)
}


summary_data <- readLines(summaryfile)
updated_data <- batch_replace(summary_data, field_var_mapping)

outsummary <- paste0(dirname(inputdir), '/Summary.callcell.csv')
writeLines(updated_data, outsummary)
