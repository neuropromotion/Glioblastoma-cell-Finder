library(data.table)
library(Seurat)
library(DoubletFinder)
library(celda)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(ggplot2)
cols <- c(
  "Stromal" = 'gray55',                
  "GBM" = 'gold1' 
)
feature_plot <- function(data, feature){
  FeaturePlot(data, feature, order=T, slot = 'data')  +
    scale_color_gradientn(
      colours = c("gray", "blue1", "deeppink1"),  
      values = c(0, 0.3, .65,  1),  
    )
}
count_plots <- function(obj, 
                        feature_lower=600,
                        feature_upper=6000,
                        count_lower=200,
                        count_upper=1e5,
                        mito=15,
                        pt.size=0){
  
  vln_plot1 <- VlnPlot(obj, features = "nFeature_RNA", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = feature_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = feature_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot2 <- VlnPlot(obj, features = "nCount_RNA", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = count_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = count_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot3 <- VlnPlot(obj, features = "mitoPercent", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = mito, linetype = "dashed", color = "red", cex=1) +
    NoLegend()
  combined_plot <- vln_plot1 | vln_plot2 | vln_plot3
  return(combined_plot)
}
source("/home/amismailov/R/chromosome_means_function.R") # chromosome_means function

deduplex <- function(obj){
  nExp <- round(0.05 * ncol(obj))
  pK <- 0.09 
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = 0.09, 
                       nExp = nExp, reuse.pANN = NULL, sct = FALSE)
  
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  cat('Cells Before deduplex: ', ncol(obj), '\n')
  obj <- subset(
    obj,
    subset = !!rlang::sym(df_col) == "Singlet"
  )
  cat('Cells Before deduplex: ', ncol(obj), '\n')
  return(obj)
} # DoubletFinder
stantard_workflow <- function(seu, n_components=15, res=0.1){
  seu <- NormalizeData(seu)                             
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)   
  seu <- ScaleData(seu, features = rownames(seu)) 
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:n_components)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:n_components)
  seu <- deduplex(seu)
  return(seu)
} 
filter_ol <- function(obj){
  cells_both <- WhichCells(
    obj,
    expression = MAG > 1 & PTPRZ1 > 1 & MOG > 1,    
    slot = "counts"                    
  )
  obj <- subset(obj, cells = cells_both, invert = TRUE)
  cat(length(cells_both), 'cells were filtered [OL filter]')
  return(obj)
} # filter doubtful ols
filter_macro <- function(obj){
  cells_both <- WhichCells(
    obj,
    expression = PTPRZ1 > 1 & PTPRC > 1,    
    slot = "counts"                    
  )
  obj <- subset(obj, cells = cells_both, invert = TRUE)
  cat(length(cells_both), 'cells were filtered [Immune filter]')
  return(obj)
} # filter doubtful macrophages

# TRAIN/VAL DATASETS 
# GBM (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89567
# GMB (10x GENOMICS)
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-whole-transcriptome-analysis-1-standard-4-0-0
# GMB (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135045


#--------------------------------------smatseq2 (GSE89567)------------------------
create_sue <-function(path="path_to/GSE89567.txt.gz"){
  dt <- fread(path)
  dt[[1]] <- trimws(gsub("[\"']", "", dt[[1]])) 
  #dt[[1]]
  num_dt <- dt[, -1, with = FALSE]
  for (j in seq_along(num_dt)) {
    set(num_dt, j = j, value = as.numeric(trimws(num_dt[[j]])))
  }
  mat <- as.matrix(num_dt)
  rownames(mat) <- dt[[1]]
  mat[1:5,1:5]
  seu <- CreateSeuratObject(counts = mat, assay = "RNA")
  return(seu)
}
seu <- create_sue()

seu$mitoPercent <- PercentageFeatureSet(seu, pattern = '^MT-')
count_plots(seu, feature_upper = 9e3, count_upper = 20e3)
seu <- subset(
  seu,
  subset = nFeature_RNA <= 9e3 &
    nCount_RNA <= 20e3 
)
seu <- stantard_workflow(seu, n_components = 5, res=0.01) 

DimPlot(seu, label=T)
FeaturePlot(seu, 'MAG')

# DEG assay: (optional)
# FAM <- FindAllMarkers(seu,
#                       logfc.threshold = 0.25,
#                       min.pct = 0.1,
#                       only.pos = TRUE,
#                       test.use = 'wilcox',
#                       layer = 'data')
# 
# top.genes = list()
# for (cl in 0:9) {
#   cluster_genes <- FAM$gene[FAM$cluster == cl][1:30]
#   curr_name <- as.character(cl)
#   top.genes[[curr_name]] = cluster_genes 
# }

# rename
new_cluster_names <- c(
  '0' = 'GBM',
  '1' = 'Stromal',
  '2' = 'Stromal')
seu <- RenameIdents(seu, new_cluster_names)

DimPlot(seu, label=F, cols = cols)

# chromosome means extraction
seu.cancer <- WhichCells(seu, idents = 'GBM')
seu.cancer <- subset(seu, cells = seu.cancer)

seu.stromal <- WhichCells(seu, idents = 'Stromal')
seu.stromal <- subset(seu, cells = seu.stromal)

cancer.barcodes <- colnames(seu.cancer)
stromal.barcodes <- colnames(seu.stromal)

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
chromosome_means = get_chromosome_means(counts)


keep <- cancer.barcodes[cancer.barcodes %in% rownames(chromosome_means)]
chromosome_means.cancer <- chromosome_means[keep, , drop = FALSE]

keep <- stromal.barcodes[stromal.barcodes %in% rownames(chromosome_means)]
chromosome_means.stromal <- chromosome_means[keep, , drop = FALSE]


write.csv(
  data.frame(chromosome_means.cancer),
  file = "chromosome_means_smartseq2.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)
write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_smartseq2_stromal.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)




#--------------------------------------10X------------------------------------------
data.dir <- "path_to_dir/filtered_feature_bc_matrix"
m <- Read10X(data.dir = data.dir)
seu <- CreateSeuratObject(counts = m, project = "GBM", min.cells = 3, min.features = 200)

seu$mitoPercent <- PercentageFeatureSet(seu, pattern = '^MT-')
seu <- subset(
  seu,
  subset = nFeature_RNA <= 7000 & nFeature_RNA > 400 &
    nCount_RNA <= 30e3 & nCount_RNA >= 1e3 & mitoPercent <= 15
)
seu <- stantard_workflow(seu, n_components = 5, res=0.1) 

DimPlot(seu)
feature_plot(seu, 'EGFR')

# optional:
# FAM <- FindAllMarkers(seu,
#                       logfc.threshold = 0.25,
#                       min.pct = 0.1,
#                       only.pos = TRUE,
#                       test.use = 'wilcox',
#                       slot = 'data')
# 

# top.genes = list()
# set.of.top.genes = vector()
# 
# for (cl in 0:5) {
#   cluster_genes <- FAM$gene[FAM$cluster == cl][1:30]   # гены текущего кластера
#   curr_name <- as.character(cl)#paste0('cluster', as.character(cl))
#   top.genes[[curr_name]] = cluster_genes
#   set.of.top.genes = c(set.of.top.genes, cluster_genes)
# }

DimPlot(seu)
feature_plot(seu, 'PTPRZ1')
new_cluster_names <- c(
  '0' = 'GBM',
  '1' = 'Stromal',
  '2' = 'Stromal',
  '3' = 'Unconfindent'
  )  
seu <- RenameIdents(seu, new_cluster_names)
seu <- subset(seu, idents = "Unconfindent", invert = TRUE) # remove unconfident cells
DimPlot(seu, label=F, cols = c('gold1', 'gray55'))

#CNV
seu.cancer <- WhichCells(seu, idents = 'GBM')
seu.cancer <- subset(seu, cells = seu.cancer)

seu.stromal <- WhichCells(seu, idents = 'Stromal')
seu.stromal <- subset(seu, cells = seu.stromal)

cancer.barcodes <- colnames(seu.cancer)
stromal.barcodes <- colnames(seu.stromal)

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
chromosome_means = get_chromosome_means(counts)


keep <- cancer.barcodes[cancer.barcodes %in% rownames(chromosome_means)]
chromosome_means.cancer <- chromosome_means[keep, , drop = FALSE]

keep <- stromal.barcodes[stromal.barcodes %in% rownames(chromosome_means)]
chromosome_means.stromal <- chromosome_means[keep, , drop = FALSE]



write.csv(
  data.frame(chromosome_means.cancer),
  file = "chromosome_means_10x_1.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)

write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_10x_1_stromal.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)

#-----------------------------------smartseq2 GSE135045-----------------------------
# data processed in main_branch.R file 
merged_filtered <- readRDS(file = "path_to/GSE135045_processed.rds")
Idents(merged_filtered) <- 'init_annot'
DimPlot(merged_filtered)

gbm <- WhichCells(merged_filtered, idents = 'GBM')
gbm <- subset(merged_filtered, cells = gbm)

stromal <- WhichCells(merged_filtered, idents = c('Macrophages/Microglia', 'Oligodendrocytes', 'T lymphocytes',
                                                  'CAFs', 'Endothelial cells'))
stromal <- subset(merged_filtered, cells = stromal)

#gbm
counts.gbm <- GetAssayData(gbm, assay = "RNA", layer = "counts")
chromosome_means.gbm = get_chromosome_means(counts.gbm)
#stromal
counts.stromal <- GetAssayData(stromal, assay = "RNA", layer = "counts")
chromosome_means.stromal = get_chromosome_means(counts.stromal)
#save
write.csv(
  data.frame(chromosome_means.gbm),
  file = "chromosome_means_MF_gbm.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)

write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_MF_stromal.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)

