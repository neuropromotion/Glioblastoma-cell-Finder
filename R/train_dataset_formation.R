library(data.table)
library(Seurat)
library(DoubletFinder)
library(celda)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
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
source("/home/amismailov/R/chromosome_means_function.R") #get_chromosome_means func
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
} # filter doubtful OLs
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

# GBM (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89567
# GMB (10x GENOMICS)
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-whole-transcriptome-analysis-1-standard-4-0-0
# GMB (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135045



#--------------------------------------smatseq2 GSE89567------------------------
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

new_cluster_names <- c(
  '0' = 'GBM',
  '1' = 'Stromal',
  '2' = 'Stromal')
seu <- RenameIdents(seu, new_cluster_names)

jpeg("vln.jpg", width = 1500, height = 1500, res=300)
DimPlot(seu, label=F, cols = c('gold1', 'gray55'))
dev.off()


# divide data into stromal and GBM part
seu.cancer <- WhichCells(seu, idents = 'GBM')
seu.cancer <- subset(seu, cells = seu.cancer)

seu.stromal <- WhichCells(seu, idents = 'Stromal')
seu.stromal <- subset(seu, cells = seu.stromal)

cancer.barcodes <- colnames(seu.cancer)
stromal.barcodes <- colnames(seu.stromal)

# get chromosome means 
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
chromosome_means = get_chromosome_means(counts)

# keep GBM cells
keep <- cancer.barcodes[cancer.barcodes %in% rownames(chromosome_means)]
chromosome_means.cancer <- chromosome_means[keep, , drop = FALSE]
# keep stromal cells
keep <- stromal.barcodes[stromal.barcodes %in% rownames(chromosome_means)]
chromosome_means.stromal <- chromosome_means[keep, , drop = FALSE]

# save
write.csv(
  data.frame(chromosome_means.cancer),
  file = "chromosome_means_smartseq2.csv",
  row.names = TRUE,            
  quote = FALSE                 
)
write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_smartseq2_stromal.csv",
  row.names = TRUE,            
  quote = FALSE                
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
FeaturePlot(seu, 'PTPRZ1')
new_cluster_names <- c(
  '0' = 'GBM',
  '1' = 'Stromal',
  '2' = 'Stromal',
  '3' = 'Unconfindent'
  )  
seu <- RenameIdents(seu, new_cluster_names)
seu <- subset(seu, idents = "Unconfindent", invert = TRUE) # delete unconfident cells
DimPlot(seu, label=F, cols = c('gold1', 'gray55'))

# divide data into stromal and GBM part
seu.cancer <- WhichCells(seu, idents = 'GBM')
seu.cancer <- subset(seu, cells = seu.cancer)
 
seu.stromal <- WhichCells(seu, idents = 'Stromal')
seu.stromal <- subset(seu, cells = seu.stromal)

cancer.barcodes <- colnames(seu.cancer)
stromal.barcodes <- colnames(seu.stromal)

# get chromosome means 
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
chromosome_means = get_chromosome_means(counts)


# keep GBM cells
keep <- cancer.barcodes[cancer.barcodes %in% rownames(chromosome_means)]
chromosome_means.cancer <- chromosome_means[keep, , drop = FALSE]

# keep stromal cells
keep <- stromal.barcodes[stromal.barcodes %in% rownames(chromosome_means)]
chromosome_means.stromal <- chromosome_means[keep, , drop = FALSE]



# save
write.csv(
  data.frame(chromosome_means.cancer),
  file = "chromosome_means_10x_1.csv",
  row.names = TRUE,           
  quote = FALSE                 
)

write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_10x_1_stromal.csv",
  row.names = TRUE,            
  quote = FALSE                
)

#-----------------------------------smartseq2 GSE135045-----------------------------
# data prepared in https://github.com/neuropromotion/CAFs-in-glioblastoma-microenvironment/blob/main/main_branch.R
merged_filtered <- readRDS(file = "path_to_RDS/merged_filtered_v2.rds")

# choose GBM and stromal cells by marker co-expression (rest cells will be advanced validation)
cancer.cells <- WhichCells(merged_filtered, expression = PTPRZ1 > 1 & SOX2 > 1 & EGFR > 1 & CD163 == 0 & CD3E == 0 & RGS5 == 0)
length(cancer.cells)
macro.cells <- WhichCells(merged_filtered, expression = PTPRC > 1 & C3 > 1 & C1QA > 1 & PTPRZ1 == 0 & CD163 > 1)
length(macro.cells)
ol.cells <- WhichCells(merged_filtered, expression = MOG > 1 & MAG > 1 & SOX2 == 0 & PTPRZ1 == 0)
length(ol.cells)

write.csv(data.frame(cell = cancer.cells), "cancer_barcodes.csv", row.names = FALSE) 
write.csv(data.frame(cell = macro.cells), "macro_barcodes.csv", row.names = FALSE) 
write.csv(data.frame(cell = ol.cells), "OL_barcodes.csv", row.names = FALSE) 

counts <- GetAssayData(merged_filtered, assay = "RNA", layer = "counts")
chromosome_means = get_chromosome_means(counts)
#save
write.csv(
  data.frame(chromosome_means),
  file = "chromosome_means_MF_full.csv",
  row.names = TRUE,            
  quote = FALSE                 
)
