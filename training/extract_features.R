library(data.table)
library(Seurat)
library(DoubletFinder)
library(celda)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
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
get_chromosome_means <- function(counts, path_to_mapped_genes='/home/amismailov/mart_export.txt'){
  row_sums <- rowSums(counts)
  counts <- counts[row_sums > 0, ]
  norm_counts <- t(t(counts) / colSums(counts) * 1000)  
  log_counts <- log2(norm_counts + 1) 
  gene_filter <- rowSums(log_counts) >= 100
  log_counts_filtered <- log_counts[gene_filter, ]
  # Remove HLA-genes
  hla_genes <- grep("HLA", rownames(log_counts_filtered), value = TRUE)
  log_counts_filtered <- log_counts_filtered[!(rownames(log_counts_filtered) %in% hla_genes), ]
  gene_chromosome_map <- read.table(path_to_mapped_genes, sep = '\t', header = T)  # Файл с генами и их хромосомами
  chromosome_means <- sapply(unique(gene_chromosome_map$Chromosome.scaffold.name), function(chrom) {
    genes <- intersect(
      gene_chromosome_map$HGNC.symbol[gene_chromosome_map$Chromosome.scaffold.name == chrom],
      rownames(log_counts_filtered)  
    )
    
    if (length(genes) > 0) {
      colMeans(log_counts_filtered[genes, , drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, ncol(log_counts_filtered))  # Если генов нет, возвращаем NA
    }
  })
  chromosome_means <- as.data.frame(chromosome_means)
  colnames(chromosome_means) <- paste0("Chr", rep(1:22))
  rownames(chromosome_means) <- rownames(t(log_counts_filtered))
  return(chromosome_means)
}
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


#decontex_workflow <- function(seu){
  counts_mat <- GetAssayData(seu, layer = "counts", assay = "RNA") # raw counts matrix extraction
  cell_md <- seu@meta.data # metadata extraction
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = cell_md
  ) 
  set.seed(42) # fix seed to reproducing 
  # Add clusters
  get_groups <- function(sobj){
    sobj <- NormalizeData(sobj, verbose = FALSE)
    sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
    sobj <- ScaleData(sobj, verbose = FALSE)
    sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
    sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
    sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
    return(sobj@meta.data[['seurat_clusters']])
  }
  add_groups <- function(sobj){
    sobj$soup_group <- get_groups(sobj)
    return(sobj)
  }
  seu <- add_groups(seu)
  sce <- decontX(sce, z = seu$soup_group) # run decontX
  # Extraction clear count and metadata
  clean_counts <- decontXcounts(sce)
  clean_md     <- as.data.frame(colData(sce))
  # create new clean Seurat-object
  seu_clean <- CreateSeuratObject(
    counts    = clean_counts,
    project   = "GBM_decontaminated",
    meta.data = clean_md
  )
  return(seu_clean)
}
#stantard_workflow <- function(seu, n_components=15, res=0.1){
  seu <- JoinLayers(seu)
  seu <- decontex_workflow(seu)
  seu <- NormalizeData(seu)                             
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)   
  seu <- ScaleData(seu, features = rownames(seu)) 
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  # seu <- RunHarmony(
  #   object = seu,
  #   group.by.vars = "Sample"
  # )
  seu <- FindNeighbors(seu, dims = 1:n_components)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, dims = 1:n_components)
  seu <- deduplex(seu)
  return(seu)
}
# GBM (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89567
# GMB (10x GENOMICS)
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-whole-transcriptome-analysis-1-standard-4-0-0
# GMB (Smart-seq2)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135045



#--------------------------------------GSE89567 (smatseq2)------------------------
dt <- fread("/mnt/jack-5/amismailov/CAF_study/GCF_2/GSE89567.txt.gz")

genes <- gsub("[\"']", "", dt[[1]])
genes <- trimws(genes)

mat <- as.matrix(dt[, -genes, with = FALSE])
rownames(mat) <- genes
storage.mode(mat) <- "double"
seu <- CreateSeuratObject(counts = mat, assay = "RNA")

seu$mitoPercent <- PercentageFeatureSet(seu, pattern = '^MT-')

seu <- subset(
  seu,
  subset = nFeature_RNA <= 8000 &
    nCount_RNA <= 20e3 & nCount_RNA >= 10e3 & mitoPercent <= 15
)
seu <- stantard_workflow(seu, n_components = 5, res=0.01) 

DimPlot(seu, label=T)
FeaturePlot(seu, 'PTPRC')

# FAM <- FindAllMarkers(seu,
#                       logfc.threshold = 0.25,
#                       min.pct = 0.1,
#                       only.pos = TRUE,
#                       test.use = 'wilcox',
#                       slot = 'data')
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
data.dir <- "/mnt/jack-5/amismailov/CAF_study/GCF_2/filtered_feature_bc_matrix"


m <- Read10X(data.dir = data.dir)
seu <- CreateSeuratObject(counts = m, project = "GBM", min.cells = 3, min.features = 200)

seu$mitoPercent <- PercentageFeatureSet(seu, pattern = '^MT-')
seu <- subset(
  seu,
  subset = nFeature_RNA <= 7000 & nFeature_RNA > 400 &
    nCount_RNA <= 30e3 & nCount_RNA >= 1e3 & mitoPercent <= 15
)
seu <- stantard_workflow(seu, n_components = 5, res=0.01) 

feature_plot(seu, 'PTPRZ1')


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
  '2' = 'Stromal'
  )  
seu <- RenameIdents(seu, new_cluster_names)

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



#-----------------------------------smartseq2 (validation dataset)-----------------------------
path = '/mnt/jack-5/amismailov/CAF_study'
dirs <- c("sample1", "sample2","sample3","sample4","sample5","sample6","sample7")
for (x in dirs) {
  file_path <- paste0(path, '/', x, '/data.txt')
  expression_matrix <- read.delim(file_path, row.names = 1)
  assign(x, CreateSeuratObject(counts = expression_matrix))
}
seu <- merge(sample1, y = c(sample2, sample3, sample4, sample5, sample6, sample7),
                add.cell.id  = dirs,
                project = 'GBM')

seu$sample <- rownames(seu@meta.data) 
seu@meta.data <- separate(seu@meta.data, col = 'sample', into = c('Sample', 'Barcode'),
                             sep = '_')
 
seu$mitoPercent <- PercentageFeatureSet(seu, pattern = '^MT-')
# filter cells
seu <- subset(seu, subset = nFeature_RNA > 600 & nCount_RNA <= 10e3 &
                            nFeature_RNA < 4500 & mitoPercent < 20)

seu <- JoinLayers(object = seu)
seu <- stantard_workflow(seu, n_components = 5, res=0.01) 

DimPlot(seu)
feature_plot(seu, 'CD163')

#TRAIN SELECTION
cancer.cells <- WhichCells(seu, expression = PTPRZ1 > 1 & SOX2 > 1 & EGFR > 1 & CD163 == 0 & CD3E == 0 & RGS5 == 0)
length(cancer.cells)
macro.cells <- WhichCells(seu, expression = PTPRC > 1 & C3 > 1 & C1QA > 1 & PTPRZ1 == 0 & CD163 > 1)
length(macro.cells)
ol.cells <- WhichCells(seu, expression = MOG > 1 & MAG > 1 & SOX2 == 0 & PTPRZ1 == 0)
length(ol.cells)

write.csv(data.frame(cell = cancer.cells), "cancer_barcodes.csv", row.names = FALSE) 
write.csv(data.frame(cell = macro.cells), "macro_barcodes.csv", row.names = FALSE) 
write.csv(data.frame(cell = ol.cells), "OL_barcodes.csv", row.names = FALSE) 


anno <- setNames(rep("Validation", ncol(seu)), colnames(seu))            # всё остальное — Cluster4 [web:154]
anno[cancer.cells] <- "GBM cells train"                                            # метки для первой группы [web:154]
anno[macro.cells] <- "Stromal cells train"                                            # метки для второй группы [web:154]
anno[ol.cells] <- "Stromal cells train"  

seu$split <- unname(anno[colnames(seu)]) 
DimPlot(seu, group.by = 'split', cols = c('gold1', 'springgreen2', 'gray'))
Idents(seu) <- 'split'

#CNV
seu.cancer <- WhichCells(seu, idents = 'GBM cells train')
seu.cancer <- subset(seu, cells = seu.cancer)

seu.stromal <- WhichCells(seu, idents = 'Stromal cells train')
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
  file = "chromosome_means_smartseq_val_cancer.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)

write.csv(
  data.frame(chromosome_means.stromal),
  file = "chromosome_means_smartseq_val_stromal.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)



#---------FOR VALIDATION (metrics measurement)
# same dataset
gbm.cells <- WhichCells(merged_filtered, idents = 'GBM') 
gbm.cells_filtered <- setdiff(gbm.cells, cancer.cells) # filter train barcodes
length(gbm.cells_filtered)

stromal.cells <- WhichCells(merged_filtered, idents = c('T lymphocytes', 
                                                        'Macrophages/Microglia', 
                                                        'Oligodendrocytes', 
                                                        'CAFs', 
                                                        'Endothelial cells'))
stromal.cells_filtered <- setdiff(stromal.cells, macro.cells)
stromal.cells_filtered <- setdiff(stromal.cells_filtered, ol.cells) # filter train barcodes

length(stromal.cells_filtered)

write.csv(data.frame(cell = gbm.cells_filtered), "cancer_barcodes_validation.csv", row.names = FALSE) 
write.csv(data.frame(cell = stromal.cells_filtered), "stromal_barcodes_validation.csv", row.names = FALSE) 



