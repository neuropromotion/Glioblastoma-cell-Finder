library(xgboost)
library(Seurat)
source("path_to/chromosome_means_function.R") # load get_chromosome_means function
predict_and_add_metadata <- function(seurat_object, model_weights_path = "path_to/xgboost.model") {
  counts <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
  chromosome_means <- get_chromosome_means(counts)
  barcodes <- rownames(chromosome_means)
  dval <- xgb.DMatrix(data = data.matrix(chromosome_means))
  model <- xgb.load(model_weights_path)
  
  pred_prob <- predict(model, dval)
  pred_cls  <- ifelse(pred_prob >= 0.5, 1L, 0L)
  
  ann <- ifelse(pred_cls == 1L, "Glioblastoma cells", "Stromal cells")
  names(ann) <- barcodes
  
  seurat_object$GBM_prob <- pred_prob  
  seurat_object$ML_annotation <- ann[Seurat::Cells(seurat_object)]
  seurat_object
}

# USAGE
#seurat.obj <- readRDS(file = "path_to_RDS/some_rds.rds") # your Seurat object
seurat.obj = predict_and_add_metadata(seurat.obj)
# Visualization:
DimPlot(seurat.obj, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(seurat.obj, features = "GBM_prob", cols = c('blue','black',"red")) # black cells - inconfident prediction
