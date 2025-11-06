library(yardstick)
library(xgboost)
library(Matrix)
library(tidymodels)
library(vip)
library(doParallel)

source("chromosome_means_function.R")# get_chromosome_means func

#-----------------------TRAIN-------------------------- 
# data processed in train_data_formation.R file 

cancers <- c("chromosome_means_smartseq2.csv", 
             "chromosome_means_10x_1.csv")

stromals <- c("chromosome_means_smartseq2_stromal.csv", 
              "chromosome_means_10x_1_stromal.csv")
df.cancer <- data.frame()
df.stromal <- data.frame() 

for (file in cancers){
  current <- read.csv(file, stringsAsFactors = FALSE)
  current$X <- NULL
  df.cancer <- rbind(df.cancer, current)
}
for (file in stromals){
  current <- read.csv(file, stringsAsFactors = FALSE)
  current$X <- NULL
  df.stromal <- rbind(df.stromal, current)
}

# merged_filtered dataset:
df.merged_filtered <- as.data.frame(read.csv("chromosome_means_MF_full.csv", stringsAsFactors = FALSE)) # Full 
gbm.barcodes <- read.csv("cancer_barcodes.csv", stringsAsFactors = FALSE) # GBM barcodes
macro.barcodes <- read.csv("macro_barcodes.csv", stringsAsFactors = FALSE) # macrophages barcodes
ol.barcodes <- read.csv("OL_barcodes.csv", stringsAsFactors = FALSE) # OL barcodes
# filter out 
mf.gbm <- df.merged_filtered[df.merged_filtered$X %in% gbm.barcodes$cell, , drop = FALSE]
mf.macro <- df.merged_filtered[df.merged_filtered$X %in% macro.barcodes$cell, , drop = FALSE]
mf.ol <- df.merged_filtered[df.merged_filtered$X %in% ol.barcodes$cell, , drop = FALSE]
# bind macrophages and OLs in stromal:
mf.stromal = rbind(mf.macro, mf.ol)
mf.gbm$X=NULL # delete barcodes
mf.stromal$X=NULL # delete barcodes

# merge 
df.cancer <- rbind(df.cancer, mf.gbm)
df.stromal <- rbind(df.stromal, mf.stromal)

dim(df.cancer) #8110x22 - unbalanced
dim(df.stromal) #4679x22
# 12789 cells for train eventually 

# target
df.cancer$target  <- factor("cancer",  levels = c("stromal","cancer"))
df.stromal$target <- factor("stromal", levels = c("stromal","cancer"))

# merge
df <- rbind(df.cancer, df.stromal)
# shuffle
df <- df[sample.int(nrow(df)), , drop = FALSE]
# split 85/15
idx_stromal <- which(df$target == "stromal")
idx_cancer  <- which(df$target == "cancer")
train_frac <- 0.80
n_stromal_train <- floor(length(idx_stromal) * train_frac)
n_cancer_train  <- floor(length(idx_cancer)  * train_frac)
set.seed(1234)
tr_stromal <- sample(idx_stromal, n_stromal_train)
tr_cancer  <- sample(idx_cancer,  n_cancer_train)

train_idx <- c(tr_stromal, tr_cancer)
valid_idx <- setdiff(seq_len(nrow(df)), train_idx)

train <- df[train_idx, , drop = FALSE]
valid <- df[valid_idx, , drop = FALSE]


feature_cols <- setdiff(names(df), "target")

X_train <- as.matrix(train[, feature_cols, drop = FALSE])
y_train <- as.integer(train$target) - 1L   # 0=stromal, 1=cancer

X_valid <- as.matrix(valid[, feature_cols, drop = FALSE])
y_valid <- as.integer(valid$target) - 1L

dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dvalid <- xgb.DMatrix(data = X_valid, label = y_valid)

# params
params <- list(
  objective = "binary:logistic",
  eval_metric = "aucpr",
  tree_method = "hist"
)
params$eta <- 0.2              # более стабильное обучение
params$max_depth <- 4          # 22 признака: глубины 3–5 обычно достаточно
params$min_child_weight <- 2   # предотвращает переобучение на малых узлах
params$gamma <- 1             # можно поднять до 0.5–1, если переобучается

params$subsample <- 0.8
params$colsample_bytree <- 0.8

params$reg_lambda <- 2         # L2
params$reg_alpha  <- 0.1       # L1 для разреживания
 


bst <- xgb.train(
  params = params, 
  data = dtrain,
  nrounds = 2000,
  watchlist = list(train = dtrain, valid = dvalid),
  early_stopping_rounds = 25,
  maximize = TRUE,
  scale_pos_weight=dim(df.cancer)[1]/dim(df.stromal)[1]
)

# Stopping. Best iteration:
# train-aucpr:0.999987	valid-aucpr:0.999613

#-------------SAVE Model weigths----------------
xgb.save(bst, "boosting/xgboost.model")

#---------------VALIDATION--------------------
# data prepared in https://github.com/neuropromotion/CAFs-in-glioblastoma-microenvironment/blob/main/main_branch.R
merged_filtered <- readRDS(file = "path_to_RDS/merged_filtered_v2.rds")
counts <- GetAssayData(merged_filtered, assay = "RNA", layer = "counts")
chromosome_means <- get_chromosome_means(counts)
barcodes <- rownames(chromosome_means)
#DMatrix
dval <- xgb.DMatrix(data = data.matrix(chromosome_means))
model <- xgb.load("boosting/xgboost.model")
pred_prob <- predict(model, dval)            
pred_cls  <- ifelse(pred_prob >= 0.5, 1L, 0L)

df_predictions <- data.frame(
  barcode = barcodes,
  predicted_class = pred_cls,
  predicted_prob = pred_prob
)
merged_filtered$GBM_prob <- pred_prob  
ann <- ifelse(df_predictions$predicted_class == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- df_predictions$barcode
merged_filtered$ML_annotation <- ann[Cells(merged_filtered)]


#DimPlot(val)
#DimPlot(merged_filtered, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(merged_filtered, features = "GBM_prob",
            cols = c('blue','black',"red"))


#---------------TEST_1--------------------
# Dataset URL:
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0
# data processed in https://github.com/neuropromotion/CAFs-in-glioblastoma-microenvironment/blob/main/Validation_1.R
val <- readRDS(file = "path_to_RDS/test_1.rds")
counts <- GetAssayData(val, assay = "RNA", layer = "counts")
chromosome_means <- get_chromosome_means(counts)

barcodes <- rownames(chromosome_means)
# DMatrix
dval <- xgb.DMatrix(data = data.matrix(chromosome_means))
model <- xgb.load("boosting/xgboost.model")
pred_prob <- predict(model, dval)            
pred_cls  <- ifelse(pred_prob >= 0.5, 1L, 0L)

df_predictions <- data.frame(
  barcode = barcodes,
  predicted_class = pred_cls,
  predicted_prob = pred_prob
)
val$GBM_prob <- pred_prob  
ann <- ifelse(df_predictions$predicted_class == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- df_predictions$barcode
val$ML_annotation <- ann[Cells(val)]


DimPlot(val, group.by = 'init_annot')
DimPlot(val, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(val, features = "GBM_prob") + 
  scale_color_gradientn(
    colours = c("blue", 'black', "red"),  # Четыре цвета
    values = c(0, 0.4, .65,  1),  # Точки привязки (нормализованные от 0 до 1)
    limits = c(0, 1)  # Соответствует col.min и col.max в DotPlot
  )

#---------------TEST_2--------------------
# Dataset URL:
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0
# data processed in https://github.com/neuropromotion/CAFs-in-glioblastoma-microenvironment/blob/main/Validation_2.R
val <- readRDS(file = "path_to_RDS/test_2.rds")
counts <- GetAssayData(val, assay = "RNA", layer = "counts")
chromosome_means <- get_chromosome_means(counts)

barcodes <- rownames(chromosome_means)
# DMatrix
dval <- xgb.DMatrix(data = data.matrix(chromosome_means))
model <- xgb.load("boosting/xgboost.model")
pred_prob <- predict(model, dval)            
pred_cls  <- ifelse(pred_prob >= 0.5, 1L, 0L)

df_predictions <- data.frame(
  barcode = barcodes,
  predicted_class = pred_cls,
  predicted_prob = pred_prob
)
val$GBM_prob <- pred_prob   
ann <- ifelse(df_predictions$predicted_class == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- df_predictions$barcode
val$ML_annotation <- ann[Cells(val)]


DimPlot(val, label=F, cols=c('red', 'blue', 'purple', 'skyblue', 'gold1'))
DimPlot(val, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(val, features = "GBM_prob") + 
  scale_color_gradientn(
    colours = c("blue", 'black', "red"),  # Четыре цвета
    values = c(0, 0.4, .65,  1),  # Точки привязки (нормализованные от 0 до 1)
    limits = c(0, 1)  # Соответствует col.min и col.max в DotPlot
  )
