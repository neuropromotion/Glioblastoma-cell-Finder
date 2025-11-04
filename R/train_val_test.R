library(yardstick)
library(xgboost)
library(Matrix)
library(tidymodels)
library(vip)
library(doParallel)

source("/home/amismailov/R/chromosome_means_function.R") # load chromosome_means function

#-----------------------TRAIN-------------------------- 
# data processed in train_dataset_formation.R file 
cancers <- c("chromosome_means_smartseq2.csv", 
             "chromosome_means_10x_1.csv", 
             "chromosome_means_MF_gbm.csv")

stromals <- c("chromosome_means_smartseq2_stromal.csv", 
              "chromosome_means_10x_1_stromal.csv",
              "chromosome_means_MF_stromal.csv")


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

dim(df.cancer) #17747    22
dim(df.stromal) #13468    22

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
  eval_metric = "auc",
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
#   [497]	train-auc:0.999996	valid-auc:0.998026

#-------------SAVE Model weigths----------------
xgb.save(bst, "boosting/xgboost.model")

#---------------TEST_1--------------------
# Dataset URL:
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0
# data processed in validation_1.R file (CAFs in glioblastoma repository)
val <- readRDS(file = "path_to/validation_1.rds")
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


#DimPlot(val)
DimPlot(val, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(val, features = "GBM_prob",
            cols = c('blue','black',"red"))

#---------------TEST_2--------------------
# Dataset URL:
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0
# data processed in validation_2.R file (CAFs in glioblastoma repository)
val <- readRDS(file = "path_to/validation_2.rds")
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


#DimPlot(val, label=F)
DimPlot(val, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
FeaturePlot(val, features = "GBM_prob",
            cols = c('blue','black',"red"))
