# Glioblastoma Cell Finder  
### Gradient Boosting Model for Cell-Type Annotation in R

This repository provides a **gradient boosting–based prediction model** for identifying glioblastoma-associated cell populations in **single-cell RNA-seq** datasets using **R**.  
The model is designed for convenient use with **Seurat** workflows.

---

## 📦 Files

| File | Description |
|------|-------------|
| `model.R` | Predtreinde model, requires weights (`xgboots.model`). Contains the main function `predict_and_add_metadata` |
| `xgboots.model` | model weights |
| `chromosome_means_function.R` | calculates chromosome means per cell |
| `train_dataset_formation.R` | TRAIN data collection pipeline |
| `train_val_test.R` | code for train/validation and independent test running |
| `gene_mapping.txt` | contains gene mapping on chromosome data (from GTF) | 


---

## Usage

The core function provided by this repository is:

```r
predict_and_add_metadata(seurat_object, weights_path='xboost.model')


