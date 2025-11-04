# Glioblastoma Cell Finder  
### Gradient Boosting Model for Cell-Type Annotation in R

This repository provides a **gradient boosting–based prediction model** for identifying glioblastoma-associated cell populations in **single-cell RNA-seq** datasets using **R**.  
The model is designed for convenient use with **Seurat** workflows.

---

## 📦 Files

| File | Description |
|------|-------------|
| `model.R` | Contains the main function `predict_and_add_metadata` |
| `weights/` *(optional)* | Directory for storing trained model weight files |

---

## 🔧 Requirements

- R (≥ 4.0)
- Seurat (≥ 4.0)
- dplyr
- data.table  
*(Additional dependencies may be required depending on model training environment)*

---

## 🚀 Usage

The core function provided by this repository is:

```r
predict_and_add_metadata(seurat_object, weights_path)
