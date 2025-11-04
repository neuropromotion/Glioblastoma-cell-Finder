Glioblastoma Cell Finder (R, Gradient Boosting)
![Status](https://img.shields.io/badge/status-active-brightgreenhttps://img.shields.io/badge[R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo= boosting model for annotating glioblastoma vs stromal cells directly on Seurat objects in R. The repository exposes a single entry point that returns the same Seurat object with two additional metadata columns for categorical labels and class probabilities.

Features
Drop‑in inference on Seurat objects with a single function call.

Adds:

ML_annotation: discrete class labels for plotting and downstream filtering.

GBM_prob: glioblastoma class probability in for continuous visualizations.

Compatible with standard Seurat plotting (DimPlot, FeaturePlot).

Installation
Ensure R (≥ 4.3) is installed.

Install required packages:

r
install.packages(c("Seurat", "xgboost"))
# optionally:
# install.packages("Matrix")
Clone the repository and make sure model.R is available in your working directory.

Quick Start
1) Load the function
r
source("model.R")  # exposes predict_and_add_metadata()
2) Run inference on a Seurat object
Inputs:

seurat.object: a Seurat object with an RNA assay and count layer (or equivalent features used during training).

model_path: path to the saved model weights (e.g., "boosting/xgboost.model").

Output:

The same Seurat object with new metadata fields:

ML_annotation: character/factor labels (e.g., "Glioblastoma cells", "Stromal cells").

GBM_prob: numeric probability for the glioblastoma class.

r
library(Seurat)

seu <- readRDS("path/to/seurat_object.rds")
seu <- predict_and_add_metadata(
  seurat.object = seu,
  model_path    = "boosting/xgboost.model"
)
3) Visualize results
Discrete class labels:

r
DimPlot(seu, group.by = "ML_annotation", label = FALSE)
Probability heatmap:

r
FeaturePlot(seu, features = "GBM_prob")
Function Contract
File: model.R

Function: predict_and_add_metadata(seurat.object, model_path)

Behavior:

Loads the trained gradient boosting model from model_path.

Computes per‑cell probabilities.

Writes two metadata columns:

ML_annotation

GBM_prob

Returns the updated Seurat object.

Data and Assumptions
The Seurat object should contain the RNA assay and the same type of features used during training (e.g., chromosome‑level summaries) with identical column order.

If feature construction is part of your pipeline, ensure the same preprocessing steps and feature ordering are applied before inference.

Troubleshooting
Mismatch in feature columns:

Ensure the inference feature construction strictly matches the training setup (same set and order of columns).

Empty or missing metadata:

Confirm that the model_path points to a valid trained model file.

Check that the object contains the assay/layer expected by the feature builder.

Suggested Project Structure
text
.
├── model.R                 # predict_and_add_metadata implementation
├── boosting/
│   └── xgboost.model       # trained model weights (example path)
├── assets/
│   └── teaser.png          # optional visuals for README
└── README.md
Example Snippet (End‑to‑End)
r
library(Seurat)
source("model.R")

seu <- readRDS("data/seurat_gbm.rds")
seu <- predict_and_add_metadata(seurat.object = seu, model_path = "boosting/xgboost.model")

DimPlot(seu, group.by = "ML_annotation", label = FALSE)
FeaturePlot(seu, features = "GBM_prob")
License
MIT. See LICENSE for details.

Acknowledgements
Seurat for the single‑cell analysis framework.

XGBoost for the gradient boosting library.
