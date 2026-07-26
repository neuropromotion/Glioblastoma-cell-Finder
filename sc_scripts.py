import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from constants import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def annotate_qc(adata: ad.AnnData) -> None:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

def apply_qc_filters(
    adata: ad.AnnData,
    verbose: bool = False,
) -> ad.AnnData:
    before = adata.n_obs
    genes_before = adata.n_vars
    adata = adata[adata.obs["n_genes_by_counts"] >= MIN_GENES].copy()

    max_genes = np.quantile(
        adata.obs["n_genes_by_counts"],
        GENE_UPPER_QUANTILE,
    )
    adata = adata[adata.obs["n_genes_by_counts"] <= max_genes].copy()

    max_counts = np.quantile(
        adata.obs["total_counts"],
        COUNT_UPPER_QUANTILE,
    )
    adata = adata[adata.obs["total_counts"] <= max_counts].copy()
    adata = adata[adata.obs["pct_counts_mt"] < MT_THRESHOLD].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

    if verbose:
        print(
            f"Cells: {before} -> {adata.n_obs} | "
            f"Genes: {genes_before} -> {adata.n_vars}\n"
            f"{GENE_UPPER_QUANTILE}th percentile n_genes = {max_genes:.0f}\n"
            f"{COUNT_UPPER_QUANTILE}th percentile total_counts = {max_counts:.0f}"
        )

    return adata

def main_pipeline(adata, n_components=15, resolution=0.1):
    adata = adata.copy() # Копия нужна, чтобы объект можно было менять без ошибок
    adata.layers["counts"] = adata.X.copy()
    # -----------------------------------------Normalization-----------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)  # CPM normalization
    sc.pp.log1p(adata)  # Log-transform

    # -----------------------------------------FindVariableFeatures-----------------------------------------------
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor='seurat',
    )
    # -----------------------------------------Scaling-----------------------------------------------
    sc.pp.scale(adata, max_value=10)  # Z-score normalization
    # -----------------------------------------PCA-----------------------------------------------
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_components)
    # -----------------------------------------Clustering-----------------------------------------------
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_components)
    sc.tl.leiden(adata, resolution=resolution, key_added='clusters')
    # -----------------------------------------UMAP-----------------------------------------------
    sc.tl.umap(adata, n_components=2)
    adata.X = adata.layers["counts"].copy()
    return adata


