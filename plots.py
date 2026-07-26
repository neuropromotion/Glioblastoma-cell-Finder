import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from constants import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde
import plotly.graph_objects as go


def plot_pca(chr_means):
    # клетки × хромосомы
    X_chr = chr_means.T

    # масштабирование признаков (хромосом)
    X_scaled = StandardScaler().fit_transform(X_chr)

    # PCA 2D
    pca = PCA(n_components=2)

    coords = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(
        coords,
        index=X_chr.index,
        columns=["PC1", "PC2"]
    )



    plt.figure(figsize=(7,6))

    plt.scatter(
        pca_df["PC1"],
        pca_df["PC2"],
        s=10
    )

    plt.xlabel(
        f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
    )

    plt.ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
    )

    plt.title("Chromosome composition PCA")
    plt.show()

def plot_3d_kde(filtered, feat_x, feat_y):
    
    x = filtered.loc[feat_x].values
    y = filtered.loc[feat_y].values

    # KDE
    kde = gaussian_kde(np.vstack([x, y]))

    # Grid
    xx, yy = np.mgrid[
        x.min():x.max():120j,
        y.min():y.max():120j
    ]

    positions = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(positions).reshape(xx.shape)

    fig = go.Figure(
        go.Surface(
            x=xx,
            y=yy,
            z=zz,
            colorscale="Viridis",
            showscale=True
        )
    )

    fig.update_layout(
        title=f"3D KDE: {feat_x} vs {feat_y}",
        scene=dict(
            xaxis_title=feat_x,
            yaxis_title=feat_y,
            zaxis_title="Density"
        )
    )

    fig.show()


'''
import numpy as np
import plotly.graph_objects as go
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata


def plot_3d_kde(filtered, x_gene, y_gene, color_feature):

    x = filtered.loc[x_gene].values
    y = filtered.loc[y_gene].values
    c = filtered.loc[color_feature].values

    # KDE
    kde = gaussian_kde(np.vstack([x, y]))

    # Grid
    xx, yy = np.mgrid[
        x.min():x.max():120j,
        y.min():y.max():120j
    ]

    positions = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(positions).reshape(xx.shape)

    # Интерполяция признака на сетку
    cc = griddata(
        points=np.column_stack([x, y]),
        values=c,
        xi=(xx, yy),
        method="linear"
    )

    # если за пределами convex hull появились NaN
    if np.isnan(cc).any():
        cc_nearest = griddata(
            np.column_stack([x, y]),
            c,
            (xx, yy),
            method="nearest"
        )
        cc = np.where(np.isnan(cc), cc_nearest, cc)

    fig = go.Figure(
        go.Surface(
            x=xx,
            y=yy,
            z=zz,
            surfacecolor=cc,
            colorscale="Plasma",
            colorbar=dict(title=color_feature),
            showscale=True
        )
    )

    fig.update_layout(
        title=f"3D KDE: {x_gene} vs {y_gene}",
        scene=dict(
            xaxis_title=x_gene,
            yaxis_title=y_gene,
            zaxis_title="Density"
        )
    )

    fig.show()

'''

def plot_3d_hist(filtered, feat_x, feat_y):


    x = filtered.loc[feat_x].values
    y = filtered.loc[feat_y].values

    H, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=25
    )

    xc = (xedges[:-1] + xedges[1:]) / 2
    yc = (yedges[:-1] + yedges[1:]) / 2

    X, Y = np.meshgrid(xc, yc)

    fig = go.Figure(
        go.Surface(
            x=X,
            y=Y,
            z=H.T,
            colorscale="Turbo"
        )
    )

    fig.update_layout(
        title=f"3D Histogram: {feat_x} vs {feat_y}",
        scene=dict(
            xaxis_title=feat_x,
            yaxis_title=feat_y,
            zaxis_title="Cell count"
        )
    )

    fig.show()