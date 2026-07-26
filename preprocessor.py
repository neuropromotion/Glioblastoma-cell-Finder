import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from constants import *
from sc_scripts import *

class Preprocessor: 
    NORM_MODES = ( 
        "cpm_shares",          # CPM → map → shares
        "cpm_agg_log_shares",  # CPM → map → log1p(features) → shares 
        "cpm_log_shares",      # CPM → log1p(genes) → map → shares
        "cpm_log",             # CPM → log1p(genes) → map
        "cpm",                 # CPM → map
    )

    def __init__(self, chr_mapping_path='gene_mapping.txt',
                       ens_map_path='ensembl_gene_mapping.csv'): 
        self.ens_map_path = ens_map_path
        self.ens_map = pd.read_csv(ens_map_path)
        
    def load_10x_data(self, path_10x): 
        path = Path(path_10x)
        features = path / "features.tsv.gz"
        if features.is_file():
            with gzip.open(features, "rt") as f:
                n_cols = len(next(f).rstrip("\n").split("\t"))
            if n_cols < 3:
                return self._read_10x_mtx_two_col_features(path)
        return sc.read_10x_mtx(path)

    @staticmethod
    def _read_10x_mtx_two_col_features(path: Path) -> ad.AnnData:
        adata = sc.read_mtx(path / "matrix.mtx.gz").T
        genes = pd.read_csv(path / "features.tsv.gz", header=None, sep="\t")
        var_names = ad.utils.make_index_unique(pd.Index(genes[1].astype(str).values))
        adata.var_names = var_names
        adata.var["gene_ids"] = genes[0].values
        barcodes = pd.read_csv(path / "barcodes.tsv.gz", header=None)
        adata.obs_names = barcodes[0].astype(str).values
        return adata

    def _log1p(self, df: pd.DataFrame) -> pd.DataFrame:
        arr = df.to_numpy(copy=True)
        arr = np.log1p(arr)
        return pd.DataFrame(arr, index=df.index, columns=df.columns)
    
    def load_txt_data(self, path_smartseq: str, sample_name = 'sample_name'):
        path_smartseq = str(path_smartseq)
        if path_smartseq.endswith('.txt.gz'):
            df = pd.read_table(path_smartseq, sep="\t", index_col=0)
        elif path_smartseq.endswith('.csv'):
            df = pd.read_csv(path_smartseq, index_col=0, nrows=15000).T
        else:
            raise ValueError(f"Unsupported file extension: {path_smartseq}")
        adata = sc.AnnData(df.T)
        adata.obs["sample"] = sample_name
        return adata
    
    def qc_filter(self, adata: ad.AnnData, verbose = False) -> ad.AnnData:
        annotate_qc(adata)
        adata = apply_qc_filters(adata, verbose = verbose)
        sc.pp.scrublet(adata)
        #print(f'Detected {adata.obs["predicted_doublet"].sum()} doublets')
        adata = adata[
            ~adata.obs["predicted_doublet"]
        ].copy()
        return adata
    
    def get_counts_matrix(self, adata: ad.AnnData) -> np.ndarray:
        df = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.obs_names,      
            columns=adata.var_names,    
        )
        return df

    def CPM(self, df: pd.DataFrame) -> pd.DataFrame:
        lib_size = df.sum(axis=1)
        cpm = df.div(lib_size, axis=0) * 1_000_000
        return cpm

    def replace_genes_names(
        self,
        data,
    ): 
        self.ens_map = self.ens_map.dropna(subset=["feature_name", "feature_id"])
        self.ens_map = self.ens_map.drop_duplicates(subset=["feature_name"])
    
        # dict: HGNC → ENSG
        ens_map = dict(zip(self.ens_map["feature_name"], self.ens_map["feature_id"]))
    
        df = data.T.copy(deep=True)
    
        # оригинальные имена генов
        #df["original_symbol"] = df.index
    
        total_genes = df.shape[0]
        #print(f"Replacing {total_genes} genes...")
    
        # маппинг index → ENSG
        df["ensembl_id"] = df.index.map(ens_map)
    
        # сколько найдено
        found = df["ensembl_id"].notna().sum()
        percent = found / total_genes * 100
    
        # удалить ненайденные
        df_clean = df.dropna(subset=["ensembl_id"]).copy()
    
        # убрать дубликаты ENSG
        df_clean = df_clean[~df_clean["ensembl_id"].duplicated()]
    
        # ENSG → индекс
        df_clean.index = df_clean.pop("ensembl_id")
        df_clean = df_clean.sort_index()
    
        #print(f"✔ Found ENSG for {found}/{total_genes} genes ({percent:.2f}%)")
        #print(f"✔ After removing duplicates: {len(df_clean)} unique ENSG")
    
        return df_clean

    def remove_HLA_genes(self, df: pd.DataFrame) -> pd.DataFrame:
        before = df.shape[1]
        df = df.loc[:, ~df.columns.str.contains("HLA-")]
        after = df.shape[1]
        #print(f"Removed {before - after} HLA genes")
        return df

    def map_genes_to_chromosomes_22(self, df: pd.DataFrame) -> pd.DataFrame:
        mapper = self.chr_mapping[
            self.chr_mapping["gene"].isin(df.index)
        ][["gene", "chromosome"]].copy()

        joined = df.join(mapper.set_index("gene"), how="inner")
        return joined.groupby("chromosome").sum(numeric_only=True)

    def add_chr7_10_derived_features(self, feat: pd.DataFrame) -> pd.DataFrame: 
        out = feat.copy()
        if "7" not in out.index:
            parts7 = [p for p in ("7p", "7q") if p in out.index]
            if not parts7:
                raise KeyError("Need '7' or '7p'/'7q'")
            out.loc["7"] = out.loc[parts7].sum(axis=0)
        if "10" not in out.index:
            parts10 = [p for p in ("10p", "10q") if p in out.index]
            if not parts10:
                raise KeyError("Need '10' or '10p'/'10q'")
            out.loc["10"] = out.loc[parts10].sum(axis=0)

        out.loc["diff_7_minus_10"] = out.loc["7"] - out.loc["10"]
        return out

    def normalize_chr_expression(self, chr_expression: pd.DataFrame) -> pd.DataFrame: 
        return chr_expression.div(chr_expression.sum(axis=0), axis=1)

    def normalize_features(self, feat: pd.DataFrame, mode: str = "cpm_shares") -> pd.DataFrame:
        if mode not in self.NORM_MODES:
            raise ValueError(f"Unknown mode {mode!r}. Choose from {self.NORM_MODES}")

        out = feat.copy()
        if mode in {"cpm_shares", "cpm_log_shares"}:
            return self.normalize_chr_expression(out)
        if mode == "cpm_agg_log_shares":
            return self.normalize_chr_expression(self._log1p(out))
        if mode in {"cpm", "cpm_log"}:
            return out
        raise ValueError(mode)


    def run_pipeline(
        self,
        path,
        _10x: bool = False,
        use_cpm: bool = True,
        mode: str = "pairwise",
        log: bool = False,
        add_qc: bool = True,
        drop_chroms: tuple[str, ...] = tuple(),
        verbose: bool = False,
    ) -> pd.DataFrame:
        if mode not in {"pairwise", "vs_background", "all_shares"}:
            raise ValueError(
                "mode must be one of: 'pairwise', 'vs_background', 'all_shares'"
            )

        if _10x:
            adata = self.load_10x_data(path)
        else:
            adata = self.load_txt_data(path)
        adata = self.qc_filter(adata, verbose=verbose)

        counts = self.remove_HLA_genes(self.get_counts_matrix(adata))
        mat_cells = self.CPM(counts) if use_cpm else counts
        genes_x_cells = self.replace_genes_names(mat_cells)

        mapper = self.chr_mapping[
            self.chr_mapping["gene"].isin(genes_x_cells.index)
        ][["gene", "chromosome"]].drop_duplicates(subset=["gene"]).copy()
        mapper["chromosome"] = mapper["chromosome"].astype(str)

        # effective n genes that actually mapped per chromosome
        n_genes = mapper.groupby("chromosome").size().astype(float)
        n_genes = n_genes[n_genes > 0]

        joined = genes_x_cells.join(mapper.set_index("gene"), how="inner")
        chrom_sum = joined.groupby("chromosome").sum(numeric_only=True)
        # length-normalize: sum / n_mapped_genes  (== mean over mapped genes)
        chrom_level = chrom_sum.div(n_genes, axis=0)

        keep = [
            str(i)
            for i in range(1, 23)
            if str(i) not in set(map(str, drop_chroms)) and str(i) in chrom_level.index
        ]
        chrom_level = chrom_level.reindex(keep)

        if chrom_level.empty:
            raise ValueError("No chromosomes left after drop_chroms / mapping.")

        if log:
            chrom_level = self._log1p(chrom_level)

        if mode == "pairwise":
            if not {"7", "10"}.issubset(chrom_level.index):
                raise ValueError("Need chromosomes 7 and 10 for mode='pairwise'")
            pair = chrom_level.loc[["7", "10"]]
            total = pair.sum(axis=0).replace(0, np.nan)
            out = pair.div(total, axis=1)
            out.loc["diff_7_minus_10"] = out.loc["7"] - out.loc["10"]
            out.loc['1'] = chrom_level.loc['1']

        elif mode == "vs_background":
            if not {"7", "10"}.issubset(chrom_level.index):
                raise ValueError("Need chromosomes 7 and 10 for mode='vs_background'")
            others = chrom_level.drop(index=["7", "10"])
            if others.empty:
                raise ValueError("No background chromosomes left for mode='vs_background'")
            bg = others.mean(axis=0).replace(0, np.nan)
            out = pd.DataFrame(
                [chrom_level.loc["7"] / bg, chrom_level.loc["10"] / bg],
                index=["7", "10"],
            )

        else:  # all_shares
            out = self.normalize_chr_expression(chrom_level)

        if add_qc:
            out.loc["n_genes"] = adata.obs.loc[out.columns, "n_genes_by_counts"]
            out.loc["total_counts"] = adata.obs.loc[out.columns, "total_counts"]
            out.loc["pct_mt"] = adata.obs.loc[out.columns, "pct_counts_mt"]
            out.loc["n_genes_norm"] = (
                out.loc["n_genes"] / out.loc["n_genes"].max()
            )
            out.loc["total_counts_norm"] = (
                out.loc["total_counts"] / out.loc["total_counts"].max()
            )
            out.drop(["n_genes", "total_counts"], axis=0, inplace=True)

        if verbose:
            print(
                f"v3 | use_cpm={use_cpm} log={log} mode={mode} | "
                f"drop={drop_chroms} | chroms={list(out.index)} | "
                f"cells={out.shape[1]}"
            )
            if "7" in out.index:
                print(
                    f"  7: mean={out.loc['7'].mean():.4f} std={out.loc['7'].std():.4f}"
                )
            if "10" in out.index:
                print(
                    f"  10: mean={out.loc['10'].mean():.4f} std={out.loc['10'].std():.4f}"
                )
        return out

 