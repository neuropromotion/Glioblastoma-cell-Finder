import numpy as np
import pandas as pd
import scipy.sparse as sp
import xgboost as xgb 

class GCF:
    def __init__(self, weights_path='/mnt/jack-5/amismailov/CAF_study/CNV/GCF_config.json'): 
        

        self.model = xgb.XGBClassifier()
        self.model.load_model(weights_path)
        
    def build_log2_cpm1000_plus1_from_adata(self, adata_or_X, genes=None, cells=None):
        is_anndata = hasattr(adata_or_X, "X") and hasattr(adata_or_X, "var_names") and hasattr(adata_or_X, "obs_names")  # [web:41]
        if is_anndata:
            adata = adata_or_X
            X = adata.layers['counts'] if ('counts' in adata.layers) else adata.X  # [web:39][web:42]
            genes = adata.var_names.to_numpy()
            cells = adata.obs_names.to_numpy()
        else:
            if genes is None or cells is None:
                raise ValueError("При передаче матрицы X необходимо указать genes и cells.")  # [web:41]
            X = adata_or_X

        if sp.issparse(X):
            Xg = X.T.tocsr()
            row_sums = np.asarray(Xg.sum(axis=1)).ravel()
            keep = row_sums > 0
            if keep.sum() != len(keep):
                Xg = Xg[keep]; genes = np.asarray(genes)[keep]
            col_sums = np.asarray(Xg.sum(axis=0)).ravel()
            col_sums[col_sums == 0] = 1.0
            Xg = Xg @ sp.diags(1000.0 / col_sums)
            Xg = Xg.tocoo(); Xg.data = np.log2(Xg.data + 1); Xg = Xg.tocsr()
            gene_sums = np.asarray(Xg.sum(axis=1)).ravel()
            keep2 = gene_sums >= 100.0
            if keep2.sum() != len(keep2):
                Xg = Xg[keep2]; genes = np.asarray(genes)[keep2]
        else:
            Xg = np.asarray(X.T)
            row_sums = Xg.sum(axis=1)
            keep = row_sums > 0
            if keep.sum() != len(keep):
                Xg = Xg[keep]; genes = np.asarray(genes)[keep]
            col_sums = Xg.sum(axis=0); col_sums[col_sums == 0] = 1.0
            Xg = Xg / col_sums * 1000.0
            Xg = np.log2(Xg + 1)
            gene_sums = Xg.sum(axis=1)
            keep2 = gene_sums >= 100.0
            if keep2.sum() != len(keep2):
                Xg = Xg[keep2]; genes = np.asarray(genes)[keep2]

        return Xg, np.asarray(genes), np.asarray(cells)  # [web:41]

    def chromosome_means_from_logmat(self, Xg, genes, cells, map_path=None):
        map_path = map_path or self.map_path
        m = pd.read_csv(map_path, sep='\t').rename(
            columns={'HGNC symbol':'symbol', 'Chromosome/scaffold name':'chr'}
        )  
        m = m.dropna(subset=['symbol','chr'])
        m['symbol'] = m['symbol'].astype(str).str.replace(r'\.\d+$','', regex=True).str.upper()
        m['chr'] = m['chr'].astype(str).str.replace('chr','', regex=False)
        m = m[m['chr'].isin([str(i) for i in range(1,23)])]

        genes = np.asarray(genes)
        mask_hla = ~pd.Series(genes).str.startswith('HLA', na=False).to_numpy()
        Xg1 = Xg[mask_hla] if sp.issparse(Xg) else Xg[mask_hla]
        genes1 = genes[mask_hla]

        gindex = pd.Index(genes1)
        chr_means = {}
        for i in range(1, 23):
            g = m.loc[m['chr'] == str(i), 'symbol'].unique()
            g_in = gindex.intersection(g)
            if len(g_in) == 0:
                chr_means[f'Chr{i}'] = np.full(len(cells), np.nan)
                continue
            rows = gindex.get_indexer(g_in)
            if sp.issparse(Xg1):
                vals = np.asarray(Xg1[rows].mean(axis=0)).ravel()
            else:
                vals = Xg1[rows].mean(axis=0)
            chr_means[f'Chr{i}'] = vals

        return pd.DataFrame(chr_means, index=cells)  # [web:41]

    def get_means(self, obj, path):
        Xg, genes, cells = self.build_log2_cpm1000_plus1_from_adata(obj)   
        return self.chromosome_means_from_logmat(Xg, genes, cells, path)
    
    def predict(self, obj, map_path='/home/amismailov/mart_export.txt'): 
        adata = obj.copy() 
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()   
            
        data = self.get_means(adata, map_path)
        data["proba"] = self.model.predict_proba(data)[:, 1]
        data["pred"] = (data["proba"] >= 0.5).astype(int)

        label = ['Glioblastoma cells' if x == 0 else 'Stromal cells' for x in data['pred'].tolist()]
        data['label'] = label
        return data['label'].to_frame()

    def require_dataframe(self, obj, name='obj'):
        if not isinstance(obj, pd.DataFrame):
            raise TypeError(f"{name} has to be pandas.DataFrame, instead it has {type(obj).__name__} type")
        if len(obj.columns) != 22:
            raise TypeError(f"{name} has to have 22 columns (22 chromosomes), instead it has {len(obj.columns)} columns")
    
    def predict_2(self, data):
        self.require_dataframe(data)
        
        data["proba"] = self.model.predict_proba(data)[:, 1]
        data["pred"] = (data["proba"] >= 0.5).astype(int)
        label = ['Glioblastoma cells' if x == 0 else 'Stromal cells' for x in data['pred'].tolist()]
        data['label'] = label
        return data['label'].to_frame()




