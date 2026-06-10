import anndata as ad
import scanpy as sc

adata = ad.read_h5ad(par["input"])
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(par["output"])
