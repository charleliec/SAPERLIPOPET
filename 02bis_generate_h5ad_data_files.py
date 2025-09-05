import os
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
import anndata as ad

os.makedirs(os.path.dirname(snakemake.output["h5ad_file"]), exist_ok=True)


gene_names = list(pd.read_csv(snakemake.input["gene_names"], sep=",", index_col=[0]).iloc[:, 0])
counts = sp.csr_matrix(scipy.io.mmread(snakemake.input["counts"], spmatrix=False))
metadata = pd.read_csv(snakemake.input["metadata"], sep=",", index_col=[0])
adata = ad.AnnData(X=counts, obs=metadata)
adata.var_names = gene_names

if snakemake.params["ref"]==False:
    query_umap = pd.read_csv(snakemake.input["umap"], sep=",", index_col=[0])
    adata.obsm["X_umap"] = np.array(query_umap)
print("Adata saved in h5ad : \n", adata)
adata.write_h5ad(filename=snakemake.output["h5ad_file"])
