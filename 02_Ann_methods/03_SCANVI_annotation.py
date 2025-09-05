print("Starting import in 03_SCANVI_annotation.py...")
import os
import re
import numpy as np
import sklearn
import scipy as sp
from scipy.sparse import csr_matrix
import pandas as pd
import anndata
from anndata import AnnData
from datetime import datetime

# Needed for the import of scanpy
date_time_launch = datetime.now().strftime("%d-%m-%Y_%Hh%Mmin%Ss")
cache_dir = f"/tmp/cache_saperlipopet_scanvi_{snakemake.params['ann_method']}_" + date_time_launch
os.environ["NUMBA_CACHE_DIR"] = os.path.join(cache_dir, "numba_cache_")
os.environ['MPLCONFIGDIR'] = os.path.join(cache_dir, 'mplconfig_cache')
os.environ["PYTORCH_KERNEL_CACHE_PATH"] = os.path.join(cache_dir, "torch_cache")
os.environ["FONTCONFIG_FILE"] = os.path.join(cache_dir, "fontconfig_cache")
os.environ["XDG_CACHE_HOME"] = cache_dir
for d in ["XDG_CACHE_HOME", "NUMBA_CACHE_DIR", "MPLCONFIGDIR", "PYTORCH_KERNEL_CACHE_PATH", "FONTCONFIG_FILE"]:
    os.makedirs(os.environ[d], exist_ok=True)
    
print("Importing scanpy...")
import scanpy as sc
import scvi
print("Done with imports.")

#scvi.settings.dl_num_workers = 79
sc.settings.verbosity = 3
os.makedirs(os.path.dirname(snakemake.output["scvi_embedding"]), exist_ok=True)

RETRAIN_MODELS = True

PATH_TO_OUTPUT_FOLDER = snakemake.params["output_dir"]

REF_UNLABELED_MAX_EPOCHS = snakemake.params["extra_params"]["r-unlab-ep"]
REF_LABELED_MAX_EPOCHS = snakemake.params["extra_params"]["r-lab-ep"]
QUERY_MAX_EPOCHS = snakemake.params["extra_params"]["q-ep"]

BATCH_KEY = snakemake.params["batch_key"]

print("loading ref...")
ref = anndata.io.read_h5ad(snakemake.input["ref_file"])
print(ref.obs["ident"].unique())
print("Ref nb genes before subsetting : ", len(ref.var_names))
sc.pp.filter_genes(ref, min_cells=100)
print("Ref nb genes after subsetting : ", len(ref.var_names))
ref_genes_set = set(ref.var_names)


print("loading query...")
query = anndata.io.read_h5ad(snakemake.input["query_file"])
query.layers["counts"] = query.X.copy()
print("Query nb genes : ", len(query.var_names))
gene_intersect = list(ref_genes_set.intersection(query.var_names))
# Subsetting query + ref
ref = ref[:, gene_intersect]
query = query[:, gene_intersect]
print("Intersect nb genes : ", len(query.var_names))
print("genes in ref but not query : ", ref_genes_set.difference(query.var_names))





#PREPROCESSING REF
print("Preprocessing reference...")

ref.layers["counts"] = ref.X.copy()
sc.pp.normalize_total(ref, target_sum=1e4)   # all genes to 10 000 counts total in each cell
sc.pp.log1p(ref)   # log(1+x)

sc.pp.highly_variable_genes(ref, n_top_genes = 5000, flavor="cell_ranger", batch_key=BATCH_KEY)


#TRAINING REFERENCE
ref.obs["scanvi_label"] = ref.obs["ident"].copy()
scvi.model.SCVI.setup_anndata(ref, batch_key=BATCH_KEY, labels_key="scanvi_label", layer="counts")
scvi_model_path = os.path.join(PATH_TO_OUTPUT_FOLDER,"Models/03_SCANVI/scvi_model_"+snakemake.params["ref_name"] +"_no"+f"_ref_unlab_max_ep={REF_UNLABELED_MAX_EPOCHS}")

try:
    if RETRAIN_MODELS:
        raise Exception("Retrain set to true.")
    scvi_ref_model = scvi.model.SCVI.load(scvi_model_path, adata=ref)
    print("SCVI Model Loaded from existing model")
except Exception as e:
    print("Error while loading SCVI model : ", e, "\n Training from 0 instead...")
    scvi_ref_model = scvi.model.SCVI(ref, n_layers=2, 
        encode_covariates=True,
        deeply_inject_covariates=False,
        use_layer_norm="both",
        use_batch_norm="none")
    scvi_ref_model.train(max_epochs=REF_UNLABELED_MAX_EPOCHS)
    scvi_ref_model.save(scvi_model_path, overwrite=True)

scvi_ref_model.view_anndata_setup()



#TRAINING REFERENCE ANNOTATED
scanvi_model_path = os.path.join(PATH_TO_OUTPUT_FOLDER,"Models/03_SCANVI/scanvi_model_"+snakemake.params["ref_name"] +"_no"+f"_ref_unlab_max_ep={REF_UNLABELED_MAX_EPOCHS}_ref_lab_max_ep={REF_LABELED_MAX_EPOCHS}")
try:
    if RETRAIN_MODELS:
        raise Exception("Retrain set to true.")
    scanvi_ref_model = scvi.model.SCANVI.load(scanvi_model_path, adata=ref)
    print("SCANVI Model Loaded from existing model")
except Exception as e:
    print("Error while loading SCANVI model : ", e, "\n Training from 0 instead...")
    scanvi_ref_model = scvi.model.SCANVI.from_scvi_model(
        scvi_ref_model,
        adata = ref,
        unlabeled_category="Unknown",
        labels_key="scanvi_label",
    )
    scanvi_ref_model.train(max_epochs=REF_LABELED_MAX_EPOCHS, n_samples_per_label=500)
    scanvi_ref_model.save(scanvi_model_path, overwrite=True)

scanvi_ref_model.view_anndata_setup()



query.obs["scanvi_label"] = "Unknown"
query.obs["orig"] = "query"
ref.obs["orig"] = "ref"
adata_full = anndata.concat([ref, query])

# ANNOTATING QUERY
#adata_full.obs = adata_full.obs[["scanvi_label", BATCH_KEY]]
scvi.model.SCANVI.setup_anndata(adata_full, layer="counts", batch_key=BATCH_KEY, labels_key="scanvi_label", unlabeled_category="Unknown")
query_model = scvi.model.SCANVI.load_query_data(
    adata_full,
    scanvi_ref_model,
    freeze_dropout = False,
)
query_model._unlabeled_indices = np.nonzero(adata_full.obs["scanvi_label"]=="Unknown")
query_model._labeled_indices = np.nonzero(adata_full.obs["scanvi_label"]!="Unknown")

query_model.train(max_epochs=QUERY_MAX_EPOCHS, plan_kwargs=dict(weight_decay=0.0))

# Saving integration with scanvi + query training
adata_full.obsm["X_scANVI_query"] = query_model.get_latent_representation(adata_full)
sc.pp.neighbors(adata_full, use_rep="X_scANVI_query", key_added="scANVI_query_neighbors")
sc.tl.umap(adata_full, key_added="X_scANVI_query_umap", neighbors_key="scANVI_query_neighbors")
sc.tl.leiden(adata_full, key_added="leiden_clustering", neighbors_key="scANVI_query_neighbors")

# PRED PROBA + LABELS
query_mask = adata_full.obs['orig'] == "query"

adata_full.obs['predictions'] = None
predictions = query_model.predict(query)
adata_full.obs.loc[query_mask, 'predictions'] = predictions

probas = query_model.predict(query, soft=True)
for col in probas.columns:
    adata_full.obs[col + '_pred_proba'] = np.nan
    adata_full.obs.loc[query_mask, col + '_pred_proba'] = probas[col].values

adata_full.write_h5ad(filename=snakemake.output["scvi_embedding"])
print("Adata saved at : ", snakemake.output["scvi_embedding"])















