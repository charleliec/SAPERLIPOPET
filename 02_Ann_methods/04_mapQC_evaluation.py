print("Starting imports in 04_mapQC_evaluation.py...")
import os
import re
import pprint
import logging
import numpy as np
import pandas as pd
import anndata
from datetime import datetime

# Needed for the import of scanpy
date_time_launch = datetime.now().strftime("%d-%m-%Y_%Hh%Mmin%Ss")
cache_dir = f"/tmp/cache_saperlipopet_mapqc_{snakemake.params['ann_method']}_" + date_time_launch
os.environ["NUMBA_CACHE_DIR"] = os.path.join(cache_dir, "numba_cache_")
os.environ['MPLCONFIGDIR'] = os.path.join(cache_dir, 'mplconfig_cache')
os.environ["PYTORCH_KERNEL_CACHE_PATH"] = os.path.join(cache_dir, "torch_cache")
os.environ["FONTCONFIG_FILE"] = os.path.join(cache_dir, "fontconfig_cache")
os.environ["XDG_CACHE_HOME"] = cache_dir
for d in ["XDG_CACHE_HOME", "NUMBA_CACHE_DIR", "MPLCONFIGDIR", "PYTORCH_KERNEL_CACHE_PATH", "FONTCONFIG_FILE"]:
    os.makedirs(os.environ[d], exist_ok=True)
import mapqc
print("Done with imports.")

logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)

os.makedirs(os.path.dirname(snakemake.output["metadata_csv"]), exist_ok=True)

BATCH_KEY = snakemake.params["batch_key"]
STUDY_KEY = snakemake.params["study_key"]

adata = anndata.io.read_h5ad(snakemake.input["scvi_embedding"])
print("study key : ", adata.obs[snakemake.params["study_key"]].value_counts())

logger.debug("embedding_scvi \n" + str(adata.obs["orig"].unique()))
logger.debug("table sample ref_q_key : \n" + str(adata.obs.groupby([BATCH_KEY, "orig"]).size().unstack(fill_value=0)))
#logger.debug("table sample ref_q_key : \n" + str(adata.obs.groupby([BATCH_KEY, "ident"]).size().unstack(fill_value=0)))



nhood_info_df, sample_dists = mapqc.run_mapqc(
    adata=adata,
    adata_emb_loc="X_scANVI_query_umap",  # our mapped embedding is in adata.obsm["X_scANVI..."]
    ref_q_key="orig",  # .obs column with reference/query information
    q_cat="query",  # category for query
    r_cat="ref",  # category for reference
    sample_key=BATCH_KEY,  # .obs column with sample information
    n_nhoods=1000,  # number of neighborhoods to use for calculation of mapQC scores
    k_min=200,  # minimum neighborhood size
    k_max=1000,  # maximum neighborhood size
    exclude_same_study=False,
    study_key=snakemake.params["study_key"],  # .obs column with study/dataset information
    grouping_key="leiden_clustering",  # .obs column with a grouping of the data (not required)
    seed=10,  # set the seed for reproducibility
    overwrite=True,
    return_nhood_info_df=True,
    return_sample_dists_to_ref_df=True
)

print(nhood_info_df)
print("Reason for being filtered : \n", nhood_info_df["filter_info"].value_counts())
logger.debug("nhood_info_df \n" + str(nhood_info_df["filter_info"].value_counts()))
logger.debug("nhood_info_df \n" + nhood_info_df.to_string())
logger.debug("smple_dists \n" + sample_dists.to_string())
logger.debug("embedding_scvi \n" + adata.obs.to_string())


control_categories = list(adata[adata.obs["orig"]=="query"].obs[BATCH_KEY].unique())
print("Control categories : ", control_categories)
stats = mapqc.evaluate(
    adata,
    case_control_key=BATCH_KEY,  # the .obs column that includes information of the case/control status of each cell (i.e. the subject it came from)
    case_cats=[],  # a list containing all case categories in your query data
    control_cats=control_categories,  # a list containing all control categories in your query data
)

pprint.pprint(stats)
logger.debug("mapqc stats \n"+pprint.pformat(stats))
print(adata.obs)
print(adata[adata.obs["orig"]=="query"].obs)

if snakemake.params["ref_evaluation"]:
    adata[adata.obs["orig"]=="query"].obs[[STUDY_KEY, BATCH_KEY, "ident", "predictions", "mapqc_score", "mapqc_score_binary"]].to_csv(snakemake.output["metadata_csv"])
else:
    adata[adata.obs["orig"]=="query"].obs[[STUDY_KEY, BATCH_KEY, "predictions", "mapqc_score", "mapqc_score_binary"]].to_csv(snakemake.output["metadata_csv"])
print("Metadata saved in : ", snakemake.output["metadata_csv"])
















