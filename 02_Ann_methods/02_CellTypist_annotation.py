import os
import re
#from snakemake.script import snakemake

import numpy as np
import sklearn
import scipy as sp
from scipy.sparse import csr_matrix
import pandas as pd

import anndata
from anndata import AnnData
import scanpy as sc

import celltypist
from celltypist import models
import scvi
from scvi.model import SCVI

sc.settings.verbosity = 0
print(sc.logging.print_header())

os.makedirs(os.path.dirname(snakemake.output["annotation"]), exist_ok=True)

PATH_TO_OUTPUT_FOLDER = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/labelTransferBenchmark/01_SeuratEvaluation/05_Output"
# INPUT_FOLDER = "00_DataPreparation"
# ANNOTATIONS_FOLDER = "01_Predictions"

# DATA_FILE_NAME = "PBMC_V2_VF1_AllGenes"
# DATASET_NAME = "Dataset3"
# PATIENTS = ["GSM5584156_1", "GSM5584156_2", "GSM5584156_3"]

print("loading ref...")
print("Inputs : " + str(snakemake.input))
try:
    MAJORITY_VOTING = snakemake.params["extra_params"]["maj-vote"] == "True"
except:
    MAJORITY_VOTING = False

try:
    USE_SGD = snakemake.params["extra_params"]["use-SGD"] == "True"
except:
    USE_SGD = False  

try:
    BALANCED_TYPES = snakemake.params["extra_params"]["balanced-types"] == "True"
    USE_SGD = True
except:
    BALANCED_TYPES = False  

gene_names = pd.read_csv(snakemake.input["ref_gene_names"], sep=",", index_col=[0])

ref_counts = csr_matrix(sp.io.mmread(snakemake.input["ref_counts"], spmatrix=False))
ref_metadata = pd.read_csv(snakemake.input["ref_metadata"], sep=",", index_col=[0])
# query_counts = csr_matrix(sp.io.mmread(os.path.join(folder_name, filename+"_query_counts.mtx"), spmatrix=False)).toarray()
# query_data = sc.read_mtx(os.path.join(folder_name, filename+"_query_counts.mtx"))
# query_metadata = pd.read_csv(os.path.join(folder_name, filename+"_query_metadata.csv"), sep=",", index_col=[0])
ref = AnnData(X=ref_counts, obs=ref_metadata)


#PREPROCESSING REF
print("Preprocessing reference...")

ref.layers["counts"] = ref.X.copy()
sc.pp.normalize_total(ref, target_sum=1e4)   # all genes to 10 000 counts total in each cell
sc.pp.log1p(ref)   # log(1+x)

sc.pp.highly_variable_genes(ref, n_top_genes = 5000, flavor="cell_ranger", batch_key="orig.ident")


#TRAINING
model = celltypist.train(ref, labels = 'ident', n_jobs = 10, feature_selection = True, use_SGD=USE_SGD, mini_batch=BALANCED_TYPES, balance_cell_type=BALANCED_TYPES)
#model_path = os.path.join(PATH_TO_OUTPUT_FOLDER,"Models/02_CellTypist/model_"+snakemake.params["reference"] +"_no"+snakemake.params["noDataset"]+".pkl")
print("Saving model...")
#model.write(model_path)


query_counts = csr_matrix(sp.io.mmread(snakemake.input["query_counts"], spmatrix=False))
query_metadata = pd.read_csv(snakemake.input["query_metadata"], sep=",", index_col=[0])
query = AnnData(X=query_counts, obs=query_metadata)
print("Query len before ann : " + str(len(query)))

query.layers["counts"] = query.X.copy()
sc.pp.normalize_total(query, target_sum=1e4)
sc.pp.log1p(query)


#TESTING
predictions = celltypist.annotate(query, model=model, majority_voting=MAJORITY_VOTING)
query = predictions.to_adata()

print("Query len after ann : " + str(len(query)))

# sc.pp.neighbors(query)
# sc.tl.umap(query)
# sc.pl.umap(query, color = ['ident', 'predicted_labels'], legend_loc = 'on data')


# conf_mat = sklearn.metrics.ConfusionMatrixDisplay.from_predictions(query.obs["ident"], query.obs["predicted_labels"], normalize="true")
# conf_mat.plot()

query.obs["predicted_labels"].to_csv(os.path.join(snakemake.output["annotation"]))

