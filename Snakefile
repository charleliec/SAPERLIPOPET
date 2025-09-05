import yaml
import os
from itertools import product
from random import sample
from datetime import datetime
import pprint

##################################################################################################################################################################
#
#               RUN SNAKEMAKE ON GPU !
#
#   If no venv created : 
#       virtualenv -p /usr/bin/python3 ~/.snakemake
#       source ~/.snakemake/bin/activate
#       pip install snakemake
#       pip install pulp==2.7.0
#
#
#   -> source ~/.snakemake/bin/activate 
#   -> snakemake --snakefile Snakefile --verbose --use-singularity --singularity-args="-B /mnt/DOSI:/mnt/DOSI --nv" --cores 4
#
#   At the end : deactivate
#
##################################################################################################################################################################


SCRIPTS_DIR = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/SAPERLIPOPET/SAPERLIPOPET"
SINGULARITY_DIR = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/SAPERLIPOPET/01_Containers"

PYTHON_SINGULARITY = os.path.join(SINGULARITY_DIR, "25-08-27_pyth-ann_quarto1.7_mapqc.sif")
R_SINGULARITY = os.path.join(SINGULARITY_DIR, "01-09-25-Rstud.sif")

print()
with open('./01_CONFIG.yaml', 'r') as f:
    config = yaml.safe_load(f)
f.close()

REFERENCE_EVALUATION = config["reference_evaluation"]

REF_PATH = config["reference_path"]
QUERY_PATH = config["query_path"]
OUTPUT_DIR = config["output_dir"]

BATCH_KEY = config["batch_key"]
STUDY_KEY = config["study_key"]
MODELS_PARAMS = config["models_params"]
UMAP_NAME = config["umap_name"]

# HEALTHY_SAMPLES = config["healthy_control_samples"]
# CASE_SAMPLES = config["case_samples"]

REF_NAME = os.path.basename(REF_PATH).rsplit(".",1)[0]
QUERY_NAME = os.path.basename(QUERY_PATH).rsplit(".",1)[0]

assert (os.path.basename(QUERY_PATH).rsplit(".",1)[1] in ["rds", "h5ad"]), 'Extension of query must be in [".rds", ".h5ad"]'
assert (os.path.basename(REF_PATH).rsplit(".",1)[1] in ["rds", "h5ad"]), 'Extension of reference must be in [".rds", ".h5ad"]'

QUERY_IS_RDS = os.path.basename(QUERY_PATH).rsplit(".",1)[1] == "rds"
REF_IS_RDS = os.path.basename(REF_PATH).rsplit(".",1)[1] == "rds"

LOGS_PATH = os.path.join(OUTPUT_DIR, "logs")
os.makedirs(LOGS_PATH, exist_ok=True)

print("References : ", REF_NAME)
print("Query : ", QUERY_NAME)
print("Annotations output dir : ", OUTPUT_DIR)
print("Models params : ", MODELS_PARAMS)


def generate_ann_params(base_method_name, params_dic, prod=True):
    ann_params = {}
    ann_names = []
    param_names = params_dic.keys()
    param_list = product(*params_dic.values()) if prod else zip(*params_dic.values())
    for parameters in param_list:
        par_names = ",".join([par_name+"="+str(param_value) for par_name, param_value in zip(param_names, parameters)])
        full_ann_name = base_method_name+":"+par_names+":"
        dic_ann_params = {par_name: param_value for par_name, param_value in zip(param_names, parameters)}
        ann_params[full_ann_name] = dic_ann_params
    return ann_params


prod = MODELS_PARAMS.pop("prod")
MODELS_PARAMS = {"r-unlab-ep": MODELS_PARAMS["ref_unlabeled_num_epochs"], "r-lab-ep": MODELS_PARAMS["ref_labeled_num_epochs"], "q-ep": MODELS_PARAMS["query_num_epochs"]}
scanvi_ann_params = generate_ann_params("03_SCANVI", MODELS_PARAMS, prod=prod)
#scanvi_ann_params = generate_ann_params("03_SCANVI", {"q-ep":[100,50, 100, 50], "r-unlab-ep":[100, 20, 200, 100], "r-lab-ep":[20, 20, 20, 10]}, prod=False)  # {"q-ep":[2, 5, 20, 80], "r-unlab-ep":[2, 5, 20, 80]}
#scanvi_ann_params = generate_ann_params("03_SCANVI", {"q-ep":[70, 100, 150], "r-unlab-ep":[50], "r-lab-ep":[15, 20, 25]}, prod=True)  # {"q-ep":[2, 5, 20, 80], "r-unlab-ep":[2, 5, 20, 80]}
#scanvi_ann_params = generate_ann_params("03_SCANVI", {"q-ep":[1], "r-unlab-ep":[1], "r-lab-ep":[1]}, prod=True)  # {"q-ep":[2, 5, 20, 80], "r-unlab-ep":[2, 5, 20, 80]}


#PY_ANN_METHODS = ["02_Celltypist:n=1:", "02_Celltypist:n=2:", "02_Celltypist:n=3:", "02_Celltypist:n=4:"] #+ list(scanvi_ann_params.keys()) 
ANN_METHODS = list(scanvi_ann_params.keys()) #+ ["02_Celltypist"]   #["02_Celltypist:maj-vote=True:", "02_Celltypist:maj-vote=False:", "02_Celltypist:maj-vote=False,use-sgd=True:"]
print("Annotation methods : ", ANN_METHODS)
    

ANN_METHODS_PARAMS = scanvi_ann_params
#ANN_METHODS_PARAMS = {"01_Seurat": {}, **scanvi_ann_params}


date_time_launch = datetime.now().strftime("%d-%m-%Y_%Hh%Mmin%Ss")


rule all:
    input:
        reports = os.path.join(OUTPUT_DIR, f"ann_tools_benchmark_on_ref={REF_NAME}_query={QUERY_NAME}" + date_time_launch + ".html"),
        annotations = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}__all_models_predictions_and_mapQC_scores.csv")

rule copy_h5ad_ref_file:
    input:
        rds_file = REF_PATH.rsplit(".",1)[0]+".h5ad"
    output:
        h5ad_file = "not_a_file_path__just_random_output_to_ignore" if REF_IS_RDS else os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "anndata.h5ad")
    shell:
        "cp {input} {output}"

rule copy_h5ad_query_file:
    input:
        rds_file = QUERY_PATH.rsplit(".",1)[0]+".h5ad"
    output:
        h5ad_file = "not_a_file_path__just_random_output_to_ignore" if QUERY_IS_RDS else os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "anndata.h5ad")
    shell:
        "cp {input} {output}"


rule generate_intermediate_ref_files:
    input:
        rds_file = REF_PATH.rsplit(".",1)[0]+".rds"
    output:
        gene_names = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "genes_names.csv"),
        counts = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "counts.mtx"),
        metadata = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "metadata.csv")
    params:
        ref_name = REF_NAME,
        ref = "True"
    singularity: R_SINGULARITY
    script: os.path.join(SCRIPTS_DIR, "02_generate_intermediate_data_files.R")

rule generate_intermediate_query_files:
    input:
        rds_file = QUERY_PATH.rsplit(".",1)[0]+".rds"
    output:
        gene_names = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "genes_names.csv"),
        counts = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "counts.mtx"),
        metadata = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "metadata.csv"),
        umap = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "umap.csv")
    params:
        query_name = QUERY_NAME,
        ref = "False",
        umap_name = UMAP_NAME
    singularity: R_SINGULARITY
    script: os.path.join(SCRIPTS_DIR, "02_generate_intermediate_data_files.R")

rule generate_h5ad_ref_files:
    input:
        gene_names = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "genes_names.csv"),
        counts = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "counts.mtx"),
        metadata = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "metadata.csv")
    output:
        h5ad_file = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "anndata.h5ad") if REF_IS_RDS else "not_a_file_path__just_random_output_bis_to_ignore" 
    params:
        ref_name = REF_NAME,
        ref = True
    singularity: PYTHON_SINGULARITY
    script: os.path.join(SCRIPTS_DIR, "02bis_generate_h5ad_data_files.py")

rule generate_h5ad_query_files:
    input:
        gene_names = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "genes_names.csv"),
        counts = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "counts.mtx"),
        metadata = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "metadata.csv"),
        umap = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "umap.csv")
    output:
        h5ad_file = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "anndata.h5ad") if QUERY_IS_RDS else "not_a_file_path__just_random_output_bis_to_ignore" 
    params:
        query_name = QUERY_NAME,
        ref = False,
        umap_name = UMAP_NAME
    singularity: PYTHON_SINGULARITY
    script: os.path.join(SCRIPTS_DIR, "02bis_generate_h5ad_data_files.py")


rule annotate_py:
    input:
        ref_file = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "anndata.h5ad"),
        query_file = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "anndata.h5ad")
    output:
        #annotation = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "{query_name}_ref={ref_name}_method={ann_method}_predicted_celltype.csv"),
        scvi_embedding = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", f"{QUERY_NAME}_ref={REF_NAME}_method={{ann_method}}_scvi_embedding_annotated_adata.h5ad")
    params:
        query_name = QUERY_NAME,
        ann_method = "{ann_method}",
        ref_name = REF_NAME,
        extra_params = lambda wildcards: ANN_METHODS_PARAMS[wildcards.ann_method],
        script_method_name = lambda wildcards: wildcards.ann_method.split(":")[0],
        batch_key = BATCH_KEY,
        output_dir = OUTPUT_DIR
    wildcard_constraints:
        ann_method="|".join(ANN_METHODS),
    singularity: PYTHON_SINGULARITY
    script:
        os.path.join(SCRIPTS_DIR, "02_Ann_methods/{params.script_method_name}_annotation.py")


rule map_qc:
    input:
        scvi_embedding = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", f"{QUERY_NAME}_ref={REF_NAME}_method={{ann_method}}_scvi_embedding_annotated_adata.h5ad")
    output:
        metadata_csv = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", f"{QUERY_NAME}_ref={REF_NAME}_method={{ann_method}}_map_qc_score.csv")
    log:
        os.path.join(LOGS_PATH, f"mapqc_logfile_{REF_NAME}_{QUERY_NAME}_{{ann_method}}" + date_time_launch + ".log")
    params:
        batch_key = BATCH_KEY,
        ann_method = "{ann_method}",
        study_key = STUDY_KEY,
        ref_evaluation = REFERENCE_EVALUATION
    wildcard_constraints:
        ann_method="|".join(ANN_METHODS),
    singularity: PYTHON_SINGULARITY
    script:
        os.path.join(SCRIPTS_DIR, "02_Ann_methods/04_mapQC_evaluation.py")



rule annotate_metaclassifier:
    input:
        metadata_csv = expand(os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", f"{QUERY_NAME}_ref={REF_NAME}_method={{ann_method}}_map_qc_score.csv"), ann_method=ANN_METHODS)
    output:
        meta_metadata = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}__all_models_predictions_and_mapQC_scores.csv")
    wildcard_constraints:
    params:
        models = ANN_METHODS
    singularity: PYTHON_SINGULARITY
    script:
          os.path.join(SCRIPTS_DIR, "02_Ann_methods/99_MetaClassifier_annotation.py")


rule generate_report:
    input:
        ref = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "ref_data", "anndata.h5ad"),
        annotations_file = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}__all_models_predictions_and_mapQC_scores.csv"),
        query = os.path.join(OUTPUT_DIR, f"{QUERY_NAME}_ref={REF_NAME}_intermediate_outputs", "query_data", "anndata.h5ad")
    output:
        report = os.path.join(SCRIPTS_DIR, f"ann_tools_benchmark_on_ref={REF_NAME}_query={QUERY_NAME}" + date_time_launch + ".html")
    log: 
        os.path.join(LOGS_PATH, f"report_logfile_{REF_NAME}_{QUERY_NAME}" + date_time_launch + ".log")
    params:
        output_dir = OUTPUT_DIR,
        scripts_dir = SCRIPTS_DIR,
        reference_evaluation = REFERENCE_EVALUATION
    singularity: PYTHON_SINGULARITY
    script:
        os.path.join(SCRIPTS_DIR, "03_generate_report.py")

rule move_report:
    input:
        os.path.join(SCRIPTS_DIR, f"ann_tools_benchmark_on_ref={REF_NAME}_query={QUERY_NAME}" + date_time_launch + ".html")

    output:
        os.path.join(OUTPUT_DIR, f"ann_tools_benchmark_on_ref={REF_NAME}_query={QUERY_NAME}" + date_time_launch + ".html")
    shell:
        "mv {input} {output}"






