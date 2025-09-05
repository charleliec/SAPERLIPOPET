<div style="display:flex;align-items:center;">
  <img src="Saperlipopet_logo.png" alt="Logo" width="500">
  <h3 style="margin: 10px 0px 0px 50px;">Subtype Annotation with Prior Evaluation of the Reference LImits and Potential Outlier Population ExTraction</h3> 
</div>
<br>


Pipeline originally developped to annotate the subtypes of NK cells. It is based on scANVI for the annotation and uses MapQC to detect new cell types that are not present in the original reference.

- [1. Usage of the pipeline](#1-usage-of-the-pipeline)
  - [1.1. Launching the pipeline](#11-launching-the-pipeline)
  - [1.2. Inputs description](#12-inputs-description)
    - [1.2.1. Configuration file](#121-configuration-file)
    - [1.2.2. Data files](#122-data-files)
  - [1.3. Potential errors or problems](#13-potential-errors-or-problems)
  - [1.4. Outputs](#14-outputs)
- [2. Advices on best use and limits](#2-advices-on-best-use-and-limits)
- [3. Results interpretation](#3-results-interpretation)
- [4. MapQC usage](#4-mapqc-usage)
- [5. Notes for future developpers](#5-notes-for-future-developpers)


## 1. Usage of the pipeline

### 1.1. Launching the pipeline

If you copied the files of the pipeline to another folder, start by updating the `SCRIPTS_DIR` variable at the top of the `Snakefile` : if the path to the `Snakefile` is `/mnt/some/path/Snakefile` it will be `/mnt/some/path`. Finally make sure that the Singularity images exist and that their paths at the top of the `Snakefile` are correct.
To run the SAPERLIPOPExt pipeline, simply create virtual environment with :

```console
virtualenv -p /usr/bin/python3 ~/.snakemake
source ~/.snakemake/bin/activate
pip install snakemake
pip install pulp==2.7.0
```

Then run the following :

```console
snakemake --snakefile Snakefile --verbose --use-singularity --singularity-args="-B /mnt/DOSI:/mnt/DOSI -B $(pwd)/.cache:/home/jovyan/.cache --nv" --cores 5
```

Depending on your GPU capacity, you can adjust the number of cores. Usually, to run 4 models, it should take around 2 to 3h. Make sure that it is running on the GPU, otherwise it could take much longer (at least 10 times longer).

To create the containers you can first build the two Docker containers which Dockerfiles are in `SAPERLIPOPET/01_Containers/` and then create the Singularity image from the Docker images using :

```console
docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp/test:/output \
--privileged -t --rm \
quay.io/singularity/docker2singularity \
myDockerImageID
```

More detailed informations are present in the `readme.txt` of each container.

### 1.2. Inputs description

**Overview of the inputs :**

* Configuration file : `01_CONFIG.yaml`. Contains the paths to the reference, the query, etc. as well as the models parameters.
* A `.rds` or `.h5ad` file containing respectively a *Seurat* object or an *Anndata* object for the reference.
* A `.rds` or `.h5ad` file containing respectively a *Seurat* object or an *Anndata* object for the query.



#### 1.2.1. Configuration file

The `01_CONFIG.yaml` file must have the following fields :

* **`reference_evaluation`** : whether the file to annotate, the query, is labeled. If it is, the labels will be used to compare with the annotation of the pipeline.
* **`reference_path`** : path to the `.rds` or `.h5ad` file containing respectively a *Seurat* object or an *Anndata* object for the reference.
* **`query_path`** : path to the `.rds` or `.h5ad` file containing respectively a *Seurat* object or an *Anndata* object for the reference.
* **`output_dir`** : path to the directory where the output files will be stored
* **`batch_key`** : name of the column in the metadata containing the samples/patient information. It is what will be used by scANVI to integrate. The pipeline runs on **raw counts**, so any prior integration will be discarded. Normally, there are several samples per dataset / study.
* **`study_key`** : name of the column in metadata containing the information about the dataset / study the cell belongs to. It should have a unique value in the query.
=> Both the `study_key` and the `batch_key` have to be the **same in the query and in the reference**.
* **`umap_name`** : name of the umap reduction in the `seurat_obj.rds`. It will be used in the following way : `seurat_obj@reductions[[umap_name]]@cell.embeddings`. If providing directly the `counts, metadata,...` set to any value (e.g. `False`).
* **`models_params`** : sets the number of models as well as their training parameters :
    * **`ref_unlabeled_num_epochs`** : list with number of training epochs of the scVI model on the reference (thus without using the labels). This parameter doesn't seem to have much impact, between 20 and 200 seems like a good number. Please note that increasing the number of epochs will increase the computing time.
    * **`ref_labeled_num_epochs`** : list with number of training epochs of the scANVI model on the reference (using the labels). Usually numbers between 10 and 40 give good results
    * **`query_num_epochs`** : list with number of training epochs of the scANVI model on both the reference and the query. Usually numbers between 50 and 300 give good results
    * **`prod`** : if set to True, the cartesian product of the `..._num_epochs` lists is used to create the models. If False, they all must have the same length and zip() will be used.
=> The number of models trained is determined by the values in the lists as well as the value of `prod`. When `prod` is `True`, there will be $\texttt{len(ref\_unlabeled\_num\_epochs)} \times \texttt{len(ref\_labeled\_num\_epochs)} \times \texttt{len(query\_num\_epochs)}$ different models. Otherwise, there will only be $\texttt{len(ref\_unlabeled\_num\_epochs)} = \texttt{len(ref\_labeled\_num\_epochs)} = \texttt{len(query\_num\_epochs)}$ models.


#### 1.2.2. Data files

There are two ways to give the reference and query data to the pipeline : either through a *Seurat* object using a `.rds` file or through an *AnnData* using a `.h5ad` file. When naming those file, respect the exact extension name : don't use uppercase or alternative notations.

If you provide a `.rds` file it must follow those specifications :

* You must have joined all the **raw counts** layers together in the `RNA$counts` slot.
* The `meta.data` must contain a column named `ident` containing the true cell identity in the reference (and when evaluating the reference in the query as well).
* The `meta.data` must contain the columns whose names are defined in `batch_key` and `study_key`. The name in the query should be **identical to that in the reference**.
* The UMAP reduction defined in the config file must exist at the following location : `seurat_obj@reductions[[`*`umap_name`*`]]@cell.embeddings`.
* The gene names nomenclature has to match between the query and the reference. It is however not required that the same genes are present in the query and in the reference.

If you provide a `.h5ad` file, it must meet the following requirements :

* The **raw counts** must be located in `adata.X`. No other layer is needed.
* The `obs` DataFrame must contain a column named `ident` containing the true cell identity in the reference (and when evaluating the reference in the query as well).
* The `obs` must contain the columns whose names are defined in `batch_key` and `study_key`. The name in the query should be **identical to that in the reference**.
* The UMAP reduction defined in the config file must exist at the following location : `adata.obsm["umap"]`. For now, using `.h5ad` as input, it is not possible to select the name of the UMAP. Set the value of `umap_name` to some default value like `"random_text"`.
* The gene names nomenclature (set in `adata.var_names`) has to match between the query and the reference. It is however not required that the same genes are present in the query and in the reference.


### 1.3. Potential errors or problems

* If scANVI takes an unreasonable amount of time to be imported ( >10 min), restart the computing cluster, it should solve the issue.
* If after running the pipeline, you realize that MapQC is not very useful because an important number of cells is grey (not sampled), consider increasing the `n_nhoods` parameter in  `02_Ann_methods/04_mapQC_evaluation.py`.

### 1.4. Outputs
The annotations will be saved in `..._all_models_predictions_and_mapQC_scores.csv` where `...` corresponds to the name of the reference and that of the query.  
A summary of the results, including UMAPs of the prediction, of the MapQC scores and so forth will also be created under `ann_tools_benchmark_on_ref=`*`{name of reference}`*`_query=`*`{name of the query}`*`_`*`{date and time}`*`.html`.
Finally, a folder named *`{QUERY_NAME}`*`_ref=`*`{REF_NAME}`*`_intermediate_outputs` containing the different files generated during the pipeline.


## 2. Advices on best use and limits

The transfer of the labels between the reference and the query, depends mainly on the qulity of the reference.
It appears that the pipeline is not capable of transfering labels to new organs, chemistry or sequencing technique. For instance, it won't be possible to annotate lung cells with only PBMC cells in the reference, or to annotate cells sequenced with V3 chemistry if in the reference all cells were sequenced with V2 or worse to annotate cells with singleNucleus data with a reference only made of singleCell data.

However, if your reference is made of PBMC cells sequenced with V2 chemistry as well as lung cells sequenced with V3 chemistry, you *should* be able to annotate lung cells with V2 chemistry or PBMC cells with V3 chemistry.
It is also advised when annotating large datasets that you manually annotate at least one sample so as to evaluate the quality of your reference. This will give you an estimate of the quality of the annotation for the other samples. If the results are not satisfactory, it is advised to annotate manually a second sample, to add it to the reference and to evaluate once again the quality of the reference. Usually, if the original reference is comprised of enough cells, around 1000 cells in the new context should suffice to achieve good annotation performance.

Sometimes, it appears that adding unlabeled cells to the reference can increase performance. This might be due to the fact that it densifies the number of points in the latent space and allows for better clustering. However, it has proven not to be always the case, and some small performance degradation have been observed.

## 3. Results interpretation

In the HTML report, there are different indicators of the quality of the annotation as well as a detection of outlier population using MapQC.

At the top, in the **Overview** tab, there is two figure showing how well the different models trained on the data agree with each other. If all the models agree with each other (score close to $1$), it is usually a sign of a reliable annotation. However, although unlikely, it could simply be that they are all biased in the same way and all make the same mistake.

Underneath, you will find the aggregated results of MapQC. MapQC is a tool capable, given a query, a reference and a mapping between the two, to evaluate which cells in the query are outliers (cell population not present in the reference). Here the mapping used is the latent space of the *scANVI* model, therefore, MapQC is run as many times as there are models.  For a cell, a MapQC score below $2$ is a sign that the cell type belongs to those present in the reference. Conversely, a score above $2$ indicates either a poor integration (therefor a poor model) or a new cell type that was not present in the reference. Cells can also be marked as *Not sampled* or *Filtered out* this will be explained in more details in the section dedicated to MapQC.

On the **Models** tab, you will find the annotation given by each model as well their respective MapQC scores.

In the `..._all_models_predictions_and_mapQC_scores.csv` file, you will find the following columns :

* One column `predictions_...` per model, giving the prediction of this model
* Similarly one column starting with `mapqc_score_...` or `binary_mapqc_score_...` per model. The binary MapQC score is useful in UMAPs to distinguish between outlier cells and non outliers as well as visualizing which cells were not sampled or filtered out.
* `final_prediction` giving the final prediction of the pipeline
* `pred_label_proportion` giving the proportion of models that "voted" for the cell type having receiving the most votes
* `proportion_difference` giving the difference in vote proportion between the first most voted type and the second most voted type (if you have 10 models, 5 that voted for cell type A, 2 that voted for cell type B and the rest all voted for other cell types, you would have a `proportion_difference` of $5-2 = 3$) 
* `mapqc_median`, `mapqc_q3`, `mapqc_max` aggregating the MapQC score
* `binary_mapqc_median`,	`binary_mapqc_q3`,	`binary_mapqc_max` are the binarized version of the above columns (containing `Nan` if not a value)


## 4. MapQC usage

If you see too many cells marked as *Not sampled*, increase `n_nhoods`. This raises the chance that each query cell is included in at least one neighborhood. You can start from $\texttt{(number\_of\_query\_cells ÷ k\_min) × 5}$ and scale upward until coverage is sufficient. The only drawback of setting this parameter too high is a longer computing time. However, the calculation of the MapQC score is usually orders of magnitude smaller than that of scANVI so it is better to set it too high than too low.

If many cells are *Filtered out*, you can try raising `k_max` to allow neighborhoods to expand further, but keep it below $\texttt{~10× k\_min}$ to avoid mixing in a single neighbourhood different cell types.

To get more hindsight on the reasons why cells were not included in the calculation of the MapQC score, you can take a look at the `logs/mapqc_logfile_...`. It mainly contains one table giving the reason why each of the neighbourhood was filtered out or not, another with the sample distances to reference and the last one being the `.obs` DataFrame of the AnnData after training of the model and thus with the embedding of the cells in the latent space.

It is also advised to read the [MapQC documentation](https://mapqc.readthedocs.io/en/latest/notebooks/mapqc_detailed.html) which is very well done.



## 5. Notes for future developpers

During development, I had quite a lot of issue using scanpy in a Singularity container. Especially, scanpy would start importing but never finish, running indefinitly. To solve this issue, I had to create cache directories in the `/tmp` of the Singularity container. So when importing scanpy in a file put the following code before :

```python
# Needed for the import of scanpy
os.makedirs("/tmp/.cache", exist_ok=True)
os.environ["NUMBA_CACHE_DIR"] = "/tmp/cache/numba_cache"
os.environ['MPLCONFIGDIR'] = '/tmp/cache/mplconfig_cache'
os.environ["PYTORCH_KERNEL_CACHE_PATH"] = "/tmp/cache/torch_cache"
os.environ["FONTCONFIG_FILE"] = "/tmp/cache/fontconfig_cache"
os.environ["XDG_CACHE_HOME"] = "/tmp/cache"
for d in ["XDG_CACHE_HOME", "NUMBA_CACHE_DIR", "MPLCONFIGDIR", "PYTORCH_KERNEL_CACHE_PATH", "FONTCONFIG_FILE"]:
    os.makedirs(os.environ[d], exist_ok=True)
```

I faced similar issues when running Quarto in a Singularity container, which explains the `03_generate_report.py` file.

Further more to start with the same Singularity environment between runs, I had to add the the `--contain` option to the command used to run the pipeline. This option creates new directory in the Singularity environment for the `/tmp` and `/home` instead of using those of the host. However, since Snakemake adds a `--home` flag, the `/home` is actually bound to the `SCRIPTS_DIR`. This was necessary for the import of *Scanpy* due to the creation of the cache directories in the `/tmp`.











































