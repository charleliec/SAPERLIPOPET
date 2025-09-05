library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(Matrix)
library(R.utils)
library(dplyr)
library(ggplot2)

library(MatrixExtra)
MatrixExtra::restore_old_matrix_behavior()

PATH_TO_OUTPUT_FOLDER = dirname(snakemake@output[["metadata"]])

seurat_obj = readRDS(file = snakemake@input[["rds_file"]])

suppressMessages({
    dir.create(file.path(PATH_TO_OUTPUT_FOLDER), recursive = TRUE)

    write.csv(seurat_obj@meta.data, snakemake@output[["metadata"]])

    counts = as.csc.matrix(t(seurat_obj$RNA$counts))
    writeMM(counts, snakemake@output[["counts"]])

    genes = Features(seurat_obj)
    write.csv(genes, snakemake@output[["gene_names"]])

    if (snakemake@params[["ref"]]=="False") {
        write.csv(seurat_obj@reductions[[snakemake@params[["umap_name"]]]]@cell.embeddings, snakemake@output[["umap"]]) 
    }
})
