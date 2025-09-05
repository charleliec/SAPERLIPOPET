print("launching R script Seurat")
#print(installed.packages())
library(Seurat)

dir.create(dirname(snakemake@output[["annotation"]]), recursive=TRUE)

PATH_TO_OUTPUT_FOLDER = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/labelTransferBenchmark/01_SeuratEvaluation/05_Output"
# INPUT_FOLDER = "00_DataPreparation"
# ANNOTATIONS_FOLDER = "01_Predictions"

# DATA_FILE_NAME = "PBMC_V2_VF1_AllGenes"
# DATASET_NAME = "Dataset3"
# PATIENTS = c("GSM5584156_1", "GSM5584156_2", "GSM5584156_3")

#source("/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/carer/labelTransferBenchmark/01_SeuratEvaluation/03_Script/01_SeuratEvaluation/00_generalDeps.R")


#Ref preprocessing

#Feature selection
# unique_patients <- unique(seuData.ref$orig.ident)
# gene_non_zero_in_patient <- matrix(0, nrow = nrow(seuData[["RNA"]]$counts), ncol = length(unique_patients),
#                     dimnames = list(rownames(seuData[["RNA"]]$counts), unique_patients))
# # Loop over patients to sum counts per gene
# for (patient in unique_patients) {
#   cells_in_pat <- colnames(seuData)[seuData$orig.ident == patient]
#   gene_non_zero_in_patient[, patient] <- Matrix::rowSums(seuData[["RNA"]]$counts[, cells_in_pat]) > 0
# }
# patient_count_per_gene = Matrix::rowSums(gene_non_zero_in_patient)
# genes_present = names(patient_count_per_gene[patient_count_per_gene > 6])
# seuData.ref = subset(seuData.ref, features=genes_present)

# Normalization + feature selection
seuData.ref = readRDS(snakemake@input[["ref"]])
seuData.ref[["RNA"]] = split(seuData.ref[["RNA"]], f=seuData.ref$orig.ident)

seuData.ref = NormalizeData(seuData.ref) #Multi batch norm
seuData.ref = FindVariableFeatures(seuData.ref)
seuData.ref = ScaleData(seuData.ref)
seuData.ref = RunPCA(seuData.ref)

# Integration of patients data
seuData.ref = IntegrateLayers(object=seuData.ref, method=CCAIntegration, orig.reduction="pca", new.reduction="integrated.cca")
seuData.ref = FindNeighbors(seuData.ref, dims = 1:30, reduction = "integrated.cca")
seuData.ref = FindClusters(seuData.ref, resolution = 0.3)
#seuData.ref = RunUMAP(seuData.ref, reduction = "integrated.cca", reduction.name = "integrated.umap", dims = 1:30)

# UMAP after integration
# labels <- unique(c(as.character(seuData.ref$ident))) #,as.character(seuData.ref$predicted.id)))
# palette <- setNames(Seurat:::DiscretePalette(length(labels), palette = "glasbey"), labels)
# (DimPlot(seuData.ref, reduction = "integrated.umap", group.by = c("ident")) + scale_color_manual(values = palette)) + 
#   (DimPlot(seuData.ref, reduction = "integrated.umap", group.by = c("orig.ident")) )


seuData.query = readRDS(snakemake@input[["query"]])

# Query preprocessing
#seuData.query[["RNA"]] = split(seuData.query[["RNA"]], f=seuData.query$orig.ident)

seuData.query = NormalizeData(seuData.query)
# Mapping to reference and annotation
seuData.anchors = FindTransferAnchors(reference = seuData.ref, query = seuData.query, dims = 1:30, reference.reduction = "integrated.cca")
predictions = TransferData(anchorset = seuData.anchors, refdata = seuData.ref$ident, dims = 1:30)
seuData.query = AddMetaData(seuData.query, metadata = predictions)
write.csv(seuData.query$predicted.id, snakemake@output[["annotation"]])


















