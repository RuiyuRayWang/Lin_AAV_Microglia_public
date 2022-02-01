#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)

query_mat <- read.table("data/GSE/AAV_MG_scRNAseq_SMARTSEQ2.tsv", header = T)
querydata <- CreateSeuratObject(query_mat, project = "MG_AAV_QUERY")
query_meta <- read.csv("misc/query_meta.csv", row.names = 1)
querydata <- AddMetaData(querydata, metadata = query_meta)
SaveH5Seurat(querydata, "data/AAV_MG_scRNAseq/cells_query.h5Seurat")

ref_mat <- read.table("data/GSE/AAV_MG_scRNAseq_10X.tsv", header = T)
refdata <- CreateSeuratObject(ref_mat, project = "MG_AAV_REF")
ref_meta <- read.csv(gzfile("misc/ref_meta.csv.gz"), row.names = 1)
refdata <- AddMetaData(refdata, metadata = ref_meta)
SaveH5Seurat(refdata, "data/AAV_MG_scRNAseq/refdata_clean.h5Seurat")

sessionInfo()
