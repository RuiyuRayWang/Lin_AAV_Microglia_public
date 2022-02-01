#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(tidyverse)

querydata <- LoadH5Seurat('data/AAV_MG_scRNAseq/cells_query.h5Seurat')  # Preprocessed, QCed Smart-seq2 microglia (Query dataset)
refdata <- LoadH5Seurat(file = 'data/AAV_MG_scRNAseq/refdata_clean.h5Seurat')  # Preprocessed, QCed 10X (Reference dataset)

## Label transfer to identify cell states
## https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

dims_use <- 1:45

anchors <- FindTransferAnchors(
  reference = refdata,
  query = querydata,
  k.filter = 50,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = dims_use
)

refdata <- RunUMAP(refdata, dims = dims_use, reduction.name = "umap", seed.use = 42, return.model = TRUE)

query.transferred <- MapQuery(
  anchorset = anchors,
  query = querydata,
  reference = refdata,
  refdata = list(
    state = "state",
    treat = "treat"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

## Filter cells with low certainty
## Prediction score lower than 0.8 are likely incorrect (Stuart et al. 2018)
query.transferred.filtered <- subset(query.transferred, subset = predicted.treat.score > 0.8)

## Computing a 'de novo' UMAP visualization
refquery <- merge(refdata, query.transferred.filtered)
refquery[['pca']] <- merge(refdata[['pca']], query.transferred.filtered[['ref.pca']])
a = 0.9922; b = 1.112
n.epochs = 600
refquery <- RunUMAP(refquery, reduction = "pca", dims = 1:45, 
                    n.epochs = n.epochs,
                    a = a, b = b
                    )

# Update metadata
aav <- WhichCells(query.transferred, expression = treatment %in% c("AAV-MG1.2","AAV-MG1.1"))
lps <- WhichCells(query.transferred, expression = treatment == "LPS")
Idents(refquery) <- "treat"
refquery <- SetIdent(refquery, cells = aav, value = "AAV transduction")
refquery <- SetIdent(refquery, cells = lps, value = "LPS (Query)")
refquery[["treat"]] <- Idents(refquery)

Idents(refquery) <- "state"
hom_pred <- WhichCells(query.transferred, expression = predicted.state == "Homeostatic")
reac_pred <- WhichCells(query.transferred, expression = predicted.state == "Reactive")
refquery <- SetIdent(refquery, cells = hom_pred, value = "Homeostatic (Predicted)")
refquery <- SetIdent(refquery, cells = reac_pred, value = "Reactive (Predicted)")
refquery[["state"]] <- Idents(refquery)

transduced <- WhichCells(querydata, expression = mScarlet > 0)
non_transduced <- WhichCells(querydata, expression = mScarlet == 0)
q_lps <- WhichCells(query.transferred.filtered, expression = treatment == "LPS")
Idents(refquery) <- "treat"
refquery <- SetIdent(refquery, cells = transduced, value = "Transduced")
refquery <- SetIdent(refquery, cells = non_transduced, value = "Non-transduced")
refquery <- SetIdent(refquery, cells = q_lps, value = "LPS(Query)")
refquery[["ident"]] <- Idents(refquery)

Idents(query.transferred.filtered) <- "transduction"
query.transferred.filtered <- SetIdent(query.transferred.filtered, cells = q_lps, value = "LPS(Query)")
query.transferred.filtered[["ident"]] <- Idents(query.transferred.filtered)
query.transferred.filtered$ident <- as.character(query.transferred.filtered$ident)

SaveH5Seurat(refquery, filename = 'data/AAV_MG_scRNAseq/refquery.h5Seurat', overwrite = T)
SaveH5Seurat(query.transferred.filtered, filename = 'data/AAV_MG_scRNAseq/cells_query_transferred.h5Seurat', overwrite = T)
Convert('data/AAV_MG_scRNAseq/cells_query_transferred.h5Seurat', dest = "h5ad", overwrite = T)

## Copy the files to 'scenic_protocol/files' for SCENIC analyses
f_loom = c('data/AAV_MG_scRNAseq/cells_query_transferred.loom', 'data/scenic_protocol/files/cells_query_transferred.loom')
if(any(file.exists(f_loom))){
  ## Overwrite loom causes pyscenic protocol to fail, remove existing file
  file.remove(f_loom)
}
SaveLoom(query.transferred.filtered, filename = f_loom[1])
file.copy(from = f_loom[1],
          to = f_loom[2])

sessionInfo()
