# Lin et al. AAV Microglia public

Characterization of microglia polarization states under targeted AAV transduction using bulk and single cell RNA sequencing.  
[**View the hosted documents on Github Pages**](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/index.html).

For full details on the evolution and screening of the AAV capsid, check out the published version of our paper:

[Lin, R. et al. Directed evolution of adeno-associated virus for efficient gene delivery to microglia. Nat Methods 1–10 (2022) doi:10.1038/s41592-022-01547-7.](https://www.nature.com/articles/s41592-022-01547-7)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Introduction

The study of microglia biology and the development of microglia-based gene therapies are in urgent need of efficient and safe vehicles for microglia transgene delivery. To address this, we developed adeno-associated virus (AAV) variants that mediate efficient in vitro and in vivo microglia transduction via directed evolution of the AAV capsid protein. To assess the effect of AAV transduction on microglia, we carried out bulk RNAseq in primary microglia and found that microglia transduced by AAV remain close to homeostatic state. Furthermore, single-cell RNA sequencing showed that the AAV-MG variants mediate safe in vivo transgene delivery without inducing microglia immune activation. These AAV variants should facilitate the applications of various genetically-encoded sensors and effectors in studying microglia-related biology and therapeutic interventions.

## List of contents

Assessment of *in vivo* AAV transduction in microglia by scRNA-seq:

  * Preprocessing
    + [Query Dataset (Smart-seq2) - Extended Figure 7](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/ExtFig7_scRNAseq_query_QC.html)
    + [Reference Dataset (10X) - Extended Figure 7](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/ExtFig7_scRNAseq_reference_QC.html)
  * Primary Analyses
    + [Label Transfer, and Comparative Analysis between Query and Reference - Extended Figure 7](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/ExtFig7_scRNAseq_label_transfer.html)
    + [Descriptive Analysis, Label Transfer, and Odds Ratio - Figure 3](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/Fig3_main_figs.html)
  * Auxiliary Analysis
    + [Correlation Analysis of Signature Genes and Transgene (mScarlet) - Extended Figure 8](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/ExtFig8_scRNAseq_corr.html)

Assessment of AAV transduction on *in vitro* cultured microglia via RNAseq:

  - [RNAseq vignette - Extended Figure 1](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/ExtFig1_rnaseq.html)

## Experiment design

Overall design:

  * RNAseq: mRNA profiles of cultured microglia treated with blank, LPS, IL4 or/and AAV variants
  * scRNAseq: 
      + mRNA profiles of cold-dissociated microglia from mice treated with Saline, LPS, or Poly(i:c), using 10X Chromium (reference dataset); 
      + mRNA profiles of cold-dissociated microglia from mice treated with LPS or stereotaxically injected with AAV, using Smart-seq2 (query dataset).

Microglia were rapidly dissociated by Dounce homogenization under cold conditions and further enriched by Percoll density gradient centrifugation. e.g. the "classical" cold-mechanical dissociation protocol (Bennett et al., 2016). Steps involving FACS sorting were omitted to minimize *ex vivo* activation.

## Data and code availability

The sequencing data used in this study can be accessed via GEO under the accession number [GSE197743](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197743).  
You may also follow the above vignettes and use `GEOquery` to download the data programmatically.

Codes reproducing the analysis and figures in the original publication are written with Rmarkdowns, rendered as html-pages and hosted on Github pages.  
If you'd like to reproduce the vignettes on your own, simply clone the repository:
```
git clone https://github.com/RuiyuRayWang/Lin_AAV_Microglia_public.git
```

## Dependencies

Some packages are required in order to build and run the vignettes. You may install them as follows:

```{r eval=FALSE}
install.packages("tidyverse")
install.packages('Seurat')
install.packages("ggpubr")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("cvms")
install.packages("epitools")
install.packages("scales")
```

```{r eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("scater")
BiocManager::install("SummarizedExperiment")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DESeq2")
```

```{r eval=FALSE}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
```

## Repository structure

```
├── .gitignore
├── LICENSE
├── README.md
├── data                                  (omitted in repo to save space)
|   ├── AAV_MG_scRNAseq                   (processed data for scRNA-seq)
|   └── GSE197743                         (generated by GEOquery, downloads data deposited on GEO)
├── misc                                  (miscellaneous. i.e. metadata)
|   ├── mg_manual.rda
|   ├── awk_extract_gtf.sh                (bash script parsing gene name conversion table)
|   ├── MM.GRCm38.102.annotation.tab
|   └── RNAseq_sample_table.csv
├── docs                                  (knitted html documents)
└── vignettes                             (rmarkdowns)
```
