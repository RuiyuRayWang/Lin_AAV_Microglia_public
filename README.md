# Lin et al. AAV Microglia public

Characterization of microglia polarization states under targeted AAV transduction using bulk and single cell RNA sequencing.

Publication:
[Lin et al.]()

## Introduction

Characterization of microglia polarization states under targeted AAV transduction using bulk and single cell RNA sequencing.

The study of microglia biology and the development of microglia-based gene therapies are in urgent need of efficient and safe vehicles for microglia transgene delivery. To address this, we developed adeno-associated virus (AAV) variants that mediate efficient in vitro and in vivo microglia transduction via directed evolution of the AAV capsid protein. To assess the effect of AAV transduction on microglia, we carried out bulk RNAseq in primary microglia and found that microglia transduced by AAV remain close to homeostatic state. Furthermore, single-cell RNA sequencing showed that the AAV-MG variants mediate safe in vivo transgene delivery without inducing microglia immune activation. These AAV variants should facilitate the applications of various genetically-encoded sensors and effectors in studying microglia-related biology and therapeutic interventions.

## List of contents

The contents of this repository has been hosted on Github Page, which can be accessed [here](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/index.html).

Assessment of AAV transduction on *in vitro* cultured microglia via RNAseq:

  + [RNAseq vignette - Figure 1](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/Fig1_rnaseq.html)

Assessment of *in vivo* AAV transduction in microglia by scRNA-seq, run in order:

  * Preprocessing
    + [Reference Dataset (10X) - Supplementary Figure 2](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/FigS2_scRNAseq_reference_QC.html)
    + [Query Dataset (Smart-seq2) - Supplementary Figure 2](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/FigS2_scRNAseq_query_QC.html)
  * Primary Analyses
    + [Label Transfer, and Comparative Analysis between Query and Reference - Supplementary Figure 2](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/FigS2_scRNAseq_label_transfer.html)
    + [Main Figures 2](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/Fig2_main_figs.html)
  * Auxiliary Analysis
    + [Correlation Analysis of Signature Genes and Transgene (mScarlet) - Supplementary Figure 3](https://ruiyuraywang.github.io/Lin_AAV_Microglia_public/FigS3_scRNAseq_corr.html)

## Experiment setups

Overall design:

  * RNAseq: mRNA profiles of cultured microglia treated with blank, LPS, IL4 or/and AAV variants
  * scRNAseq: 
      + mRNA profiles of cold-dissociated microglia from mice treated with Saline, LPS, or Poly(i:c), using 10X Chromium (reference dataset); 
      + mRNA profiles of cold-dissociated microglia from mice treated with LPS or stereotaxically injected with AAV, using Smart-seq2 (query dataset).

Microglia were rapidly dissociated by Dounce homogenization under cold conditions and further enriched by Percoll density gradient centrifugation. e.g. the "classical" cold-mechanical dissociation protocol (Bennett et al., 2016). Steps involving FACS sorting were omitted to minimize *ex vivo* activation.

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
├── docs
├── notebooks
├── rscripts
└── vignettes                             (rmarkdowns)
```