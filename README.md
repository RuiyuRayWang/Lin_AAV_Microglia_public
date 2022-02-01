# Lin et al. AAV Microglia public
Code for analyzing scRNA-seq data (in review) in Lin et al.

## Introduction

The study of microglia biology and the development of microglia-based gene therapies are in urgent need of efficient and safe vehicles for microglia transgene delivery. To address this, we developed adeno-associated virus (AAV) variants that mediate efficient in vitro and in vivo microglia transduction via directed evolution of the AAV capsid protein. To assess the effect of AAV transduction on microglia, we carried out bulk RNAseq in primary microglia and found that microglia transduced by AAV remain close to homeostatic state. Furthermore, single-cell RNA sequencing showed that the AAV-MG variants mediate safe in vivo transgene delivery without inducing microglia immune activation. These AAV variants should facilitate the applications of various genetically-encoded sensors and effectors in studying microglia-related biology and therapeutic interventions.

## Experiment details

Overall design:

    - RNAseq: mRNA profiles of cultured microglia treated with blank, LPS, IL4 or/and AAV variants
    - scRNAseq: 
        * mRNA profiles of cold-dissociated microglia from mice treated with Saline, LPS, or Poly(i:c), using 10X Chromium (reference dataset); 
        * mRNA profiles of cold-dissociated microglia from mice treated with LPS or stereotaxically injected with AAV, using Smart-seq2 (query dataset).

To minimize "*ex vivo*" activation, microglia were rapidly dissociated by Dounce homogenization under cold conditions and further enriched by Percoll gradient density centrifugation. i.e. the "classical" cold-mechanical dissociation protocol (Bennett et al., 2016).

## Repository structure

```
├── .gitignore
├── LICENSE
├── README.md
├── data                   (not shown)
|   ├── AAV_MG_RNAseq      (data for RNAseq)
|   ├── AAV_MG_scRNAseq    (data for scRNA-seq)
|   └── GSE                (processed data deposited on GEO)
├── misc                   (miscellaneous. i.e. metadata)
|   ├── query_meta.csv
|   ├── ref_meta.csv.gz
|   └── RNAseq_sample_table.csv
├── notebooks
├── rscripts
└── vignettes              (rmarkdowns to generate plots in the original publication)
```