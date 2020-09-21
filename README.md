# gse121212
 repo for reproducing gse121212 study

# Links
The publication: [Atopic Dermatitis Is an IL-13–Dominant Disease with Greater Molecular Heterogeneity Compared to Psoriasis](https://www.jidonline.org/article/S0022-202X(19)30007-7/fulltext#appsec1)

Raw data repository: [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE121212&o=acc_s%3Aa)

# Data
| No. | Comparision                                              | Read R object                                                                                                                 |
| :---|:---------------------------------------------------------| :-----------------------------------------------------------------------------------------------------------------------------|
| 1   | aLSvsCO (acute lesion vs healthy control)                | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/aLS_CO_DGE.rds")) |
| 2   | aLSvsNL (acute lesion vs non-lesion (subject matched)    | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/aLS_NL_DGE.rds")) |
| 3   | cLSvsCO (chronic lesion vs healthy control)              | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/cLS_CO_DGE.rds")) |
| 4   | cLSvvNL (chronic lesion vs non-lesion (subject matched)  | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/cLS_NL_DGE.rds")) |
| 5   | NLvsCO (non-lesion vs healthy control)                   | readRDS(url(https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/NL_CO_DGE.rds"))   |


# Materials and methods
## Bioinformatics tools
[Generaly Applicable Gene-set Enrichment for Pathway Analysis (GAGE)](http://bioconductor.org/packages/release/bioc/html/gage.html)

[Enrichr](https://cran.r-project.org/web/packages/enrichR/index.html)
