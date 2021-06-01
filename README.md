# gse121212
 repo for reproducing gse121212 study

# Links
The publication: [Atopic Dermatitis Is an IL-13–Dominant Disease with Greater Molecular Heterogeneity Compared to Psoriasis](https://www.jidonline.org/article/S0022-202X(19)30007-7/fulltext#appsec1)

Raw data repository: [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE121212&o=acc_s%3Aa)

# Data
| No. | Type     | Comparision                                             | Read R object                                                                                                                        |
| :---|:---------|:--------------------------------------------------------| :------------------------------------------------------------------------------------------------------------------------------------|
| 1   | DGE list |aLSvsCO (acute lesion vs healthy control)                | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/aLS_CO_DGE.rds"))        |
| 2   | DGE list |aLSvsNL (acute lesion vs non-lesion (subject matched)    | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/aLS_NL_DGE.rds"))        |
| 3   | DGE list |cLSvsCO (chronic lesion vs healthy control)              | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/cLS_CO_DGE.rds"))        |
| 4   | DGE list |cLSvsNL (chronic lesion vs non-lesion (subject matched)  | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/cLS_NL_DGE.rds"))        |
| 5   | DGE list |NLvsCO (non-lesion vs healthy control)                   | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/ab5db865c53d44ee381294bbaec34c31013e1f52/data/rds/NL_CO_DGE.rds"))         |
| 6   | data list|raw counts (filtering out lowly expressed genes)         | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/9d921a28370fe4c2222e02722cde59f61e8d4f51/data/rds/gse121212_list_raw.rds"))|
| 1*  | DGE list |aLSvsCO_all (acute lesion vs healthy control, all genes) | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/13d9d300bfd98adc70f85315eab7d7d619bad3bc/data/rds/aLS_CO_DGE_all.rds"))    |
| 2*  | DGE list |aLSvsNL_all (acute lesion vs non-lesion (subject matched, all genes)| readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/9e4484a64028b51d987d6157ba0112b48c0f7323/data/rds/aLS_NL_DGE_all.rds"))        |
| 3*  | DGE list |cLSvsCO_all (chronic lesion vs healthy control, all genes)   | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/13d9d300bfd98adc70f85315eab7d7d619bad3bc/data/rds/cLS_CO_DGE_all.rds"))    |
| 4*  | DGE list |cLSvsNL_all (chronic lesion vs non-lesion (subject matched, all genes)  | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/fa5ca19185fdb4c429a03e61c2b0f495c8e5fd23/data/rds/cLS_NL_DGE.rds"))        |
| 5*  | DGE list |NLvsCO (non-lesion vs healthy control, all genes)        | readRDS(url("https://rawcdn.githack.com/tuhulab/gse121212/13d9d300bfd98adc70f85315eab7d7d619bad3bc/data/rds/NL_CO_DGE_all.rds"))     |

1-5: AD Differentially expressed genes (|log2FC|>1, adj. p<0.05); aLS - acute lesion; CO - healthy control; cLS - chronic lesion; NL - non-lesional
6: raw counts, metadata, gene annotation
1*-5*: AD (all genes)

# Materials and methods
## Bioinformatics tools
[Generaly Applicable Gene-set Enrichment for Pathway Analysis (GAGE)](http://bioconductor.org/packages/release/bioc/html/gage.html)

[Enrichr](https://cran.r-project.org/web/packages/enrichR/index.html)
