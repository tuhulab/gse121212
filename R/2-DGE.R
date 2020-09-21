# DGE analysis (GSE121212)
library(dplyr)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(detectCores()-1))

## Load data
data_list <- readRDS("data/rds/gse121212_list_raw.rds")

## Manipulate data
index_ad_co <- which(data_list$metadata$group %in% c("AD", "CTRL"))
data_list_AD_CO <- data_list
data_list_AD_CO$metadata <- 
  data_list$metadata[index_ad_co,] %>% 
  mutate(skin_type = skin_type %>% 
           stringr::str_replace("healthy", "CO")  %>% 
           stringr::str_replace("non-lesional", "NL") %>% 
           stringr::str_replace("lesional", "aLS") %>% 
           stringr::str_replace("chronic_lesion", "cLS") %>% 
           factor(levels = c("CO", "NL", "cLS", "aLS"))) 
data_list_AD_CO$counttable <- data_list_AD_CO$counttable[,index_ad_co]
## DGE analysis
DGE <- DESeqDataSetFromMatrix(countData = data_list_AD_CO$counttable,
                              colData = data_list_AD_CO$metadata,
                              design = ~ Sex + skin_type)
DGE_result <- DESeq(DGE)
resultsNames(DGE_result)
