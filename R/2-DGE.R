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

extract_significant_genes <- function(DESeq2_DGE_result = ..., 
                                      coef_name = ...,
                                      lfc_shrink_method = "apeglm",
                                      df_gene_annotation = ...){
  lfcshrink_result <- lfcShrink(DESeq2_DGE_result, coef=coef_name, type=lfc_shrink_method)
  significant_result <- bind_cols(df_gene_annotation, lfcshrink_result %>% as.data.frame()) %>% 
    filter(padj < .05, abs(log2FoldChange) > 1)
}

aLS_CO_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_result, 
                                                coef_name = "skin_type_aLS_vs_CO", 
                                                df_gene_annotation = data_list_AD_CO$gene_annotation)
cLS_CO_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_result, 
                                                coef_name = "skin_type_cLS_vs_CO", 
                                                df_gene_annotation = data_list_AD_CO$gene_annotation)
NL_CO_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_result, 
                                               coef_name = "skin_type_NL_vs_CO", 
                                               df_gene_annotation = data_list_AD_CO$gene_annotation)
saveRDS(aLS_CO_significant, "data/rds/aLS_CO_DGE.rds")
saveRDS(cLS_CO_significant, "data/rds/cLS_CO_DGE.rds")
saveRDS(NL_CO_significant, "data/rds/NL_CO_DGE.rds")

## DGE analysis - aLSvsNL_aLSvscLS_subject_matched
index_ad <- which(data_list$metadata$group == "AD")
data_list_AD <- data_list
data_list_AD$metadata <- 
  data_list$metadata[index_ad,] %>% 
  mutate(skin_type = skin_type %>% 
           stringr::str_replace("healthy", "CO")  %>% 
           stringr::str_replace("non-lesional", "NL") %>% 
           stringr::str_replace("lesional", "aLS") %>% 
           stringr::str_replace("chronic_lesion", "cLS") %>% 
           factor(levels = c("NL", "cLS", "aLS"))) 
data_list_AD$counttable <- data_list_AD$counttable[,index_ad]

DGE_sub_match <- DESeqDataSetFromMatrix(countData = data_list_AD$counttable,
                                        colData = data_list_AD$metadata,
                                        design = ~ subject_id + skin_type)
DGE_sub_match_result <- DESeq(DGE_sub_match)
resultsNames(DGE_sub_match_result)
aLS_NL_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_sub_match_result, 
                                                coef_name = "skin_type_aLS_vs_NL", 
                                                df_gene_annotation = data_list_AD$gene_annotation)
cLS_NL_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_sub_match_result, 
                                                coef_name = "skin_type_cLS_vs_NL", 
                                                df_gene_annotation = data_list_AD$gene_annotation)