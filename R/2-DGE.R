# DGE analysis (GSE121212)
library(dplyr)
library(DESeq2)
library(BiocParallel)

register(MulticoreParam(detectCores()-1))

## Load data
data_list <- readRDS("data/rds/gse121212_list_raw.rds")

## Manipulate data
index_ad_pso <- which(data_list$metadata$group %in% c("PSO", "CTRL"))

data_list_PSO_CO <- data_list
data_list_PSO_CO$metadata <- 
  data_list$metadata[index_ad_pso,] %>% 
  mutate(skin_type = skin_type %>% 
           stringr::str_replace("healthy", "CO")  %>% 
           stringr::str_replace("non-lesional", "NL") %>% 
           stringr::str_replace("lesional", "LS") %>% 
           factor(levels = c("CO", "NL", "LS"))) 
data_list_PSO_CO$counttable <- data_list_PSO_CO$counttable[,index_ad_pso]

## DGE analysis
DGE <- DESeqDataSetFromMatrix(countData = data_list_PSO_CO$counttable,
                              colData = data_list_PSO_CO$metadata,
                              design = ~ Sex + skin_type)
DGE_result <- DESeq(DGE)
resultsNames(DGE_result)

## DGE-Extract significant genes -----
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
                                                df_gene_annotation = data_list_PSO_CO$gene_annotation)
cLS_CO_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_result, 
                                                coef_name = "skin_type_cLS_vs_CO", 
                                                df_gene_annotation = data_list_PSO_CO$gene_annotation)
NL_CO_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_result, 
                                               coef_name = "skin_type_NL_vs_CO", 
                                               df_gene_annotation = data_list_PSO_CO$gene_annotation)
saveRDS(aLS_CO_significant, "data/rds/aLS_CO_DGE.rds")
saveRDS(cLS_CO_significant, "data/rds/cLS_CO_DGE.rds")
saveRDS(NL_CO_significant, "data/rds/NL_CO_DGE.rds")

## DGE-Extract all genes ----------------
extract_all_genes <- function(DESeq2_DGE_result = ..., 
                              coef_name = ...,
                              lfc_shrink_method = "apeglm",
                              df_gene_annotation = ...){
  lfcshrink_result <- lfcShrink(DESeq2_DGE_result, coef=coef_name, type=lfc_shrink_method)
  all_result <- bind_cols(df_gene_annotation, lfcshrink_result %>% as.data.frame())
}
PSO_LS_CO_all <- extract_all_genes(DESeq2_DGE_result = DGE_result, coef_name = "skin_type_LS_vs_CO", df_gene_annotation = data_list_PSO_CO$gene_annotation)
PSO_NL_CO_all <- extract_all_genes(DESeq2_DGE_result = DGE_result, coef_name = "skin_type_NL_vs_CO", df_gene_annotation = data_list_PSO_CO$gene_annotation)
saveRDS(PSO_LS_CO_all, "data/rds/PSO/PSO_LS_CO_DGE_all.rds")
saveRDS(PSO_NL_CO_all, "data/rds/PSO/PSO_NL_CO_DGE_all.rds")


## DGE analysis - PSO_LSvsNL_subject_matched
index_pso <- which(data_list$metadata$group == "PSO")
data_list_PSO <- data_list
data_list_PSO$metadata <- 
  data_list$metadata[index_pso,] %>% 
  mutate(skin_type = skin_type %>% 
           stringr::str_replace("healthy", "CO")  %>% 
           stringr::str_replace("non-lesional", "NL") %>% 
           stringr::str_replace("lesional", "LS") %>% 
           factor(levels = c("NL", "LS"))) 
data_list_PSO$counttable <- data_list_PSO$counttable[,index_pso]

DGE_sub_match <- DESeqDataSetFromMatrix(countData = data_list_PSO$counttable,
                                        colData = data_list_PSO$metadata,
                                        design = ~ subject_id + skin_type)
DGE_sub_match_result <- DESeq(DGE_sub_match)
resultsNames(DGE_sub_match_result)
LS_NL_significant <- extract_significant_genes(DESeq2_DGE_result = DGE_sub_match_result, 
                                                coef_name = "skin_type_LS_vs_NL", 
                                                df_gene_annotation = data_list_PSO$gene_annotation)

LS_NL_all <- extract_all_genes(DESeq2_DGE_result = DGE_sub_match_result, 
                               coef_name = "skin_type_LS_vs_NL", 
                               df_gene_annotation = data_list_PSO$gene_annotation)
saveRDS(LS_NL_all, "data/rds/PSO/PSO_LS_NL_DGE_all.rds")
