# DGE analysis (GSE121212)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(stringr)
library(SummarizedExperiment)
library(tidySummarizedExperiment)
library(tidybulk)
library(purrr)
register(MulticoreParam(future::availableCores()-4))
## Load data
data_list <- 
  readRDS("data/rds/gse121212_list_raw.rds")
## SummarizedExperiement
col_data <- 
  data_list$metadata %>% 
  mutate(skin_type = case_when(skin_type == "lesional" ~ "LS",
                               skin_type == "chronic_lesion" ~ "LS",
                               skin_type == "healthy" ~ "CO",
                               skin_type == "non-lesional" ~ "NL"),
         skin_type = paste(group, skin_type, sep = "-") %>% str_remove("CTRL-")) %>% 
  mutate(skin_type = forcats::fct_relevel(skin_type, c("CO", "PSO-NL", "PSO-LS", "AD-NL", "AD-LS")))

se <- SummarizedExperiment(assays = list(counts = data_list$counttable %>% as.matrix),
                           colData = col_data)
names(se) <- data_list$gene_annotation$gene_id
se <- se[rownames(se) %>% stringr::str_detect("\\d{1,}P$|\\d{1,}P\\d{1,}$|\\.|-AS\\d{1}|-DT", negate = TRUE),] # gene cleaner
se <- se %>% keep_abundant(factor_of_interest = skin_type)
DGE_design <- 
  tibble(
    C1 = c("AD-LS", "AD-LS", "PSO-LS", "PSO-LS", "AD-LS", "AD-NL", "AD-NL", "PSO-NL"),
    C2 = c("AD-NL", "CO", "PSO-NL", "CO", "PSO-LS", "PSO-NL", "CO", "CO")
  ) %>% 
  mutate(
    C1_t = str_extract(C1, "AD|PSO|CO"),
    C2_t = str_extract(C2, "AD|PSO|CO"),
    contrast = paste0(C1, "vs", C2),
    subject_matched = (C1_t == C2_t),
    se_exp = map2(C1, C2, function(c1, c2){se[,se$skin_type %in% c(c1, c2)]}),
    design_f = ifelse(subject_matched, 
                      "~ subject_id + skin_type",
                      "~ Sex + skin_type"),
    deseq_dataset = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula())),
    deseq_res = map(deseq_dataset, ~ DESeq(.x, parallel = TRUE)),
    res_name = map_chr(deseq_res, ~ grep("skin_type", resultsNames(.x), value = TRUE)),
    res_shrink = map2(deseq_res, res_name, ~ lfcShrink(.x, coef = .y, type = "apeglm", parallel = TRUE)),
    res_shrink = map(res_shrink, ~ as_tibble(.x, rownames = "gene_name")))

saveRDS(DGE_design %>% select(contrast, res_shrink), "data/DGE_full_model.rds")
