# load libraries
library(dplyr)

# generate data list
data_list <- list()
count_data <- readr::read_csv("data/count_cleaned.csv")
counts <- count_data[, -1:-6]

# filter out lowly-expressed genes
mean_expression <- counts %>% as.matrix() %>% matrixStats::rowMeans2()
gene_keep <- mean_expression >= 1 
sum(gene_keep) # no of genes kept

data_list$gene_annotation <- count_data[gene_keep, 1:6]
data_list$counttable <- count_data[gene_keep, -1:-6]

# construct metadata
sra_run_table <- readr::read_csv("data/SraRunTable.txt") # %>% dplyr::select(Run, Skin_Type)
biosample <- readr::read_lines("data/biosample_result.txt")
group <- biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>% 
  stringr::str_match("AD|PSO|CTRL") %>% as.character()
subject_id_df <- 
  biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>%
  stringr::str_match("(AD|PSO|CTRL)\\_\\d{3}")
subject_id <- subject_id_df[,1]
skin_type <- biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>%
  stringr::str_match("non-lesional|chronic_lesion|lesional|healthy") %>% as.character()
sra_run_table_final <-
  sra_run_table %>% rename(SAMN_ID = BioSample,
                       SRR_ID = Run) %>%
  dplyr::select(SRR_ID, SAMN_ID) %>% 
  left_join(tibble(SAMN_ID = biosample %>% 
                     stringr::str_match("SAMN\\d{1,};") %>% 
                     purrr::discard(is.na) %>% stringr::str_remove(";"), 
                   group, subject_id, skin_type)) %>% 
  dplyr::select(-SAMN_ID)

appendix_table2 <- readr::read_csv("data/metadata_Table2.csv") %>% 
  rename(FLG_mutation = `FLG Mutation`, 
         severity = `Severity\n(ScorAD/PASI)`,
         biopsy_site_non_lesional = `biopsy site non-\nlesional`,
         biopsy_site_lesional = `biopsy site lesional`,
         palmer_hyperlinearity = `Palmar hyperlinearity`,
         keratosis_pilaris = `Keratosis\npilaris`) %>% 
  mutate(FLG_mutation = FLG_mutation %>% stringr::str_replace_all(" ","_"),
         biopsy_site_lesional = biopsy_site_lesional %>% stringr::str_replace_all(" ","_"),
         biopsy_site_non_lesional = biopsy_site_non_lesional %>% stringr::str_replace_all(" ","_"))

data_list$metadata <- tibble(SRR_ID = 
                               data_list$counttable %>% 
                               colnames() %>% 
                               stringr::str_match("SRR\\d{1,}") %>% 
                               as.character()) %>% left_join(sra_run_table_final) %>% 
  left_join(appendix_table2, by = c("subject_id" = "ID.pheno"))


data_list %>% saveRDS("data/rds/gse121212_list_raw.rds")



### GSE121212
data_list <- readr::read_rds("data/rds/gse121212_list_raw.rds")

library(SummarizedExperiment)
library(tidybulk)
library(tidySummarizedExperiment)
library(stringr)

se <- SummarizedExperiment(assays = list(counts = data_list$counttable),
                           colData = data_list$metadata)
names(se) <- data_list$gene_annotation$gene_id
g_non_pseudo <- 
  rownames(se) %>% 
  gprofiler2::gconvert() %>% 
  filter(!description %>% str_detect("pseudogene")) %>% 
  pull(input) %>% unique()
se_filterG <- se[g_non_pseudo,] %>% scale_abundance()
saveRDS(se_filterG, "data/rds/se.rds")
