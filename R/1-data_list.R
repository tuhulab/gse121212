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
data_list$metadata <- tibble(BAM_ID = 
                               data_list$counttable %>% 
                               colnames() %>% 
                               stringr::str_match("SRR\\d{1,}\\.sra\\.bam"))


meta_data <- readr::read_csv("data/SraRunTable.txt") # %>% dplyr::select(Run, Skin_Type)
biosample <- readr::read_lines("data/biosample_result.txt")

group <- biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>% 
  stringr::str_match("AD|PSO|CTRL")

subject_id_df <- 
  biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>%
  stringr::str_match("(AD|PSO|CTRL)\\_\\d{3}")

subject_id <- subject_id_df[,1]
skin_type <- biosample[biosample %>% stringr::str_detect("\\d{1,}\\.")] %>% 
  stringr::str_remove("\\d{1,}\\.") %>%
  stringr::str_match("non-lesional|chronic_lesion|lesional|healthy")

meta_data_final <-
  meta_data %>% rename(SAMN_ID = BioSample,
                       SRR_ID = Run) %>%
  dplyr::select(SRR_ID, SAMN_ID) %>% 
  left_join(tibble(SAMN_ID = biosample %>% 
                     stringr::str_match("SAMN\\d{1,};") %>% 
                     purrr::discard(is.na) %>% stringr::str_remove(";"), 
                   group, subject_id, skin_type)) %>% 
  dplyr::select(-SAMN_ID)