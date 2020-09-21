library(dplyr)
library(ggplot2)
fastq <- readr::read_csv("ids/fastq.id", col_names = FALSE)
fastq_srr_id <- fastq$X1 %>% stringr::str_extract("SRR\\d{1,}") %>% unique()

all_srr_id <- readr::read_csv("ids/SRR_Acc_List.txt", col_names = FALSE) %>% pull(X1)
all_srr_id[!(all_srr_id %in% fastq_srr_id)] %>% data.frame() %>% readr::write_csv("ids/missing_srr_id.txt", col_names = FALSE)

all_srr_id %>% paste0(".sra") %>% data.frame() %>% readr::write_csv("ids/ids.txt", col_names = FALSE)

# Check MALAT1 gene
counttable <- readr::read_csv("data/count_cleaned.csv")
which(counttable$gene_id == "MALAT1")

tibble(sum_expression_level = counttable[, -1:-6] %>% rowSums(),
       gene_name = counttable$gene_id) %>% arrange(-sum_expression_level) %>% mutate(rank = 1:nrow(.)) %>% View()
