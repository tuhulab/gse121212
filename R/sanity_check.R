library(dplyr)

fastq <- readr::read_csv("ids/fastq.id", col_names = FALSE)
fastq_srr_id <- fastq$X1 %>% stringr::str_extract("SRR\\d{1,}") %>% unique()

all_srr_id <- readr::read_csv("ids/SRR_Acc_List.txt", col_names = FALSE) %>% pull(X1)
all_srr_id[!(all_srr_id %in% fastq_srr_id)] %>% data.frame() %>% readr::write_csv("ids/missing_srr_id.txt", col_names = FALSE)
