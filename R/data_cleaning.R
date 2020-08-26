# Tidying up GSE121212 countable -----------------------
library(dplyr)

counttable <- readr::read_tsv("counts/gene_name/counts.txt", skip=1)

counttable_annotate <- 
  counttable %>% 
  rename(gene_id = Geneid) %>% 
  mutate(Chr = Chr %>% str_extract("chr\\w{1,}(?=;)|chr\\w{1,}$")) %>% 
  mutate(Strand = Strand %>% str_extract("\\+|\\-")) %>% 
  mutate(Start = Start %>% str_extract("\\d{1,}(?=;)|\\d{1,}$"), 
         End = End %>% str_extract("\\d{1,}(?=;)|\\d{1,}$"))

counttable_annotate %>% write_csv("data/count_cleaned.csv")