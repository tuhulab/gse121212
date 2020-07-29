library(dplyr)

fastq <- readr::read_csv("ids/fastq.id", col_names = FALSE)
fastq_srr_is <- fastq$X1