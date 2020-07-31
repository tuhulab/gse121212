#!/bin/bash

# feature counting
module load tools subread/1.6.2
mkdir -p /home/projects/ku_00015/people/tuhu/gse121212/counts
featureCounts -T 32 -t exon -g gene_id -a /home/projects/ku_00015/data/rna-seq/reference/gtf/gencode.v34.annotation.gtf -o /home/projects/ku_00015/people/tuhu/gse121212/counts/counts.txt /home/projects/ku_00015/data/gse121212_data/mapping/hisat2/*.bam

