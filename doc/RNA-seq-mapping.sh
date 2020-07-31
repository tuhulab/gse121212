#!/bin/bash

# load module
echo "Step 0 - Load modules".
module load tools parallel/20200522 anaconda2/4.4.0 hisat2/2.2.0 subread/1.6.2

# define ids
ID_PATH=/home/projects/ku_00015/people/tuhu/gse121212/ids/ids.txt

# mapping
echo "Step 2 - Mapping"
mkdir -p /home/projects/ku_00015/data/gse121212_data/mapping/hisat2
mkdir -p /home/projects/ku_00015/data/gse121212_data/mapping/hisat2-report
cat ${ID_PATH} | parallel "hisat2 -x /home/projects/ku_00015/data/rna-seq/reference/genome-GRCh38/genome -1  /home/projects/ku_00015/data/gse121212_data/reads/{}_1.fastq -2 /home/projects/ku_00015/data/gse121212_data/reads/{}_2.fastq -S /home/projects/ku_00015/data/gse121212_data/mapping/hisat2/{}.bam --summary-file /home/projects/ku_00015/data/gse121212_data/mapping/hisat2-report/{}-report.summary"

