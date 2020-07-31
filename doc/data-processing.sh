#!/bin/bash

# download data
module load sratoolkit/2.9.6-1
prefetch --option-file /home/people/tuhu/Downloads/SRR_Acc_List.txt

# fasterq-dump data
cd /home/projects/ku_00015/data/sratool_data/sra
module load parallel/20200522 sratoolkit/2.9.6-1
cat /home/people/tuhu/Downloads/SRR_Acc_List.txt | parallel "fasterq-dump /home/projects/ku_00015/data/sratool_data/sra/{}.sra"

# download missing data
module load sratoolkit/2.10.7 parallel/20200522
prefetch --option-file /home/projects/ku_00015/people/tuhu/gse121212/ids/missing_srr_id.txt
cd /home/projects/ku_00015/data/sratool_data/sra
cat /home/projects/ku_00015/people/tuhu/gse121212/ids/missing_srr_id.txt | parallel "fasterq-dump /home/projects/ku_00015/data/sratool_data/sra/{}.sra"

# load module
echo "Step 0 - Load modules".
module load tools fastp/0.19.4 parallel/20200522 anaconda2/4.4.0 hisat2/2.2.0 subread/1.6.2

# initiate - ids
ID_PATH=/home/projects/ku_00015/people/tuhu/gse121212/ids/ids.txt

# QC(1) - fastqc
echo "Step 0 - FastQC"
module load perl/5.30.2 jdk/14 fastqc/0.11.8 parallel/20200522
mkdir /home/projects/ku_00015/data/gse121212_data/qc
cat ${ID_PATH} |  parallel "fastqc /home/projects/ku_00015/data/gse121212_data/reads/{}_1.fastq /home/projects/ku_00015/data/gse121212_data/reads/{}_2.fastq -o /home/projects/ku_00015/data/gse121212_data/qc/"

# QC(2) - multiqc
module load anaconda3/4.4.0
cd /home/projects/ku_00015/data/gse121212_data/
multiqc qc/*_fastqc.zip -o qc/

# trimming
# echo "Step 1 - Trimming"
# mkdir -p /home/projects/ku_00015/data/gse121212_data/trimmed/fastp
# mkdir -p /home/projects/ku_00015/data/gse121212_data/trimmed/fastp-report
# cat ${ID_PATH} |  parallel "fastp -i /home/projects/ku_00015/data/gse121212_data/reads/{}_1.fastq -I /home/projects/ku_00015/data/gse121212_data/reads/{}_2.fastq -o /home/projects/ku_00015/data/gse121212_data/trimmed/fastp/{}_1.fastq.gz -O /home/projects/ku_00015/data/gse121212_data/trimmed/fastp/{}_2.fastq.gz -w 1 -h /home/projects/ku_00015/data/gse121212_data/trimmed/fastp-report/{}-report.html"

# mapping
echo "Step 2 - Mapping"
mkdir -p /home/projects/ku_00015/data/gse121212_data/mapping/hisat2
mkdir -p /home/projects/ku_00015/data/gse121212_data/mapping/hisat2-report
cat ${ID_PATH} | parallel "hisat2 -x /home/projects/ku_00015/data/rna-seq/reference/genome-GRCh38/genome -1  /home/projects/ku_00015/data/gse121212_data/reads/{}_1.fastq -2 /home/projects/ku_00015/data/gse121212_data/reads/{}_2.fastq -S /home/projects/ku_00015/data/gse121212_data/mapping/hisat2/{}.bam --summary-file /home/projects/ku_00015/data/gse121212_data/mapping/hisat2-report/{}-report.summary"

#use qsub to submit mapping task
qsub -W group_list=ku_00015 -A ku_00015 -l nodes=1:ppn=40:thinnode,walltime=24:00:00 -m n /home/projects/ku_00015/people/tuhu/gse121212/doc/RNA-seq-mapping.sh

# feature counting
module load tools subread/1.6.2
mkdir -p /home/projects/ku_00015/people/tuhu/gse121212/counts
featureCounts -T 32 -t exon -g gene_id -a /home/projects/ku_00015/data/rna-seq/reference/gtf/gencode.v34.annotation.gtf -o /home/projects/ku_00015/people/tuhu/gse121212/counts/counts.txt /home/projects/ku_00015/data/gse121212_data/mapping/hisat2/*.bam

qsub -W group_list=ku_00015 -A ku_00015 -l nodes=1:ppn=40:thinnode,walltime=24:00:00 -m n /home/projects/ku_00015/people/tuhu/gse121212/doc/RNA-seq-featureCounts.sh