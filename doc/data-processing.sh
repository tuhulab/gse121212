#!/bin/bash

# download data
module load sratoolkit/2.9.6-1
prefetch --option-file /home/people/tuhu/Downloads/SRR_Acc_List.txt
fastq-dump --gzip /home/projects/ku_00015/data/sratool_data/sra/SRR8052*.sra # convert sra file to fastq file

# load module
echo "Step 0 - Load modules".
module load tools fastp/0.19.4 parallel/20200522 anaconda2/4.4.0 hisat2/2.2.0 subread/1.6.2

# initiate - ids
TEST_ID_PATH=/home/projects/ku_00015/people/tuhu/multiomics-ad-transcriptomics/data/ids/test_ids.txt
ID_PATH=/home/projects/ku_00015/people/tuhu/multiomics-ad-transcriptomics/data/ids/ids_exclude87.txt

# QC(1) - fastqc
echo "Step 0 - FastQC"
module load perl/5.30.2 jdk/14 fastqc/0.11.8 parallel/20200522
cat ${ID_PATH} |  parallel "fastqc /home/projects/ku_00015/data/rna-seq/read/{}_1.fastq.gz /home/projects/ku_00015/data/rna-seq/read/{}_2.fastq.gz -o /home/projects/ku_00015/data/rna-seq/qc/"
cd /home/projects/ku_00015/data/rna-seq
fastqc read/*.fastq -o qc/

# QC(2) - multiqc
module unload fastqc/0.11.8 perl/5.30.2
module load anaconda3/4.4.0
cd /home/projects/ku_00015/data/rna-seq
multiqc qc/*_fastqc.zip -o qc/

# QC(3) - multi-multi QC

# trimming
echo "Step 1 - Trimming"
mkdir -p /home/projects/ku_00015/data/rna-seq/trimmed/fastp
mkdir -p /home/projects/ku_00015/data/rna-seq/trimmed/fastp-report
cat ${ID_PATH} |  parallel "fastp -i /home/projects/ku_00015/data/rna-seq/read/{}_1.fastq.gz -I /home/projects/ku_00015/data/rna-seq/read/{}_2.fastq.gz -o /home/projects/ku_00015/data/rna-seq/trimmed/fastp/{}_1.fastq.gz -O /home/projects/ku_00015/data/rna-seq/trimmed/fastp/{}_2.fastq.gz -w 1 -h /home/projects/ku_00015/data/rna-seq/trimmed/fastp-report/{}-report.html"

# mapping
echo "Step 2 - Mapping"
mkdir -p /home/projects/ku_00015/data/rna-seq/mapping/hisat2
mkdir -p /home/projects/ku_00015/data/rna-seq/mapping/hisat2-report
cat ${ID_PATH} | parallel "hisat2 -x /home/projects/ku_00015/data/rna-seq/reference/genome-GRCh38/genome -1  /home/projects/ku_00015/data/rna-seq/trimmed/fastp/{}_1.fastq.gz -2 /home/projects/ku_00015/data/rna-seq/trimmed/fastp/{}_2.fastq.gz -S /home/projects/ku_00015/data/rna-seq/mapping/hisat2/{}.bam --summary-file /home/projects/ku_00015/data/rna-seq/mapping/hisat2-report/{}-report.summary"

hisat2 -p 20 -x /home/projects/ku_00015/data/rna-seq/reference/genome-GRCh38/genome -1  /home/projects/ku_00015/data/rna-seq/trimmed/fastp/NG-23827_01_CO_34_BI_NN_02_G3_lib390556_6741_1_1.fastq.gz -2 /home/projects/ku_00015/data/rna-seq/trimmed/fastp/NG-23827_01_CO_34_BI_NN_02_G3_lib390556_6741_1_2.fastq.gz -S /home/projects/ku_00015/data/rna-seq/mapping/hisat2/NG-23827_01_CO_34_BI_NN_02_G3_lib390556_6741_1.bam --summary-file /home/projects/ku_00015/data/rna-seq/mapping/hisat2-report/NG-23827_01_CO_34_BI_NN_02_G3_lib390556_6741_1-report.summary

# mapping - sort 
mkdir -p /home/projects/ku_00015/data/rna-seq/mapping/hisat2-sort-bam
module load parallel/20200522 samtools/1.9
cat ${TEST_ID_PATH} | parallel "samtools sort -o /home/projects/ku_00015/data/rna-seq/mapping/hisat2-sort-bam/{}.bam /home/projects/ku_00015/data/rna-seq/mapping/hisat2/{}.bam"

# Mapping stat by picard tools
module load java/1.8.0 picard-tools/2.20.2 parallel/20200522 gcc/9.3.0 intel/perflibs/2019_update5 R/3.6.1
mkdir -p /home/projects/ku_00015/data/rna-seq/mapping/hisat2-picard-CollectInsertSizeMetrics
cat ${TEST_ID_PATH} | parallel "java -jar /services/tools/picard-tools/2.20.2/picard.jar CollectInsertSizeMetrics \I=/home/projects/ku_00015/data/rna-seq/mapping/hisat2/{}.bam \O=/home/projects/ku_00015/data/rna-seq/mapping/hisat2-picard-CollectInsertSizeMetrics/{}.txt \H=/home/projects/ku_00015/data/rna-seq/mapping/hisat2-picard-CollectInsertSizeMetrics/{}.pdf"

# feature counting
mkdir -p /home/projects/ku_00015/data/rna-seq/counts
featureCounts -T 32 -t exon -g gene_id -a /home/projects/ku_00015/data/rna-seq/reference/gtf/hg38.ensGene.gtf -o /home/projects/ku_00015/data/rna-seq/counts/counts.txt /home/projects/ku_00015/data/rna-seq/mapping/hisat2/*.bam

# test one mapping tool