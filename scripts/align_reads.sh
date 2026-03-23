#!/bin/bash

REF="/home/yourname/xoo_project_data/reference.fna"
DATA="/home/yourname/xoo_project_data/fastq"

for sample in $DATA/*_R1.fastq; do
    base=$(basename $sample _R1.fastq)

    bwa mem $REF \
        ${DATA}/${base}_R1.fastq \
        ${DATA}/${base}_R2.fastq | \
    samtools sort -o /home/yourname/xoo_project_data/bam/${base}.bam
done
