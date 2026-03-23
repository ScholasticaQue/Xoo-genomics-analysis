#!/bin/bash

REF="/path/to/reference.fna"
DATA="/path/to/your/bam_files"

for bam in $DATA/*.bam; do
    sample=$(basename $bam .bam)

    snippy --cpus 8 \
           --outdir results/snippy_$sample \
           --ref $REF \
           --bam $bam
done
