#!/bin/bash

THREADS=4
GENOME_SIZE=5m

# Pre-subsampled Illumina (DO THIS ONCE BEFORE RUNNING SCRIPT)
ILLUMINA_R1=BX01_1_sub.fastq
ILLUMINA_R2=BX01_2_sub.fastq

# Set Java memory
export _JAVA_OPTIONS="-Xmx6G"

# Loop through all Nanopore datasets
for READS in *.fastq.gz
do
    SAMPLE=$(basename $READS .fastq.gz)

    echo "==============================="
    echo "Processing $SAMPLE"
    echo "==============================="

    mkdir -p $SAMPLE
    cd $SAMPLE

    LOG="${SAMPLE}.log"

    # -----------------------------
    # STEP 1 — Flye (SKIP IF DONE)
    # -----------------------------
    if [ ! -f flye/assembly.fasta ]; then
        echo "[Flye] Running assembly..." | tee -a $LOG
        flye --nano-raw ../$READS \
             --out-dir flye \
             --genome-size $GENOME_SIZE \
             --threads $THREADS >> $LOG 2>&1
    else
        echo "[Flye] Skipped (already exists)" | tee -a $LOG
    fi

    ASSEMBLY=flye/assembly.fasta

    # -----------------------------
    # STEP 2 — Racon
    # -----------------------------
    if [ ! -f racon2.fasta ]; then
        echo "[Racon] Round 1" | tee -a $LOG
        minimap2 -t $THREADS -x map-ont $ASSEMBLY ../$READS > aln1.paf
        racon -t $THREADS ../$READS aln1.paf $ASSEMBLY > racon1.fasta

        echo "[Racon] Round 2" | tee -a $LOG
        minimap2 -t $THREADS -x map-ont racon1.fasta ../$READS > aln2.paf
        racon -t $THREADS ../$READS aln2.paf racon1.fasta > racon2.fasta

        rm aln1.paf aln2.paf racon1.fasta
    else
        echo "[Racon] Skipped" | tee -a $LOG
    fi

    # -----------------------------
    # STEP 3 — Pilon
    # -----------------------------
    if [ ! -f final_polished.fasta ]; then
        echo "[Pilon] Running polishing..." | tee -a $LOG

        bwa index racon2.fasta

        bwa mem -t $THREADS racon2.fasta ../$ILLUMINA_R1 ../$ILLUMINA_R2 > aln.sam

        samtools view -Sb aln.sam | samtools sort -o aln.sorted.bam
        samtools index aln.sorted.bam

        pilon --genome racon2.fasta \
              --bam aln.sorted.bam \
              --output final_polished >> $LOG 2>&1

        rm aln.sam
    else
        echo "[Pilon] Skipped" | tee -a $LOG
    fi

    # -----------------------------
    # STEP 4 — BUSCO
    # -----------------------------
    if [ ! -d busco ]; then
        echo "[BUSCO] Running..." | tee -a $LOG

        busco \
            -i final_polished.fasta \
            -o busco \
            -l xanthomonadales_odb10 \
            -m genome \
            -c $THREADS \
            -f >> $LOG 2>&1
    else
        echo "[BUSCO] Skipped" | tee -a $LOG
    fi

    cd ..
done

echo "=== ALL STRAINS COMPLETED ==="
