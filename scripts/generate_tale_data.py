#!/usr/bin/env python3

import os
import glob
import subprocess
import pandas as pd

# ===============================
# CONFIG
# ===============================
GENOME_DIR = "genomes"
OUTDIR = "TALE_project"

RICE_GENOME = "rice_genome.fa"
RICE_GFF = "rice_genes.gff3"
RICE_BED = "rice_genes.bed"
GENOME_FILE = "rice.genome"

THREADS = 8

# ===============================
# HELPER
# ===============================
def run(cmd, step):
    print(f"\n🚀 Running: {step}")
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

# ===============================
# STEP 1 — AnnoTALE
# ===============================
os.makedirs(f"{OUTDIR}/tale_raw", exist_ok=True)

for genome in glob.glob(f"{GENOME_DIR}/*.fna"):
    sample = os.path.basename(genome).replace(".fna", "")
    out = f"{OUTDIR}/tale_raw/{sample}"

    cmd = f"annotateTALEs.py --genome {genome} --outdir {out}"
    run(cmd, f"AnnoTALE ({sample})")

# ===============================
# STEP 2 — Combine TALEs
# ===============================
os.makedirs(f"{OUTDIR}/tale_clusters", exist_ok=True)

run(
    f"cat {OUTDIR}/tale_raw/*/*.faa > {OUTDIR}/tale_clusters/all_tales.faa",
    "Combine TALEs"
)

# ===============================
# STEP 3 — CD-HIT clustering
# ===============================
run(
    f"cd-hit -i {OUTDIR}/tale_clusters/all_tales.faa "
    f"-o {OUTDIR}/tale_clusters/tales_90.faa "
    f"-c 0.9 -n 5 -d 0 -T {THREADS}",
    "CD-HIT clustering"
)

# ===============================
# STEP 4 — Extract RVDs
# ===============================
os.makedirs(f"{OUTDIR}/rvds", exist_ok=True)

rvd_out = f"{OUTDIR}/rvds/preditale_input.tsv"

with open(rvd_out, "w") as out:
    for file in glob.glob(f"{OUTDIR}/tale_raw/*/*.txt"):
        with open(file) as f:
            for line in f:
                if "RVD" in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        out.write(f"{parts[0]}\t{parts[-1]}\n")

print("RVD extraction complete")

# ===============================
# STEP 5 — PrediTALE
# ===============================
os.makedirs(f"{OUTDIR}/pred_targets", exist_ok=True)

run(
    f"PrediTALE -i {rvd_out} -g {RICE_GENOME} "
    f"-o {OUTDIR}/pred_targets/all_targets.tsv",
    "PrediTALE prediction"
)

# ===============================
# STEP 6 — Filter targets
# ===============================
os.makedirs(f"{OUTDIR}/filtered_targets", exist_ok=True)

run(
    f"awk '$6 < 1e-6' {OUTDIR}/pred_targets/all_targets.tsv "
    f"> {OUTDIR}/filtered_targets/high_conf.tsv",
    "Filter high-confidence targets"
)

# ===============================
# STEP 7 — Convert to BED
# ===============================
os.makedirs(f"{OUTDIR}/mapping", exist_ok=True)

df = pd.read_csv(
    f"{OUTDIR}/filtered_targets/high_conf.tsv",
    sep="\t",
    header=None
)

df.columns = [
    "region","position","strand","score",
    "sequence","pvalue","rvds","tale"
]

def parse_region(region):
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)

coords = df["region"].apply(parse_region)

df["chr"] = coords.apply(lambda x: x[0])
df["start"] = coords.apply(lambda x: x[1])
df["end"] = coords.apply(lambda x: x[2])

df["start"] = df["start"] + df["position"] - 1
df["end"] = df["start"] + 20

bed_file = f"{OUTDIR}/mapping/targets.bed"
df[["chr","start","end","tale"]].to_csv(
    bed_file, sep="\t", header=False, index=False
)

print("BED file created")

# ===============================
# STEP 8 — Map to genes
# ===============================
mapped_genes = f"{OUTDIR}/mapping/mapped_targets.tsv"

run(
    f"bedtools intersect -a {bed_file} -b {RICE_GFF} "
    f"-wa -wb > {mapped_genes}",
    "Map to genes"
)

# ===============================
# STEP 9 — Promoter generation
# ===============================
promoters = f"{OUTDIR}/mapping/promoters.bed"

run(
    f"bedtools flank -i {RICE_BED} -g {GENOME_FILE} "
    f"-l 1000 -r 0 -s > {promoters}",
    "Generate promoters"
)

# ===============================
# STEP 10 — Map to promoters
# ===============================
mapped_promoters = f"{OUTDIR}/mapping/mapped_promoters.tsv"

run(
    f"bedtools intersect -a {bed_file} -b {promoters} "
    f"-wa -wb > {mapped_promoters}",
    "Map to promoters"
)

print("\n🎉 FULL PIPELINE COMPLETE")
