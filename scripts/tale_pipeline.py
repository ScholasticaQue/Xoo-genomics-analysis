import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import networkx as nx

# ===============================
# LOAD DATA
# ===============================
df = pd.read_csv("data/all_targets.tsv", sep="\t", comment="#", header=None)

df.columns = [
    "region", "position", "strand", "score",
    "sequence", "pvalue", "rvds", "tale"
]

# ===============================
# FILTER HIGH CONFIDENCE
# ===============================
df = df[df["pvalue"] < 1e-5]

df.to_csv("results/high_conf_targets.tsv", sep="\t", index=False)

print("High confidence sites:", len(df))

# ===============================
# EXTRACT COORDINATES
# ===============================
def parse_region(region):
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)

coords = df["region"].apply(parse_region)
df["chr"] = coords.apply(lambda x: x[0])
df["start"] = coords.apply(lambda x: x[1])
df["end"] = coords.apply(lambda x: x[2])

# approximate binding site position
df["start"] = df["start"] + df["position"] - 1
df["end"] = df["start"] + 20

# ===============================
# SAVE BED FILE
# ===============================
bed = df[["chr", "start", "end", "tale"]]
bed.to_csv("results/targets.bed", sep="\t", header=False, index=False)

print("BED file written")

# ===============================
# BASIC STATISTICS
# ===============================
tale_counts = df["tale"].value_counts()
gene_like_regions = df["region"].nunique()

print("Unique TALEs:", df["tale"].nunique())
print("Unique regions:", gene_like_regions)

tale_counts.to_csv("results/tale_target_counts.csv")

# ===============================
# FIGURE 1 — TALE TARGET DISTRIBUTION
# ===============================
plt.figure()
tale_counts.hist(bins=50)
plt.xlabel("Targets per TALE")
plt.ylabel("Frequency")
plt.title("TALE Target Distribution")
plt.savefig("figures/fig1_tale_distribution.png", dpi=300)

# ===============================
# SIMULATED GENE MAPPING (placeholder)
# replace later with BEDTOOLS output
# ===============================
mapped = pd.read_csv("results/mapped_targets.tsv", sep="\t", header=None)
mapped.columns = ["chr","start","end","tale","gff_chr","source","feature","gff_start","gff_end","score","strand","frame","attributes"]

mapped["gene"] = mapped["attributes"].str.extract('ID=([^;]+)')

df = mapped
# ===============================
# FIGURE 2 — TOP TARGETED GENES
# ===============================
gene_counts = df["gene"].value_counts().head(20)

plt.figure()
gene_counts.plot(kind="bar")
plt.ylabel("Hits")
plt.title("Top Targeted Regions")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("figures/fig2_top_targets.png", dpi=300)

# ===============================
# FIGURE 3 — HEATMAP
# ===============================
import seaborn as sns

top_genes = df["gene"].value_counts().head(30).index
heatmap_df = pd.crosstab(df["tale"], df["gene"])

heatmap_df = heatmap_df[top_genes]

plt.figure(figsize=(10,8))
sns.heatmap(heatmap_df, cmap="viridis")
plt.title("TALE–Target Heatmap")
plt.savefig("figures/fig3_heatmap.png", dpi=300)

# ===============================
# FIGURE 4 — NETWORK
# ===============================
G = nx.Graph()

for _, row in df.iterrows():
    G.add_edge(row["tale"], row["gene"])

plt.figure(figsize=(12,12))
nx.draw(G, node_size=5, with_labels=False)
plt.savefig("figures/fig4_network.png", dpi=300)

print("All figures generated")

# ===============================
# CORE vs ACCESSORY TARGETS
# ===============================
gene_presence = df.groupby("gene")["tale"].nunique()

core = gene_presence[gene_presence > 5]
accessory = gene_presence[gene_presence <= 5]

print("Core targets:", len(core))
print("Accessory targets:", len(accessory))

gene_presence.to_csv("results/gene_presence.csv")
