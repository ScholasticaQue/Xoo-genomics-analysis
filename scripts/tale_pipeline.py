import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import os

# ===============================
# SETUP
# ===============================
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# ===============================
# LOAD RAW TARGET DATA
# ===============================
df = pd.read_csv(
    "data/all_targets.tsv",
    sep="\t",
    comment="#",
    header=None,
    names=["region","position","strand","score","sequence","pvalue","rvds","tale"]
)

print("Total raw rows:", len(df))

# ===============================
# CLEAN TYPES
# ===============================
df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
df["score"] = pd.to_numeric(df["score"], errors="coerce")

df = df.dropna(subset=["pvalue"])

print("After cleaning:", len(df))

# ===============================
# FILTER HIGH CONFIDENCE
# ===============================
df = df[df["pvalue"] < 1e-6]

print("High-confidence sites:", len(df))
print("Max p-value:", df["pvalue"].max())

# ===============================
# PARSE COORDINATES
# ===============================
def parse_region(region):
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)

coords = df["region"].apply(parse_region)

df["chr"] = coords.apply(lambda x: x[0])
df["start"] = coords.apply(lambda x: x[1])

# ===============================
# FIX: REAL TALE BINDING LENGTH
# ===============================
df["length"] = df["rvds"].str.count("-") + 1

df["start"] = df["start"] + df["position"] - 1
df["end"] = df["start"] + df["length"]

# ===============================
# SAVE BED FILE
# ===============================
df[["chr","start","end","tale"]].to_csv(
    "results/targets.bed", sep="\t", header=False, index=False
)

print("BED file written")

mapped = pd.read_csv("results/mapped_promoters.tsv", sep="\t", header=None)

mapped.columns = [
    "chr","start","end","tale",
    "p_chr","p_start","p_end",
    "gene","dot","strand"
]

# keep only valid genes
mapped = mapped.dropna(subset=["gene"])

# remove duplicates
mapped = mapped.drop_duplicates(subset=["tale","gene"])

print("Unique TALE–gene interactions:", len(mapped))

# ===============================
# BASIC STATS
# ===============================
gene_counts = mapped["gene"].value_counts()
tale_counts = mapped["tale"].value_counts()

print("Unique genes:", mapped["gene"].nunique())
print("Unique TALEs:", mapped["tale"].nunique())

# ===============================
# CORE vs ACCESSORY (IMPROVED)
# ===============================
n_tales = mapped["tale"].nunique()

# dynamic threshold (10% of TALEs OR minimum 5)
core_threshold = max(5, int(0.1 * n_tales))

gene_presence = mapped.groupby("gene")["tale"].nunique()

core = gene_presence[gene_presence >= core_threshold]
accessory = gene_presence[gene_presence < core_threshold]

print(f"Core threshold: {core_threshold}")
print("Core genes:", len(core))
print("Accessory genes:", len(accessory))

gene_presence.to_csv("results/gene_presence.csv")

# ===============================
# FIGURE 1 — DISTRIBUTION
# ===============================
plt.figure()
gene_presence.hist(bins=50)
plt.xlabel("Number of TALEs targeting gene")
plt.ylabel("Frequency")
plt.title("Distribution of TALE Targeting per Gene")
plt.savefig("figures/fig1_distribution.png", dpi=300)

# ===============================

# FIGURE 2 — TOP TARGET GENES
# ===============================
top_genes = gene_counts.head(20)

plt.figure(figsize=(10,5))
top_genes.plot(kind="bar")
plt.ylabel("Number of TALEs")
plt.xticks(rotation=90)
plt.title("Top 20 Most Targeted Genes")
plt.tight_layout()
plt.savefig("figures/fig2_top_genes.png", dpi=300)

# ===============================
# FIGURE 3 — HEATMAP (CLEANED)
# ===============================
top30 = gene_counts.head(30).index

heatmap_df = pd.crosstab(mapped["tale"], mapped["gene"])
heatmap_df = heatmap_df[top30]

# remove sparse rows/cols
heatmap_df = heatmap_df.loc[
    heatmap_df.sum(axis=1) > 2,
    heatmap_df.sum(axis=0) > 2
]

plt.figure(figsize=(10,8))
sns.heatmap(heatmap_df, cmap="viridis")
plt.title("TALE–Gene Interaction Heatmap")
plt.savefig("figures/fig3_heatmap.png", dpi=300)

# ===============================
# FIGURE 4 — NETWORK (FILTERED)
# ===============================
# ===============================
# FIGURE 4 — NETWORK (FIXED)
# ===============================
G = nx.Graph()

for _, row in mapped.iterrows():
    if row["gene"] in top30:
        G.add_edge(row["tale"], row["gene"])

# layout
pos = nx.spring_layout(G, k=0.5, seed=42)

# separate nodes safely
tales = [n for n in G.nodes if "tempTALE" in n]
genes = [n for n in G.nodes if n not in tales]

plt.figure(figsize=(10,10))

# TALE nodes
nx.draw_networkx_nodes(
    G, pos,
    nodelist=tales,
    node_size=30,
    node_color="skyblue",
    label="TALEs"
)

# gene nodes
nx.draw_networkx_nodes(
    G, pos,
    nodelist=genes,
    node_size=80,
    node_color="orange",
    label="Genes"
)

# edges
nx.draw_networkx_edges(G, pos, alpha=0.3)

plt.title("TALE–Gene Interaction Network (Top Genes)")
plt.legend()
plt.axis("off")

plt.savefig("figures/fig4_network.png", dpi=300)
# ===============================
# SAVE FINAL TABLE
# ===============================
mapped.to_csv("results/final_tale_gene_interactions.tsv", sep="\t", index=False)

print("Pipeline complete")

