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
# FIX BINDING LENGTH
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

# ===============================
# LOAD PROMOTER-MAPPED TARGETS
# ===============================
mapped = pd.read_csv("results/mapped_promoters.tsv", sep="\t", header=None)

mapped.columns = [
    "chr","start","end","tale",
    "p_chr","p_start","p_end",
    "gene","dot","strand"
]

mapped = mapped.dropna(subset=["gene"])
mapped = mapped.drop_duplicates(subset=["tale","gene"])

print("Unique TALE–gene interactions:", len(mapped))

# ===============================
# ADD GENE ANNOTATION (CRITICAL)
# ===============================
annot = pd.read_csv("data/annotation.tsv", sep="\t", header=None)
annot.columns = ["gene", "gene_name"]

mapped = mapped.merge(annot, on="gene", how="left")

print("Annotated genes:", mapped["gene_name"].notna().sum())

# ===============================
# BASIC STATS
# ===============================
gene_counts = mapped["gene"].value_counts()
tale_counts = mapped["tale"].value_counts()

print("Unique genes:", mapped["gene"].nunique())
print("Unique TALEs:", mapped["tale"].nunique())

# ===============================
# CORE vs ACCESSORY
# ===============================
n_tales = mapped["tale"].nunique()
core_threshold = max(5, int(0.1 * n_tales))

gene_presence = mapped.groupby("gene")["tale"].nunique()

core = gene_presence[gene_presence >= core_threshold]
accessory = gene_presence[gene_presence < core_threshold]

print(f"Core threshold: {core_threshold}")
print("Core genes:", len(core))
print("Accessory genes:", len(accessory))

# SAVE CORE DATA
core.to_csv("results/core_genes_list.csv")
gene_presence.to_csv("results/gene_presence.csv")

# ===============================
# KEY GENE EXTRACTION (FIXED)
# ===============================
keywords = ["SWEET", "NAC", "WRKY", "BZIP"]

key_genes_df = mapped[
    mapped["gene_name"].str.contains("|".join(keywords), case=False, na=False)
]

key_genes_df.to_csv("results/key_gene_interactions.tsv", sep="\t", index=False)

key_gene_counts = key_genes_df["gene"].value_counts()
key_gene_counts.to_csv("results/key_gene_counts.csv")

print("Key gene interactions:", len(key_genes_df))
print("Unique key genes:", key_genes_df["gene"].nunique())

# ===============================
# CORE + KEY OVERLAP
# ===============================
core_genes = set(core.index)
key_genes = set(key_genes_df["gene"])

core_key = core_genes.intersection(key_genes)

pd.Series(list(core_key)).to_csv("results/core_key_genes.txt", index=False)

print("Core key genes:", len(core_key))

# ===============================
# CORE GENE CATEGORIZATION
# ===============================
core_df = pd.DataFrame({"gene": list(core_genes)})
core_df["category"] = "unknown"

for k in keywords:
    core_df.loc[
        core_df["gene"].isin(
            mapped[mapped["gene_name"].str.contains(k, case=False, na=False)]["gene"]
        ),
        "category"
    ] = k

core_df.to_csv("results/core_gene_categories.csv", index=False)

# ===============================
# TOP CORE GENES
# ===============================
core_ranked = gene_presence.loc[core.index].sort_values(ascending=False)
core_ranked.head(20).to_csv("results/top_core_genes.csv")

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
# FIGURE 2 — TOP GENES
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
# FIGURE 3 — HEATMAP
# ===============================
top30 = gene_counts.head(30).index

heatmap_df = pd.crosstab(mapped["tale"], mapped["gene"])
heatmap_df = heatmap_df[top30]

heatmap_df = heatmap_df.loc[
    heatmap_df.sum(axis=1) > 2,
    heatmap_df.sum(axis=0) > 2
]

plt.figure(figsize=(10,8))
sns.heatmap(heatmap_df, cmap="viridis")
plt.title("TALE–Gene Interaction Heatmap")
plt.savefig("figures/fig3_heatmap.png", dpi=300)

# ===============================
# FIGURE 4 — NETWORK (FIXED)
# ===============================
# ===============================
# FIGURE 4 — NETWORK (FIXED SAFE)
# ===============================
G = nx.Graph()

for _, row in mapped.iterrows():
    if row["gene"] in top30:
        G.add_edge(row["tale"], row["gene"])

pos = nx.spring_layout(G, k=0.6, iterations=100, seed=42)

tales = [n for n in G.nodes if "tempTALE" in n]
genes = [n for n in G.nodes if n not in tales]

# SAFE filtering (FIX)
graph_nodes = set(G.nodes())

key_nodes = [g for g in key_genes_df["gene"].unique() if g in graph_nodes]
core_nodes = [g for g in core_genes if g in graph_nodes]

plt.figure(figsize=(10,10))

# TALE nodes
nx.draw_networkx_nodes(G, pos, nodelist=tales, node_size=30)

# normal genes
nx.draw_networkx_nodes(
    G, pos,
    nodelist=[g for g in genes if g not in key_nodes],
    node_size=60
)

# key genes (red)
nx.draw_networkx_nodes(
    G, pos,
    nodelist=key_nodes,
    node_size=120
)

# core genes (outlined)
nx.draw_networkx_nodes(
    G, pos,
    nodelist=core_nodes,
    node_size=150,
    node_color="none",
    edgecolors="black",
    linewidths=1.5
)

nx.draw_networkx_edges(G, pos, alpha=0.3)

plt.title("TALE–Gene Interaction Network (Biological Highlights)")
plt.axis("off")

plt.savefig("figures/fig4_network.png", dpi=300)
# ===============================
# SAVE FINAL TABLE
# ===============================
mapped.to_csv("results/final_tale_gene_interactions.tsv", sep="\t", index=False)

print("Pipeline complete")
