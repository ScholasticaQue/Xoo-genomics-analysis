import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import os
import re

# ===============================
# HELPER FUNCTION
# ===============================
def extract_strain(tale_name):
    match = re.search(r'_(.*?)-', tale_name)
    return match.group(1) if match else "unknown"

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
# CLEAN
# ===============================
df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
df["score"] = pd.to_numeric(df["score"], errors="coerce")
df = df.dropna(subset=["pvalue"])

print("After cleaning:", len(df))

# ===============================
# FILTER
# ===============================
df = df[df["pvalue"] < 1e-6].copy()

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
# BINDING LENGTH
# ===============================
df["length"] = df["rvds"].str.count("-") + 1
df["start"] = df["start"] + df["position"] - 1
df["end"] = df["start"] + df["length"]

# ===============================
# SAVE BED
# ===============================
df[["chr","start","end","tale"]].to_csv(
    "results/targets.bed", sep="\t", header=False, index=False
)

print("BED file written")

# ===============================
# LOAD PROMOTER TARGETS
# ===============================
mapped = pd.read_csv("results/mapped_promoters.tsv", sep="\t", header=None)

mapped.columns = [
    "chr","start","end","tale",
    "p_chr","p_start","p_end",
    "gene","dot","strand"
]

mapped = mapped.dropna(subset=["gene"]).copy()
mapped = mapped.drop_duplicates(subset=["tale","gene"])

print("Unique TALE–gene interactions:", len(mapped))

# ===============================
# ADD ANNOTATION
# ===============================
annot = pd.read_csv("data/annotation.tsv", sep="\t", header=None)
annot.columns = ["gene", "gene_name"]

mapped = mapped.merge(annot, on="gene", how="left")

print("Annotated genes:", mapped["gene_name"].notna().sum())

# ===============================
# ADD STRAIN INFO
# ===============================
mapped["strain"] = mapped["tale"].apply(extract_strain)

# ===============================
# BASIC STATS
# ===============================
gene_counts = mapped["gene"].value_counts()
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

core.to_csv("results/core_genes_list.csv")
gene_presence.to_csv("results/gene_presence.csv")

# ===============================
# KEY GENES
# ===============================
keywords = ["SWEET", "NAC", "WRKY", "BZIP"]

key_genes_df = mapped[
    mapped["gene_name"].str.contains("|".join(keywords), case=False, na=False)
].copy()

key_genes_df.to_csv("results/key_gene_interactions.tsv", sep="\t", index=False)

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
# 🔥 NEW: STRAIN × GENE MATRIX
# ===============================
strain_gene_matrix = pd.crosstab(mapped["strain"], mapped["gene"])

strain_gene_matrix.to_csv("results/strain_gene_matrix.csv")

# ===============================
# 🔥 NEW: TALE × GENE MATRIX
# ===============================
tale_gene_matrix = pd.crosstab(mapped["tale"], mapped["gene"])
tale_gene_matrix.to_csv("results/tale_gene_matrix.csv")

# ===============================
# 🔥 FULL TARGET ANALYSIS
# ===============================
full_targets = pd.read_csv("results/mapped_targets.tsv", sep="\t", header=None)

print("mapped_targets columns:", full_targets.shape[1])

if full_targets.shape[1] == 12:
    full_targets.columns = [
        "chr","start","end","tale",
        "p_chr","source","feature","g_start","g_end",
        "dot","strand","info"
    ]
elif full_targets.shape[1] == 13:
    full_targets.columns = [
        "chr","start","end","tale",
        "p_chr","source","feature","g_start","g_end",
        "score","strand","phase","info"
    ]
else:
    raise ValueError("Unexpected columns")

# extract gene id
full_targets["gene"] = full_targets["info"].str.extract(r"gene_id=([^;]+)")

# merge annotation
full_targets = full_targets.merge(annot, on="gene", how="left")

# ===============================
# SWEET ANALYSIS
# ===============================
sweet_full = full_targets[
    full_targets["gene_name"].str.contains("SWEET", case=False, na=False)
].copy()

print("FULL SWEET interactions:", len(sweet_full))
print("Unique SWEET genes:", sweet_full["gene"].nunique())

sweet_full["strain"] = sweet_full["tale"].apply(extract_strain)

# counts
tale_sweet_counts = sweet_full["tale"].value_counts()
strain_counts = sweet_full["strain"].value_counts()

tale_sweet_counts.to_csv("results/tale_sweet_counts_full.csv")
strain_counts.to_csv("results/strain_sweet_counts_full.csv")

# ===============================
# 🔥 NEW: SWEET NETWORK TABLE
# ===============================
sweet_network = sweet_full[["strain","tale","gene","gene_name"]]
sweet_network.to_csv("results/sweet_network_table.tsv", sep="\t", index=False)

# ===============================
# FIGURES
# ===============================
plt.figure()
gene_presence.hist(bins=50)
plt.title("Distribution of TALE Targeting per Gene")
plt.savefig("figures/fig1_distribution.png", dpi=300)

# ===============================
# HEATMAP
# ===============================
top30 = gene_counts.head(30).index

heatmap_df = pd.crosstab(mapped["tale"], mapped["gene"])
heatmap_df = heatmap_df[top30]

heatmap_df = heatmap_df.loc[
    heatmap_df.sum(axis=1) > 2,
    heatmap_df.sum(axis=0) > 2
]

plt.figure(figsize=(10,8))
sns.heatmap(heatmap_df)
plt.title("TALE–Gene Interaction Heatmap")
plt.savefig("figures/fig3_heatmap.png", dpi=300)

# ===============================
# NETWORK
# ===============================
G = nx.Graph()

for _, row in mapped.iterrows():
    if row["gene"] in top30:
        G.add_edge(row["tale"], row["gene"])

pos = nx.spring_layout(G, seed=42)

plt.figure(figsize=(10,10))
nx.draw(G, pos, node_size=20)
plt.title("TALE–Gene Network")
plt.savefig("figures/fig4_network.png", dpi=300)

print("Pipeline complete")
