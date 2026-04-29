import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

# ===============================
# LOAD DATA
# ===============================
df = pd.read_csv("results/mapped_targets.tsv", sep="\t", header=None)

df.columns = [
    "chr","start","end","tale",
    "gff_chr","source","feature",
    "gff_start","gff_end","score",
    "strand","frame","attributes"
]

# ===============================
# KEEP ONLY GENES
# ===============================
df = df[df["feature"] == "gene"]

# ===============================
# EXTRACT GENE ID + NAME
# ===============================
df["gene_id"] = df["attributes"].str.extract('gene_id=([^;]+)')
df["gene_name"] = df["attributes"].str.extract('Name=([^;]+)')

df["gene"] = df["gene_id"].fillna(df["gene_name"])

df = df.dropna(subset=["gene"])

print("Total mapped interactions:", len(df))

# ===============================
# REMOVE DUPLICATES (CRITICAL)
# ===============================
df = df.drop_duplicates(subset=["tale", "gene"])

print("Unique TALE–gene interactions:", len(df))

# ===============================
# TARGET COUNTS
# ===============================
gene_counts = df["gene"].value_counts()
tale_counts = df["tale"].value_counts()

gene_counts.to_csv("results/gene_target_counts.csv")
tale_counts.to_csv("results/tale_target_counts.csv")

# ===============================
# CORE vs ACCESSORY (REAL)
# ===============================
gene_presence = df.groupby("gene")["tale"].nunique()

core = gene_presence[gene_presence >= 5]
accessory = gene_presence[gene_presence < 5]

print("Core genes:", len(core))
print("Accessory genes:", len(accessory))

gene_presence.to_csv("results/gene_presence.csv")

# ===============================
#  FIGURE 1 — TOP TARGET GENES
# ===============================
top = gene_counts.head(20)

plt.figure(figsize=(10,5))
top.plot(kind="bar")
plt.ylabel("Number of TALEs targeting gene")
plt.title("Top TALE Target Genes")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("figures/fig1_top_genes.png", dpi=300)

# ===============================
#  FIGURE 2 — DISTRIBUTION
# ===============================
plt.figure()
gene_presence.hist(bins=50)
plt.xlabel("Number of TALEs per gene")
plt.ylabel("Frequency")
plt.title("Gene Targeting Distribution")
plt.savefig("figures/fig2_distribution.png", dpi=300)

# ===============================
#  FIGURE 3 — HEATMAP
# ===============================
top_genes = gene_counts.head(30).index

heatmap_df = pd.crosstab(df["tale"], df["gene"])
heatmap_df = heatmap_df[top_genes]

# reduce noise
heatmap_df = heatmap_df.loc[
    heatmap_df.sum(axis=1) > 2,
    heatmap_df.sum(axis=0) > 2
]

plt.figure(figsize=(10,8))
sns.heatmap(heatmap_df, cmap="viridis")
plt.title("TALE–Gene Interaction Heatmap")
plt.savefig("figures/fig3_heatmap.png", dpi=300)

# ===============================
#  FIGURE 4 — NETWORK
# ===============================
G = nx.Graph()

for _, row in df.iterrows():
    if row["gene"] in top_genes:
        G.add_edge(row["tale"], row["gene"])

plt.figure(figsize=(12,12))
nx.draw(G, node_size=15, with_labels=False)
plt.savefig("figures/fig4_network.png", dpi=300)

print("Final gene-level analysis complete")
