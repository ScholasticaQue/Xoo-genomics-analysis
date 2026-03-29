library(ggplot2)

data <- read.csv(snakemake@input[[1]], row.names=1)

pca <- prcomp(data, scale.=TRUE)

df <- as.data.frame(pca$x)

pdf(snakemake@output[[1]])
ggplot(df, aes(PC1, PC2)) +
  geom_point(size=3) +
  theme_minimal()
dev.off()
