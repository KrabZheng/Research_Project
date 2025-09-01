library(readr)
library(dplyr)
library(ggplot2)

setwd("/Users/SofTicE/Downloads")

df <- read_tsv("GSE81089_FPKM_cufflinks.tsv.gz")

glimpse(data)   

## 1.  Install/load the packages you need  ------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("biomaRt", "dplyr"), ask = FALSE)

library(biomaRt)
library(dplyr)


## 3.  Pull the mapping from Ensembl via biomaRt ------------------------
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = unique(df$Ensembl_gene_id),
    mart       = mart
)

## 4.  Join the mapping back to your data frame ------------------------
df <- df %>%
      left_join(mapping,
                by = c("Ensembl_gene_id" = "ensembl_gene_id")) %>%
      rename(hgnc_symbol = hgnc_symbol)       # keeps the name tidy

ggplot(subset(df, hgnc_symbol == "HNF1B"), aes(x = sample_id, y = HNF1B)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "TP53 expression across samples",
         x = "Sample ID",
         y = "FPKM")

# 1. load needed packages
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. pick out just the three genes
genes_of_interest <- c("GAPDH", "UGT1A10", "NFE2L2")

df_sub <- df %>%
  filter(hgnc_symbol %in% genes_of_interest)

# 3. pivot into long format
df_long <- df_sub %>%
dplyr::select(-Ensembl_gene_id) %>%
  pivot_longer(
    cols      = -hgnc_symbol,
    names_to  = "Sample",
    values_to = "Expression"
  ) %>%
  # make sure the genes show up in the order we want
  mutate(hgnc_symbol = factor(hgnc_symbol, levels = genes_of_interest))

# 4. draw the heatmap
ggplot(df_long, aes(x = Sample, y = hgnc_symbol, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = median(df_long$Expression, na.rm = TRUE),
    name     = "Expr"
  ) +
  # here we compute the sample order by NFE2L2, lowest → highest
  scale_x_discrete(
    limits = df_long %>%
      filter(hgnc_symbol == "NFE2L2") %>%
      arrange(Expression) %>%
      pull(Sample)
  ) +
  labs(
    x     = "Sample (ordered by NFE2L2 expr.)",
    y     = "Gene",
    title = "Expression of HNF1B, UGT1A10 & NFE2L2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )