library(tidyverse)
library(readxl)

# Set raw data paths
gene_path <- "data/SCAN-B_GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv" #transformed gene expression data

clin_path <- "data/SCAN-B_GSE96058/41523_2022_465_MOESM2_ESM.xlsx" #newest clinical data for SCAN-B

ID_links1_path <- "data/SCAN-B_GSE96058/GSE96058-GPL11154_series_matrix.txt" #file #1 for IDs to link clinical to gene expression

ID_links2_path <- "data/SCAN-B_GSE96058/GSE96058-GPL18573_series_matrix.txt" #file #2 for IDs to link clinical to gene expression

# Import data
SB_gene <- read.delim(gene_path, sep = ",", header = T, check.names = F)

SB_clin <- read_xlsx(clin_path, sheet = "SCANB.9206", col_names = T)

ID_links1 <- read.delim(ID_links1_path, sep = "\t", header = T, check.names = F)

ID_links2 <- read.delim(ID_links2_path, sep = "\t", header = T, check.names = F)

# Build key to link clinical data to gene expression and remove replicates

