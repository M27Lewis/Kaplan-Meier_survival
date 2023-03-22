library(tidyverse)

# Set raw data paths
MB_path <- "data/METABRIC_cBioPortal/data_mrna_agilent_microarray.txt"

MB_z_path <- "data/METABRIC_cBioPortal/data_mrna_agilent_microarray_zscores_ref_all_samples.txt" #Z-scored data from cBioPortal

# Import data
MB_data <- read.delim(MB_path, sep = "\t", header = T, check.names = F)

MB_z_data <- read.table(MB_z_path, sep = "\t", header = T, check.names = F)

# Convert gene names column to rownames and remove Entrez IDs
MB_data <- MB_data |> 
  column_to_rownames(var = "Hugo_Symbol") |> 
  select(!Entrez_Gene_Id)

MB_z_data <- MB_z_data |>
  column_to_rownames(var = "Hugo_Symbol") |> 
  select(!Entrez_Gene_Id)


# Export processed data as RDS file
saveRDS(MB_data, "data/METABRIC_cBioPortal/METABRIC_cBioPortal_raw_data.rds")

saveRDS(MB_z_data, "data/METABRIC_cBioPortal/METABRIC_cBioPortal_z-scored_data.rds")
