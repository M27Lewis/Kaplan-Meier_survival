library(tidyverse)

# Set raw data paths
MB_clin_path <- "data/METABRIC_cBioPortal/data_clinical_patient.txt"

# Import data
MB_clin_data <- read.delim(MB_clin_path, sep = "\t", header = F)

