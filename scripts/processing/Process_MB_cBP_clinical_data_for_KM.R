library(tidyverse)
library(janitor)

# Set raw data paths
MB_clin_path <- "data/METABRIC_cBioPortal/data_clinical_patient.txt"

# Import data
MB_clin_data <- read.delim(MB_clin_path, sep = "\t", header = F)

# Prepare data for KM plot
MB_clin_data <- MB_clin_data |> 
  slice(5:n()) |> #Remove first 4 rows
  janitor::row_to_names(row_number = 1) |> #Set first row as column names
  select(PATIENT_ID, 
         OS_MONTHS, 
         OS_STATUS, 
         CLAUDIN_SUBTYPE, 
         VITAL_STATUS, 
         RFS_STATUS, 
         RFS_MONTHS) |> #Select columns needed for analysis
  mutate_all(~na_if(., '')) |> #Recode empty values as NAs
  na.omit()
  


