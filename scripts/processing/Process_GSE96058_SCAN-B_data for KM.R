library(tidyverse)
library(readxl)

# Set raw data paths
gene_path <- "data/SCAN-B_GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv" #transformed gene expression data

clin_path <- "data/SCAN-B_GSE96058/41523_2022_465_MOESM2_ESM.xlsx" #newest clinical data for SCAN-B

ID_links1_path <- "data/SCAN-B_GSE96058/GSE96058-GPL11154_series_matrix_reduced.txt" #file #1 for IDs to link clinical to gene expression

ID_links2_path <- "data/SCAN-B_GSE96058/GSE96058-GPL18573_series_matrix_reduced.txt" #file #2 for IDs to link clinical to gene expression

# Import data
SB_gene <- read.delim(gene_path, sep = ",", header = T, check.names = F)

SB_clin <- read_xlsx(clin_path, sheet = "SCANB.9206", col_names = T)

ID_links1 <- read.delim(ID_links1_path, sep = "\t", header = F, check.names = F)

ID_links2 <- read.delim(ID_links2_path, sep = "\t", header = F, check.names = F)

# Build key to link clinical data to gene expression and remove replicates
ID_links2 <- ID_links2 |> 
  select(-V1) #remove the first column from links2

full_links <- cbind(ID_links1, ID_links2) #combine the two datasets together

full_links$V1 <- make.unique(full_links$V1, sep = "-") #make V1 column values unique for changing to rownames

full_links <- column_to_rownames(full_links, var = "V1") #make column V1 the rownames

full_links_clean <- full_links[,-grep("repl", full_links[1,])] #remove columns that indicate replicates in the sample title row

rownames(full_links_clean) <- gsub("!", "", rownames(full_links_clean)) #remove !s from rownames

full_links_clean <- full_links_clean[order(as.numeric(substring(full_links_clean[1,], 2)), na.last = F)] #reorder columns so the F-ID numbers are in order

full_links_clean <- t(full_links_clean) #transpose the dataframe
row.names(full_links_clean) <- NULL #remove the rownames that are distracting

full_links_clean <- as.data.frame(full_links_clean) |> 
  rename(scanb_external_id = Sample_characteristics_ch1) #rename the scanb id column

full_links_clean$scanb_external_id <- gsub("scan-b external id: ", "", full_links_clean$scanb_external_id) #remove extra text from scanb id column so I can match to other datasets

full_links_clean$short_id <- substr(full_links_clean$scanb_external_id, 1, 23)

id_key <- full_links_clean |> 
  select(Sample_title, short_id, scanb_external_id) #make simple ID key file

# Modify the clinical data and append it to the ID links
SB_clin <- SB_clin |> 
  mutate(scanb_external_id = paste(Patient, Case, GEX.assay, sep = ".")) #create scanb id column by combining patient, case, and assay entries together

SB_clin$scanb_external_id <- paste(substring(SB_clin$scanb_external_id, 1, nchar(SB_clin$scanb_external_id) - 5), 
                                   substring(SB_clin$scanb_external_id, nchar(SB_clin$scanb_external_id) - 3, nchar(SB_clin$scanb_external_id)), sep = "") #remove the 5th from last character (#2) from each entry in scanb ID so they match the key IDs

SB_clin <- SB_clin |> 
  mutate(short_id = paste(Patient, Case, Sample, sep = ".")) #create scanb id column by combining patient, case, and assay entries together

full_clin <- left_join(id_key, SB_clin, by = "short_id") #combine ID key with clinical data, based on the ID key

full_clin <- full_clin |> 
  distinct(Sample_title, .keep_all = T) #remove duplicated sample IDs that were rerun, but have the same clinical data in them

red_clin <- full_clin |> 
  select(Sample_title, short_id, scanb_external_id.x, scanb_external_id.y, `Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)`, Size.mm, TreatGroup, ClinGroup, SSP.PAM50, SSP.Subtype, NCN.PAM50, NCN.Subtype, DRFi_days, DRFi_event, OS_days, OS_event, RFi_days, RFi_event, BCFi_days, BCFi_event)

red_clin <- red_clin |> 
  drop_na(scanb_external_id.y) |> #remove patients with NAs in the scanb_external_id.y column, having no clinical data. 14 patients
  drop_na(OS_days) #remove patients with NAs in all survival data, thus not useful for KM. 1 patient.

saveRDS(red_clin, "data/SCAN-B_GSE96058/SB_reduced_clinical_data.rds")

# Subset the gene expression data to prepare for KM analysis
clin_ids <- unique(red_clin$Sample_title) #list of patient IDs that I can use to subset the gene expression data

rownames(SB_gene) <- SB_gene[,1] #assign first column (gene names) as rownames

SB_gene <- SB_gene[,-1]

SB_gene <- as.data.frame(SB_gene) |> 
  select(any_of(clin_ids))

sb_genez <- as.data.frame(scale(SB_gene))

saveRDS(sb_genez, "data/SCAN-B_GSE96058/SB_subset_gene_expression_z-scored.rds")
