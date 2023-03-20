
MB_path <- "data/METABRIC_cBioPortal/data_mrna_agilent_microarray.txt"

MB_z_path <- "data/METABRIC_cBioPortal/data_mrna_agilent_microarray_zscores_ref_all_samples.txt"

MB_data <- read.delim(MB_path, sep = "\t", header = T)

MB_z_data <- read.delim(MB_z_path, sep = "\t", header = T)

