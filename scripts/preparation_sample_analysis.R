library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)

## Data is already processed -> log2() transformation

raw_metab <- readRDS("data/SummarizedExperiment_extraction_method_cohort_metabolomics.RDS")
raw_prot <- readRDS("data/SummarizedExperiment_extraction_method_cohort_proteomics.RDS")


# annotation tumor and healthy tissue
# raw_metab
metab_count <- as.data.frame(raw_metab@assays@data@listData[[1]])
colnames(metab_count) <- raw_metab@colData@listData[["Sample.Description"]]

metab_count <- unique(metab_count)
metabolomics_TUvsNG <- as.data.frame(metab_count)


# raw_prot
prot_count <- as.data.frame(raw_prot@assays@data@listData[[1]])
colnames(prot_count) <- raw_prot@colData@listData[["Pseudonym"]]
prot_count <- tibble::rownames_to_column(prot_count, "ID")


prot_count <- prot_count %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")

prot_count <- prot_count[prot_count$ID != "NA",]
prot_count <- unique(prot_count)

proteomics_TUvsNG <- as.data.frame(prot_count[,-1])
rownames(proteomics_TUvsNG) <- prot_count$ID

# additionally differentiate between autoSP3 and MTBE
autoSP3orMTBE <- raw_prot@colData@listData[["Processing"]]
sum(autoSP3orMTBE == "MTBE_SP3")
# first 18 are MTBE_SP3, autoSP3 are only 17
# some samples do not have healthy or cancerous partner with sample preparation
# Thereby, these will be combined so it ends up with 20 samples that are pair-wise combined just like in metabolomics 
# We will do these manually so it is exactly the same as metabolomics
proteomics_TUvsNG <- proteomics_TUvsNG[,c(6, 5, 12, 11, 16, 15, 14, 13, 2, 1, 27, 26, 24, 23, 18, 17, 9, 8, 4, 3)]
# "IJ17TV_NG"  "IJ17TV_TU" "8U04FR_NG"  "8U04FR_TU" these samples are from autoSP3 and not MTBE_SP3
colnames(proteomics_TUvsNG) <- gsub("[.]1", "", colnames(proteomics_TUvsNG))
colnames(proteomics_TUvsNG)[1] <- "8JDLQY_TU"

# save raw matrices
saveRDS(metabolomics_TUvsNG, "./results/metabolomics_TUvsNG.Rda")
saveRDS(proteomics_TUvsNG, "./results/proteomics_TUvsNG.Rda")

## computation for z-transformation
# log-FC

# metab
metabolomics_NG <- metabolomics_TUvsNG[,seq(2, 20, by = 2)]
metabolomics_TU <- metabolomics_TUvsNG[,seq(1, 20, by = 2)]
logFC_metabolomics_TUvsNG <- metabolomics_TU - metabolomics_NG
colnames(logFC_metabolomics_TUvsNG) <- gsub("_TU", "", colnames(logFC_metabolomics_TUvsNG))

# prot
proteomics_NG <- proteomics_TUvsNG[,seq(2, 20, by = 2)]
proteomics_TU <- proteomics_TUvsNG[,seq(1, 20, by = 2)]
logFC_proteomics_TUvsNG <- proteomics_TU - proteomics_NG
colnames(logFC_proteomics_TUvsNG) <- gsub("_TU", "", colnames(logFC_proteomics_TUvsNG))


## z-transformation
#metab
mean_metabolomics_FC <- rowMeans(as.matrix(logFC_metabolomics_TUvsNG), na.rm = T)
sd_metabolomics_FC <- rowSds(as.matrix(logFC_metabolomics_TUvsNG), na.rm = T)
logFC_metabolomics_z <- as.data.frame((logFC_metabolomics_TUvsNG - mean_metabolomics_FC)/sd_metabolomics_FC)

#prot
mean_proteomics_FC <- rowMeans(as.matrix(logFC_proteomics_TUvsNG), na.rm = T)
sd_proteomics_FC <- rowSds(as.matrix(logFC_proteomics_TUvsNG), na.rm = T)
logFC_proteomics_z <- as.data.frame((logFC_proteomics_TUvsNG - mean_proteomics_FC)/sd_proteomics_FC)

# save
saveRDS(logFC_metabolomics_z, "./results/logFC_metabolomics_z.Rda")
saveRDS(logFC_proteomics_z, "./results/logFC_proteomics_z.Rda")


