library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

## Data is already processed -> log2() transformation

raw_metab <- readRDS("data/SummarizedExperiment_extraction_method_cohort_metabolomics.RDS")
raw_prot <- readRDS("data/SummarizedExperiment_extraction_method_cohort_proteomics.RDS")


# annotation tumor and healthy tissue
# raw_metab
raw_metab_count <- as.data.frame(raw_metab@assays@data@listData[[1]])
colnames(raw_metab_count) <- raw_metab@colData@listData[["Sample.Description"]]
raw_metab_count <- tibble::rownames_to_column(raw_metab_count, "ID")

raw_metab_count <- unique(raw_metab_count)
raw_metabolomics_TUvsNG <- as.data.frame(raw_metab_count[,-1])
row.names(raw_metabolomics_TUvsNG) <- raw_metab_count$ID
metabolomics_TUvsNG <- raw_metabolomics_TUvsNG

# raw_prot
raw_prot_count <- as.data.frame(raw_prot@assays@data@listData[[1]])
colnames(raw_prot_count) <- raw_prot@colData@listData[["Pseudonym"]]
raw_prot_count <- tibble::rownames_to_column(raw_prot_count, "ID")


raw_prot_count <- raw_prot_count %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")
raw_prot_count <- raw_prot_count[raw_prot_count$ID != "NA",]
raw_prot_count <- unique(raw_prot_count)

raw_proteomics_TuvsNG <- as.data.frame(raw_prot_count[,-1])
row.names(raw_proteomics_TuvsNG) <- raw_prot_count$ID
# additionally differentiate between autoSP3 and MTBE
autoSP3orMTBE <- raw_prot@colData@listData[["Processing"]]
sum(autoSP3orMTBE == "MTBE_SP3")
# first 18 are MTBE_SP3, autoSP3 are only 17
# some samples do not have healthy or cancerous partner with sample preparation
# Thereby, these will be combined so it ends up with 20 samples that are pair-wise combined just like in metabolomics 
# We will do these manually so it is exactly the same as metabolomics
proteomics_TuvsNG <- raw_proteomics_TuvsNG[, c(6, 5, 12, 11, 16, 15, 14, 13, 2, 1, 27, 26, 24, 23, 18, 17, 9, 8, 4, 3)]
# "IJ17TV_NG"  "IJ17TV_TU" "8U04FR_NG"  "8U04FR_TU" these samples are from autoSP3 and not MTBE_SP3
colnames(proteomics_TuvsNG) <- gsub("[.]1", "", colnames(proteomics_TuvsNG))
colnames(proteomics_TuvsNG)[1] <- "8JDLQY_TU"


## computation for z-transformation
# log-FC

# metab
metabolomics_NG <- metabolomics_TUvsNG[,seq(2, 20, by = 2)]
metabolomics_TU <- metabolomics_TUvsNG[,seq(1, 20, by = 2)]
logFC_metabolomics_TUvsNG <- metabolomics_TU - metabolomics_NG
colnames(logFC_metabolomics_TUvsNG) <- gsub("_TU", "", colnames(logFC_metabolomics_TUvsNG))

# prot
proteomics_NG <- proteomics_TuvsNG[,seq(2, 20, by = 2)]
proteomics_TU <- proteomics_TuvsNG[,seq(1, 20, by = 2)]
logFC_proteomics_TUvsNG <- proteomics_TU - proteomics_NG
colnames(logFC_proteomics_TUvsNG) <- gsub("_TU", "", colnames(logFC_proteomics_TUvsNG))

## z-transformation
#metab
mean_metabolomics <- rowMeans(as.matrix(logFC_metabolomics_TUvsNG), na.rm = T)
sd_metabolomics <- rowSds(as.matrix(logFC_metabolomics_TUvsNG), na.rm = T)

metabolomics_z <- as.data.frame((logFC_metabolomics_TUvsNG - mean_metabolomics)/sd_metabolomics)

mean_proteomics <- rowMeans(as.matrix(logFC_proteomics_TUvsNG), na.rm = T)
sd_proteomics <- rowSds(as.matrix(logFC_proteomics_TUvsNG), na.rm = T)

proteomics_z <- as.data.frame((logFC_proteomics_TUvsNG - mean_proteomics)/sd_proteomics)

