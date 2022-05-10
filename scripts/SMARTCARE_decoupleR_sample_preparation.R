library(readr)
library(readxl)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)

raw_metab <- readRDS("data/SummarizedExperiment_extraction_method_cohort_metabolomics.RDS")
raw_prot <- readRDS("data/SummarizedExperiment_extraction_method_cohort_proteomics.RDS")

# raw_metab
metab_count <- as.data.frame(raw_metab@assays@data@listData[[1]])
colnames(metab_count) <- raw_metab@colData@listData[["Sample.Description"]]

metab_count <- unique(metab_count)
metabolomics_TUvsNG <- as.data.frame(metab_count)

# For pathway enrichment analysis
PEA_metab <- metabolomics_TUvsNG[, c(3:8,15 ,16)]

# raw_prot
# MTBESP3 = Prot1
# autoSP3 = Prot2
prot_count <- as.data.frame(raw_prot@assays@data@listData[[1]])
colnames(prot_count) <- raw_prot@colData@listData[["Pseudonym"]]
prot_count <- tibble::rownames_to_column(prot_count, "ID")

prot_count <- prot_count %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")

prot_count <- prot_count[prot_count$ID != "NA",]
prot_count <- unique(prot_count)
prot_count <- prot_count %>%
  column_to_rownames("ID")

# Prot1
PEA_prot_MTBE <- prot_count[, c(12, 11, 16, 15, 14, 13, 18, 17)]
# Prot 2
PEA_prot_autoSP3 <- prot_count[, c(34, 33, 30, 29, 35, 28, 32, 31)]
