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
PEA_metab <- metabolomics_TUvsNG[, c(3:8,15,16)]

# logFC metab
PEA_metab_NG <- PEA_metab[,seq(2, 8, by = 2)]
PEA_metab_TU <- PEA_metab[,seq(1, 8, by = 2)]
logFC_PEA_metab <- PEA_metab_TU - PEA_metab_NG
colnames(logFC_PEA_metab) <- gsub("_TU", "", colnames(logFC_PEA_metab))

# z-score metab
mean_PEA_metab_FC <- rowMeans(as.matrix(logFC_PEA_metab), na.rm = T)
sd_PEA_metab_FC <- rowSds(as.matrix(logFC_PEA_metab), na.rm = T)
logFC_PEA_metab_z <- as.data.frame((logFC_PEA_metab - mean_PEA_metab_FC)/sd_PEA_metab_FC)

saveRDS(logFC_PEA_metab_z, "./results/single_patient/logFC_PEA_metab_z.Rda")


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
# logFC prot1
PEA_prot_MTBE_NG <- PEA_prot_MTBE[,seq(2, 8, by = 2)]
PEA_prot_MTBE_TU <- PEA_prot_MTBE[,seq(1, 8, by = 2)]
logFC_PEA_prot_MTBE <- PEA_prot_MTBE_TU - PEA_prot_MTBE_NG
colnames(logFC_PEA_prot_MTBE) <- gsub("_TU", "", colnames(logFC_PEA_prot_MTBE))

# z-score prot1
mean_PEA_prot_MTBE_FC <- rowMeans(as.matrix(logFC_PEA_prot_MTBE), na.rm = T)
sd_PEA_prot_MTBE_FC <- rowSds(as.matrix(logFC_PEA_prot_MTBE), na.rm = T)
logFC_PEA_prot_MTBE_z <- as.data.frame((logFC_PEA_prot_MTBE - mean_PEA_prot_MTBE_FC)/sd_PEA_prot_MTBE_FC)

saveRDS(logFC_PEA_prot_MTBE_z, "./results/single_patient/logFC_PEA_prot_MTBE_z.Rda")


# Prot 2
PEA_prot_autoSP3 <- prot_count[, c(34, 33, 30, 29, 35, 28, 32, 31)]
colnames(PEA_prot_autoSP3) <- gsub("[.]1", "", colnames(PEA_prot_autoSP3))
# logFC prot2
PEA_prot_autoSP3_NG <- PEA_prot_autoSP3[,seq(2, 8, by = 2)]
PEA_prot_autoSP3_TU <- PEA_prot_autoSP3[,seq(1, 8, by = 2)]
logFC_PEA_prot_autoSP3 <- PEA_prot_autoSP3_TU - PEA_prot_autoSP3_NG
colnames(logFC_PEA_prot_autoSP3) <- gsub("_TU", "", colnames(logFC_PEA_prot_autoSP3))

# z-score prot2
mean_PEA_prot_autoSP3_FC <- rowMeans(as.matrix(logFC_PEA_prot_autoSP3), na.rm = T)
sd_PEA_prot_autoSP3_FC <- rowSds(as.matrix(logFC_PEA_prot_autoSP3), na.rm = T)
logFC_PEA_prot_autoSP3_z <- as.data.frame((logFC_PEA_prot_autoSP3 - mean_PEA_prot_autoSP3_FC)/sd_PEA_prot_autoSP3_FC)

saveRDS(logFC_PEA_prot_autoSP3_z, "./results/single_patient/logFC_PEA_prot_autoSP3_z.Rda")



