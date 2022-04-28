library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)



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

#raw_metabolomics_TUvsNG <- as.data.frame(raw_metabolomics_TUvsNG) %>%
#            dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x))
# maybe after log-Transformation

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

#raw_proteomics_TuvsNG <- as.data.frame(raw_proteomics_TuvsNG) %>%
#  dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x))
# replace NA values in dataframe with 0 or remove rows?


