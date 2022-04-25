library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(SummarizedExperiment)

raw_metab <- readRDS("data/SummarizedExperiment_extraction_method_cohort_metabolomics.RDS")
raw_prot <- readRDS("data/SummarizedExperiment_extraction_method_cohort_proteomics.RDS")


# annotation tumor and healthy tissue
# raw_metab
metab_names <- as.list(raw_metab@colData@listData[["tissue"]])
raw_metab_count <- as.data.frame(raw_metab@assays@data@listData[[1]])
colnames(raw_metab_count) <- paste(colnames(raw_metab_count), metab_names, sep = "_")

# raw_prot
raw_prot_count <- as.data.frame(raw_prot@assays@data@listData[[1]])
