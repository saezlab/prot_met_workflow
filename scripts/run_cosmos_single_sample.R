library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(dplyr)
library(tidyr)

data("meta_network")
meta_network <- meta_network[which(meta_network$source != meta_network$target),]

# use single patients of logFC_metabolomics_z and logFC_proteomics_z
# 3 patients: e.g. "1FF2F9", "36AT2O" and "C5FQXS"

patient <-"1FF2F9"

# input is z-score of logFC

# single patient input
# metab
met_cosmos <- as.data.frame(matrix(NA, ncol = 2, nrow = dim(logFC_metabolomics_z)[1]))
colnames(met_cosmos) <- c("ID", "logFC_z")
met_cosmos$ID <- logFC_metabolomics_z[,1]
met_cosmos$logFC_z<- logFC_metabolomics_z[,which(colnames(logFC_metabolomics_z) == patient)]

# filter top 50 metab
met_cosmos <- met_cosmos %>%
  arrange(desc(abs(logFC_z))) %>%
  slice(1:50)


# prot
prot_cosmos <- as.data.frame(matrix(NA, ncol = 2, nrow = dim(logFC_proteomics_TUvsNG)[1]))
colnames(prot_cosmos) <- c("ID", "logFC_z")
prot_cosmos$ID <- logFC_proteomics_TUvsNG[,1]
prot_cosmos$logFC_z<- logFC_proteomics_TUvsNG[,which(colnames(logFC_proteomics_TUvsNG) == patient)]

# filter top 100 prot
prot_cosmos <- prot_cosmos %>%
  arrange(desc(abs(logFC_z))) %>%
  slice(1:100)


Biocrates_metabolite_identifier <- read_excel("support/Biocrates_metabolite_identifier.xlsx")
Biocrates_metabolite_identifier$feature <- gsub("[()/: -]",".",Biocrates_metabolite_identifier$feature)

met_cosmos1 <- merge(met_cosmos, Biocrates_metabolite_identifier, by.x = "ID", by.y = "feature")
# why less -> x1.Met.His and X3.Met.His left out





met_cosmos <- met_cosmos %>%
  mutate(HMDB = strsplit(as.character(HMDB), "/")) %>%
  unnest(HMDB) %>%
  filter(HMDB != "") 





mapping <- met_cosmos[,c(1,6)]

cosmos_met_input <- unlist(met_cosmos[,which(colnames(met_cosmos) == patient)])
names(cosmos_met_input) <- met_cosmos$HMDB


prot_cosmos <- prot_cosmos %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")




write_csv(SIF, file = paste("results/",paste(patient, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/",paste(patient, "_ATT.csv",sep = ""), sep = ""))
