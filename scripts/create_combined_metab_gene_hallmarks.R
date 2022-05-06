library(readr)
library(readxl)
library(cosmosR)
library(GSEABase)

source("scripts/support_functions.R")

pathways <- import_gmt("support/h.all.v7.5.1.symbols.gmt")

data("meta_network")
meta_network <- meta_network[which(meta_network$source != meta_network$target),]

metabs <- meta_network[grepl("Metab__HMDB", meta_network$source) | grepl("Metab__HMDB", meta_network$target),]
metabs <- metabs[,-2]
metabs <- metabs[grepl("Gene", metabs$source) | grepl("Gene", metabs$target),]

metabs_reactant <- metabs[grepl("Metab__HMDB",metabs$source),]
metabs_products <- metabs[grepl("Metab__HMDB",metabs$target),]

names(metabs_reactant) <- c("metab","gene") 
names(metabs_products) <- c("gene","metab")

metabs <- as.data.frame(rbind(metabs_reactant, metabs_products))
metabs$gene <- gsub("Gene.*__","",metabs$gene)
metabs$metab <- gsub("_[a-z]$","",metabs$metab)
metabs$metab <- gsub("Metab__","",metabs$metab)

metabs <- merge(metabs,pathways)
metabs <- metabs[,-1]
names(metabs) <- c("feature","term")
names(pathways) <- c("feature","term")

combined_metab_gene_hallmarks <- as.data.frame(rbind(pathways,metabs))
combined_metab_gene_hallmarks <- unique(combined_metab_gene_hallmarks)

write_csv(combined_metab_gene_hallmarks, file = "support/combined_metab_gene_hallmarks.csv")
