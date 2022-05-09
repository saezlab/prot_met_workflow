library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)
library(dplyr)
library(tidyr)

data("meta_network")
meta_network <- meta_network[which(meta_network$source != meta_network$target),]

# use single patients of logFC_metabolomics_z and logFC_proteomics_z
# 3 patients: e.g. "1FF2F9", "36AT2O" and "K6R512"

patient <-"K6R512"

# input is z-score of logFC

# single patient input
# metab
met_cosmos <- as.data.frame(matrix(NA, ncol = 2, nrow = dim(logFC_metabolomics_z)[1]))
colnames(met_cosmos) <- c("ID", "logFC_z")
met_cosmos$ID <- rownames(logFC_metabolomics_z)
met_cosmos$logFC_z<- logFC_metabolomics_z[,which(colnames(logFC_metabolomics_z) == patient)]

# filter top 50 metab
met_cosmos <- met_cosmos %>%
  arrange(desc(abs(logFC_z))) %>%
  slice(1:50)


# prot
prot_cosmos <- as.data.frame(matrix(NA, ncol = 2, nrow = dim(logFC_proteomics_z)[1]))
colnames(prot_cosmos) <- c("ID", "logFC_z")
prot_cosmos$ID <- rownames(logFC_proteomics_z)
prot_cosmos$logFC_z<- logFC_proteomics_z[,which(colnames(logFC_proteomics_z) == patient)]

# filter top 100 prot
prot_cosmos <- prot_cosmos %>%
  arrange(desc(abs(logFC_z))) %>%
  slice(1:100)


Biocrates_metabolite_identifier <- read_excel("support/Biocrates_metabolite_identifier.xlsx")
Biocrates_metabolite_identifier$feature <- gsub("[()/: -]",".",Biocrates_metabolite_identifier$feature)

met_cosmos <- merge(met_cosmos, Biocrates_metabolite_identifier, by.x = "ID", by.y = "feature")
# why less -> x1.Met.His and X3.Met.His left out for 1FF2F9. Is this a problem?

met_cosmos <- met_cosmos %>%
  mutate(HMDB = strsplit(as.character(HMDB), "/")) %>%
  unnest(HMDB) %>%
  filter(HMDB != "")

mapping <- met_cosmos[,c(1,6)]

cosmos_met_input <- unlist(met_cosmos$logFC_z)
names(cosmos_met_input) <- met_cosmos$HMDB

cosmos_met_input <- prepare_metab_inputs(cosmos_met_input, c("c","m"))

cosmos_prot_input <- unlist(prot_cosmos$logFC_z)
names(cosmos_prot_input) <- prot_cosmos$ID

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cplex")

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
my_options$solverPath <- "cplex/cplex.exe"
# my_options$solverPath <- "cbc/cbc-osx/cbc" #or cbc solver executable
my_options$solver <- "cplex" #or cbc
# my_options$solver <- "cbc"
my_options$timelimit <- 1800
my_options$mipGAP <- 0.05
my_options$threads <- 6

metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(cosmos_met_input, meta_network)
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(cosmos_prot_input, meta_network)

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network,
                                                      signaling_data = sig_input,
                                                      metabolic_data = metab_input,
                                                      maximum_network_depth = 4,
                                                      CARNIVAL_options = my_options)


my_options$timelimit <- 600

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

SIF <- SIF[which(SIF$Weight != 0),]

ATT <- merge(ATT, proteomics_DE_t[,c(3,6)], all.x = T, by.x = "Nodes", by.y = "ID")
ATT$Nodes <- gsub(",","_",ATT$Nodes)

write_csv(SIF, file = paste("results/single_patient/",paste(patient, "_SIF.csv",sep = ""), sep = ""))
write_csv(ATT, file = paste("results/single_patient/",paste(patient, "_ATT.csv",sep = ""), sep = ""))

# include backwards run

my_options$timelimit <- 1800

test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network,
                                                       signaling_data = sig_input,
                                                       metabolic_data = metab_input,
                                                       maximum_network_depth = 4,
                                                       CARNIVAL_options = my_options)

my_options$timelimit <- 600

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)

formatted_res_back <- format_COSMOS_res(test_result_back)

SIF_back <- formatted_res_back[[1]]
ATT_back <- formatted_res_back[[2]]

SIF_back <- SIF_back[which(SIF_back$Weight != 0),]


ATT_back <- merge(ATT_back, proteomics_DE_t[,c(3,6)], all.x = T, by.x = "Nodes", by.y = "ID")
ATT_back$Nodes <- gsub(",","_",ATT_back$Nodes)


write_csv(SIF_back, file = paste("results/single_patient/",paste(patient, "_SIF_back.csv",sep = ""), sep = ""))
write_csv(ATT_back, file = paste("results/single_patient/",paste(patient, "_ATT_back.csv",sep = ""), sep = ""))

SIF_full <- as.data.frame(rbind(SIF,SIF_back))
SIF_full <- unique(SIF_full)

ATT_full <- as.data.frame(rbind(ATT,ATT_back))
ATT_full <- unique(ATT_full)


ATT_full <- as.data.frame(ATT_full)

write_csv(SIF_full, file = paste("results/single_patient/",paste(patient, "_SIF_full.csv",sep = ""), sep = ""))
write_csv(ATT_full, file = paste("results/single_patient/",paste(patient, "_ATT_full.csv",sep = ""), sep = ""))

