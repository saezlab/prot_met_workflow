library(readr)
library(readxl)
library(cosmosR)
library(decoupleR)

data("meta_network")
meta_network <- meta_network[which(meta_network$source != meta_network$target),]

metabolomics_DE_t <- as.data.frame(read_delim("data/metabolomics_DE_t_TUvsNG.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE))
proteomics_DE_t <- as.data.frame(read_delim("data/proteomics_DE_t_TUvsNG_autoSP3.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE))

Biocrates_metabolite_identifier <- read_excel("support/Biocrates_metabolite_identifier.xlsx")
Biocrates_metabolite_identifier$feature <- gsub("[()/: -]",".",Biocrates_metabolite_identifier$feature)

metabolomics_DE_t <- merge(metabolomics_DE_t, Biocrates_metabolite_identifier, by.x = "ID", by.y = "feature")

library(dplyr)
library(tidyr)

metabolomics_DE_t <- metabolomics_DE_t %>%
  mutate(HMDB = strsplit(as.character(HMDB), "/")) %>%
  unnest(HMDB) %>%
  filter(HMDB != "")

mapping <- metabolomics_DE_t[,c(1,12)]

cosmos_met_input <- metabolomics_DE_t$t
names(cosmos_met_input) <- metabolomics_DE_t$HMDB

proteomics_DE_t <- proteomics_DE_t %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")

cosmos_prot_input <- proteomics_DE_t$t
names(cosmos_prot_input) <- proteomics_DE_t$ID


cosmos_met_input <- cosmos_met_input[which(abs(cosmos_met_input) > 3)]
cosmos_prot_input <- cosmos_prot_input[which(abs(cosmos_prot_input) > 6)]

cosmos_met_input <- prepare_metab_inputs(cosmos_met_input, c("c","m"))

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
                                                      maximum_network_depth = 5,
                                                      CARNIVAL_options = my_options)

my_options$timelimit <- 1200

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

SIF <- SIF[which(SIF$Weight != 0),]

ATT <- merge(ATT, proteomics_DE_t[,c(3,6)], all.x = T, by.x = "Nodes", by.y = "ID")
ATT$Nodes <- gsub(",","_",ATT$Nodes)

write_csv(SIF, file = paste("results/final_run/", "SIF_n6.csv", sep = ""))
write_csv(ATT, file = paste("results/final_run/", "ATT_n6.csv", sep = ""))


# include backwards run

my_options$timelimit <- 1800


test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network,
                                                       signaling_data = sig_input,
                                                       metabolic_data = metab_input,
                                                       maximum_network_depth = 4,
                                                       CARNIVAL_options = my_options)

my_options$timelimit <- 1800

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)

formatted_res_back <- format_COSMOS_res(test_result_back)

SIF_back <- formatted_res_back[[1]]
ATT_back <- formatted_res_back[[2]]

SIF_back <- SIF_back[which(SIF_back$Weight != 0),]


ATT_back <- merge(ATT_back, proteomics_DE_t[,c(3,6)], all.x = T, by.x = "Nodes", by.y = "ID")
ATT_back$Nodes <- gsub(",","_",ATT_back$Nodes)


write_csv(SIF_back, file = paste("results/final_run/",paste("SIF_back.csv",sep = ""), sep = ""))
write_csv(ATT_back, file = paste("results/final_run/",paste("ATT_back.csv",sep = ""), sep = ""))

SIF_full <- as.data.frame(rbind(SIF,SIF_back))
SIF_full <- unique(SIF_full)

ATT_full <- as.data.frame(rbind(ATT,ATT_back))
ATT_full <- unique(ATT_full)


ATT_full <- as.data.frame(ATT_full)

write_csv(SIF_full, file = paste("results/final_run/",paste("SIF_full.csv",sep = ""), sep = ""))
write_csv(ATT_full, file = paste("results/final_run/",paste("ATT_full.csv",sep = ""), sep = ""))


