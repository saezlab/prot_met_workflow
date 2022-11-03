library(readr)
library(readxl)
library(ocean)
library(visNetwork)
library(cosmosR)

data("meta_network")
meta_network <- meta_network[which(meta_network$source != meta_network$target),]

metabolomics_DE_t <- as.data.frame(read_delim("data/metabolomics_DE_t_TUvsNG.txt", 
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

t_table <- metabolomics_DE_t[,c(12,5)]
names(t_table) <- c("ID","MTBE_SP3")

t_table <- as.data.frame(t_table)

t_table <- t_table[t_table$ID %in% mapping_table$metab,]

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))

##Prepare the metabolic enzyme sets
penalty_min <- 7 #minimum 1 and integer
penalty_max <- 9 #maximum 9 and integer

######## SUBNETWORKS

# View(unique(recon2_redhuman$pathway))

all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)

sub_network <- translate_complexes(sub_network)

sub_network_nocofact <- remove_cofactors(sub_network)

sub_network_nocofact <- compress_transporters(sub_network_nocofact = sub_network_nocofact)

sub_network_nocofact <- split_transaminases(sub_network_nocofact = sub_network_nocofact)

enzymes <- unique(sub_network_nocofact$attributes$V1)
enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]

sub_forest <- forestMaker(enzymes, sub_network_nocofact$reaction_network, branch_length = c(1,1), remove_reverse = T)

###################

reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,   
                                            forest = sub_forest,
                                            measured_metabolites = t_table$KEGG)

reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)

penalty <- 8 #has be between penalty_min and penalty_max and integer

regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, filter_imbalance = c(0,1))

t_table <- t_table[!duplicated(t_table$KEGG),] #‘cpd:C00047’ is duplicated

##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 1000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES

mapping_table$metab <- gsub("histamineextracellular","histamine",mapping_table$metab)
mapping_table$metab <- gsub("sperminec10h30n4","spermine",mapping_table$metab)
translated_results <- translate_results(regulons_df = regulons_df, t_table = t_table, mapping_table = mapping_table)

translated_regulons_df <- translated_results$regulons_df
translated_regulons_df$ID <- paste(translated_regulons_df$set, gsub("_[a-z]$","",translated_regulons_df$targets), sep = "___")
translated_regulons_df <- translated_regulons_df[,-c(1,2)]

translated_regulons_df <- translated_regulons_df %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
translated_regulons_df <- as.data.frame(translated_regulons_df)

translated_regulons_df$set <- gsub("___.*","",translated_regulons_df$ID)
translated_regulons_df$targets <- gsub(".*___","",translated_regulons_df$ID)
translated_regulons_df <- translated_regulons_df[,c(3,4,2)]

##PROTEO
proteomics_DE_t_sp3 <- as.data.frame(read_delim("data/proteomics_DE_t_TUvsNG_autoSP3.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE))[,-c(1,2)]
names(proteomics_DE_t_sp3)[1] <- "ID"

proteomics_DE_t_sp3 <- proteomics_DE_t_sp3 %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")

proteomics_DE_t_MTBE <- as.data.frame(read_delim("data/proteomics_DE_t_TUvsNG_MTBE_SP3.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE))[,-c(1,2)]
names(proteomics_DE_t_MTBE)[1] <- "ID"

proteomics_DE_t_MTBE <- proteomics_DE_t_MTBE %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID) %>%
  filter(ID != "")

summarised_mean_NES_df <- mean_NES_df

#if only forward reactions are considered
# summarised_mean_NES_df <- summarised_mean_NES_df[!grepl("_reverse",summarised_mean_NES_df$KEGG),]

summarised_mean_NES_df$MTBE_SP3 <- ifelse(grepl("_reverse",summarised_mean_NES_df$KEGG),
                                                     summarised_mean_NES_df$MTBE_SP3 * -1,
                                                     summarised_mean_NES_df$MTBE_SP3)

summarised_mean_NES_df <- summarised_mean_NES_df[!grepl("_gluakg",summarised_mean_NES_df$KEGG),]
summarised_mean_NES_df <- summarised_mean_NES_df[!grepl("_glugln",summarised_mean_NES_df$KEGG),]

summarised_mean_NES_df$KEGG <- gsub("[>.].*","",summarised_mean_NES_df$KEGG)
summarised_mean_NES_df$KEGG <- gsub("_reverse","",summarised_mean_NES_df$KEGG)

summarised_mean_NES_df <- summarised_mean_NES_df %>% group_by(KEGG) %>% summarise_each(funs(mean(., na.rm = TRUE)))
summarised_mean_NES_df <- as.data.frame(summarised_mean_NES_df)

summarised_mean_NES_df <- summarised_mean_NES_df[,c(1,2)]

names(summarised_mean_NES_df) <- c("ID","ocean_NES")

summarised_mean_NES_df <- summarised_mean_NES_df %>%
  mutate(ID = strsplit(as.character(ID), "_")) %>%
  unnest(ID) %>%
  filter(ID != "")
summarised_mean_NES_df <- as.data.frame(summarised_mean_NES_df)
summarised_mean_NES_df <- merge(summarised_mean_NES_df, proteomics_DE_t_sp3[,c(1,4)])
summarised_mean_NES_df <- merge(summarised_mean_NES_df, proteomics_DE_t_MTBE[,c(1,4)], by = "ID")
names(summarised_mean_NES_df) <- c("ID","ocEAn","SP3","MTBE")

summarised_mean_NES_df <- summarised_mean_NES_df %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
summarised_mean_NES_df <- as.data.frame(summarised_mean_NES_df)

cor.test(summarised_mean_NES_df$ocEAn, summarised_mean_NES_df$SP3, method = "kendall")
cor.test(summarised_mean_NES_df$ocEAn, summarised_mean_NES_df$MTBE, method = "kendall")

plots <- plotMetaboliteContribution(enzyme = 'DLD_PDHX_PDHA1_PDHB_DLAT', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'SDHA_SDHD_SDHC_SDHB', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'PCK2', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'CS', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'IDH3B_IDH3G_IDH3A', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 1, vis.height = 2000) %>%
  visSave(file = "results/ocean_network.html")

sub_network_nocofact_tca <- sub_network_nocofact$reaction_network
tca_prots <- c("OGDH_DLD_PDH_DLST>1052",
               "SUCLG1_SUCLA2_reverse",
               "SDHA_SDHD_SDHC_SDHB",
               "FH>1055",
               "ME2",
               "DLD_PDHX_PDHA1_PDHB_DLAT",
               "CS",
               "ACO1",
               "IDH3B_IDH3G_IDH3A",
               "MDH2",
               "PKM>1130",
               "transporter>198",
               "CRAT")

sub_network_nocofact_tca <- sub_network_nocofact_tca[sub_network_nocofact_tca$source %in% tca_prots | sub_network_nocofact_tca$target %in% tca_prots,]
sub_network_nocofact_tca_att <- sub_network_nocofact$attributes
sub_network_nocofact_tca_att <- sub_network_nocofact_tca_att[sub_network_nocofact_tca_att$V1 %in% unique(c(sub_network_nocofact_tca$source,sub_network_nocofact_tca$target)),]

sub_network_nocofact_tca <- list(sub_network_nocofact_tca, sub_network_nocofact_tca_att)

plot_reaction_network(sub_network_nocofact_tca, t_table, mean_NES_df, column_index = 1, vis.height = 2000) %>%
  visSave(file = "results/ocean_network_tca.html")

tca_prots_clean <- gsub(">.*","",tca_prots)
tca_prots_clean <- gsub("_reverse*","",tca_prots_clean)
tca_prots_clean <- unlist(do.call(strsplit,list(tca_prots_clean,"_")))

####################
summarised_mean_NES_df_tca <- mean_NES_df[mean_NES_df$KEGG %in% tca_prots,]

#if only forward reactions are considered
# summarised_mean_NES_df <- summarised_mean_NES_df[!grepl("_reverse",summarised_mean_NES_df$KEGG),]

summarised_mean_NES_df_tca$MTBE_SP3 <- ifelse(grepl("_reverse",summarised_mean_NES_df_tca$KEGG),
                                              summarised_mean_NES_df_tca$MTBE_SP3 * -1,
                                              summarised_mean_NES_df_tca$MTBE_SP3)

summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca[!grepl("_gluakg",summarised_mean_NES_df_tca$KEGG),]
summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca[!grepl("_glugln",summarised_mean_NES_df_tca$KEGG),]

summarised_mean_NES_df_tca$KEGG <- gsub("[>.].*","",summarised_mean_NES_df_tca$KEGG)
summarised_mean_NES_df_tca$KEGG <- gsub("_reverse","",summarised_mean_NES_df_tca$KEGG)

summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca %>% group_by(KEGG) %>% summarise_each(funs(mean(., na.rm = TRUE)))
summarised_mean_NES_df_tca <- as.data.frame(summarised_mean_NES_df_tca)

summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca[,c(1,2)]

names(summarised_mean_NES_df_tca) <- c("ID","ocean_NES")

summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca %>%
  mutate(ID = strsplit(as.character(ID), "_")) %>%
  unnest(ID) %>%
  filter(ID != "")
summarised_mean_NES_df_tca <- as.data.frame(summarised_mean_NES_df_tca)
summarised_mean_NES_df_tca <- merge(summarised_mean_NES_df_tca, proteomics_DE_t_sp3[,c(1,4)])
summarised_mean_NES_df_tca <- merge(summarised_mean_NES_df_tca, proteomics_DE_t_MTBE[,c(1,4)], by = "ID")
names(summarised_mean_NES_df_tca) <- c("ID","ocEAn_score","SP3_tvalue","MTBE_tvalue")
summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca[complete.cases(summarised_mean_NES_df_tca),]

summarised_mean_NES_df_tca <- summarised_mean_NES_df_tca %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
summarised_mean_NES_df_tca <- as.data.frame(summarised_mean_NES_df_tca)

to_hm <- summarised_mean_NES_df_tca
to_hm <- to_hm[!duplicated(to_hm[,1]),]
row.names(to_hm) <- to_hm[,1]
to_hm <- to_hm[,-1]

pheatmap::pheatmap(to_hm[,c(3,1,2)], display_numbers = T, cluster_cols = F)

cor.test(to_hm$ocEAn_score, to_hm$SP3_tvalue, method = "pearson")
cor.test(to_hm$ocEAn_score, to_hm$MTBE_tvalue, method = "pearson")

plots <- plotMetaboliteContribution(enzyme = 'SDHA_SDHD_SDHC_SDHB', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'PKM', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)

plots <- plotMetaboliteContribution(enzyme = 'MDH2', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 3, nLabels =  15)

plot(plots$scatter)
