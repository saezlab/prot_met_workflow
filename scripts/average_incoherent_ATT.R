library(readr)
library(dplyr)
library(readxl)


ATT <- read.csv("results/final_run/ATT.csv")
ATT_back <- read.csv("results/final_run/ATT_back.csv")

ATT_full <- as.data.frame(rbind(ATT,ATT_back))
ATT_full <- unique(ATT_full)
ATT_full_which <- ATT_full
ATT_full_which[which(duplicated(ATT_full_which$Nodes)), 1]

ATT_full <- ATT_full %>%
  group_by(Nodes) %>%
  summarise_at(vars("ZeroAct", "UpAct", "DownAct", "AvgAct", "measured", "Activity", "t"), mean)


write_csv(ATT_full, file = paste("results/final_run/",paste("ATT_full_adjusted.csv",sep = ""), sep = ""))
