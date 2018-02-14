library(data.table)
library(plyr)
library(xlsx)

bc_final <- read.csv("breast_cancer.csv", sep = ",", header = T, stringsAsFactors = F)
bc_remain <- read.csv("remain.csv", sep = ",", header = T, stringsAsFactors = F)
bc_final <- bc_final[bc_final$GPL %in% bc_remain$GPL, ]
GSE <- unique(bc_final$GSE)
write.table(x = GSE, file = "GSE.txt", sep = "\n", row.names = F, col.names = F)

target_bc <- fread("breast_cancer.csv", sep = ",", stringsAsFactors = F, header = T)
target_bc <- target_bc[, colnames(target_bc)[21:50] := NULL]
target_bc_pre <- read.csv("breast_cancer_2016_12.csv", header = T, stringsAsFactors = F)
bc_all <- union(unique(target_bc_pre$GPL), bc_remain$GPL)
bc_add <- setdiff(bc_all, unique(target_bc_pre$GPL))
target_bc <- target_bc[target_bc$GPL %in% bc_add, ]
target_bc[is.na(target_bc)] <- ""
target_bc <- as.data.frame(target_bc)
write.csv(x = target_bc, file = "bc_add.csv", sep = ",", row.names = F)

i <- 0
for (summary in unique(target_bc$GSE_Summary)) {
  subbc <- subset(target_bc, target_bc$GSE_Summary == summary)
  label <- c()
  label[1] <- ifelse(any(subbc$Ch1_characteristic != ""), 1, 0)
  label[2] <- ifelse(any(subbc$Ch2_characteristic != ""), 1, 0)
  i <- i + 1
  if (label[1] == 0 && label[2] == 0) {next}
  reshapeForBCData(data = subbc, ch1 = label[1], ch2 = label[2])
  write.csv(x = cbind(subbc, chframe), file = paste(i, ".csv", sep = ""), sep = ",", row.names = F, col.names = T)
}

colnamebc <- read.table("../数据收集模板.csv", sep = ",", stringsAsFactors = F)[, 1]
target_bc <- data.table(target_bc)
target_bc <- target_bc[, colnamebc[21:54] := NA]
target_bc <- as.data.frame(target_bc)

library(tidyr)
bc_cellLine <- FillBCData(target_data = target_bc, target_col = "cell_line", target_col_position = 21)
target_bc[, 21] <- bc_cellLine[, 21]
bc_cellType <- FillBCData(target_data = target_bc, target_col = "cell_type", target_col_position = 22)
target_bc[, 22] <- bc_cellType[, 22]
bc_subtype <- FillBCData(target_data = target_bc, target_col = "tumor_subtype", target_col_position = 23)
target_bc[, 23] <- bc_subtype[, 23]
bc_age <- FillBCData(target_data = target_bc, target_col = "^age", target_col_position = 24)
target_bc[, 24] <- bc_age[, 24]
bc_gender <- FillBCData(target_data = target_bc, target_col = "gender", target_col_position = 25)
target_bc[, 25] <- bc_gender[, 25]
bc_race1 <- FillBCData(target_data = target_bc, target_col = "race", target_col_position = 26)
bc_race2 <- FillBCData(target_data = target_bc, target_col = "ethnicity", target_col_position = 26)
bc_race <- cbind(race1 = bc_race1[, 26, drop = F], race2 = bc_race2[, 26, drop = F])
bc_race <- unite(bc_race, col = "race_total", race:race, sep = " ")
target_bc[, 26] <- bc_race
bc_stage <- FillBCData(target_data = target_bc, target_col = "stage", target_col_position = 27)
target_bc[, 27] <- bc_stage[, 27]
bc_tstage1 <- FillBCData(target_data = target_bc, target_col = "t_stage", target_col_position = 28)
bc_tstage2 <- FillBCData(target_data = target_bc, target_col = "tumor_stage", target_col_position = 28)
bc_tstage3 <- FillBCData(target_data = target_bc, target_col = "stage_t", target_col_position = 28)
bc_tstage <- cbind(bc_tstage1[, 28, drop = F], bc_tstage2[, 28, drop = F], bc_tstage3[, 28, drop = F])
bc_tstage <- unite(bc_tstage, col = "tstage_total", tStage:tStage, sep = " ")
target_bc[, 28] <- bc_tstage
bc_nstage1 <- FillBCData(target_data = target_bc, target_col = "n_stage", target_col_position = 29)
bc_nstage2 <- FillBCData(target_data = target_bc, target_col = "node_stage", target_col_position = 29)
bc_nstage3 <- FillBCData(target_data = target_bc, target_col = "stage_n", target_col_position = 29)
bc_nstage <- cbind(bc_nstage1[, 29, drop = F], bc_nstage2[, 29, drop = F], bc_nstage3[, 29, drop = F])
bc_nstage <- unite(bc_nstage, col = "nstage_total", nStage:nStage, sep = " ")
target_bc[, 29] <- bc_nstage
bc_mstage1 <- FillBCData(target_data = target_bc, target_col = "m_stage", target_col_position = 30)
bc_mstage2 <- FillBCData(target_data = target_bc, target_col = "metastasis_stage", target_col_position = 30)
bc_mstage3 <- FillBCData(target_data = target_bc, target_col = "stage_m", target_col_position = 30)
bc_mstage <- cbind(bc_mstage1[, 30, drop = F], bc_mstage2[, 30, drop = F], bc_mstage3[, 30, drop = F])
bc_mstage <- unite(bc_mstage, col = "mstage_total", mStage:mStage, sep = " ")
target_bc[, 30] <- bc_mstage
bc_tgrade <- FillBCData(target_data = target_bc, target_col = "grade", target_col_position = 31)
target_bc[, 31] <- bc_tgrade[, 31]
bc_site <- FillBCData(target_data = target_bc, target_col = "site", target_col_position = 32)
target_bc[, 32] <- bc_site[, 32]
bc_vital <- FillBCData(target_data = target_bc, target_col = "vital", target_col_position = 33)# no record
target_bc[, 33] <- bc_vital[, 33]
bc_osTime <- FillBCData(target_data = target_bc, target_col = "survival", target_col_position = 34)
target_bc[, 34] <- bc_osTime[, 34]
bc_rfs_event <- FillBCData(target_data = target_bc, target_col = "rfs", target_col_position = 35)
target_bc[, 35] <- bc_rfs_event[, 35]
bc_rfs_time <- FillBCData(target_data = target_bc, target_col = "rfs", target_col_position = 36)
target_bc[, 36] <- bc_rfs_time[, 36]
bc_pfs_event <- FillBCData(target_data = target_bc, target_col = "pfs", target_col_position = 37)# no record
bc_pfs_time <- FillBCData(target_data = target_bc, target_col = "pfs", target_col_position = 38)# no record
bc_dfs_event <- FillBCData(target_data = target_bc, target_col = "dfs", target_col_position = 39)
target_bc[, 39] <- bc_dfs_event[, 39]
bc_dfs_time <- FillBCData(target_data = target_bc, target_col = "dfs", target_col_position = 40)
target_bc[, 40] <- bc_dfs_time[, 40]
bc_dss_event <- FillBCData(target_data = target_bc, target_col = "dss", target_col_position = 41)
target_bc[, 41] <- bc_dss_event[, 41]
bc_dss_time <- FillBCData(target_data = target_bc, target_col = "dss", target_col_position = 42)
target_bc[, 42] <- bc_dss_time[, 42]
bc_treatment <- FillBCData(target_data = target_bc, target_col = "treat", target_col_position = 43)
target_bc[, 43] <- bc_treatment[, 43]
bc_disease <- FillBCData(target_data = target_bc, target_col = "disease", target_col_position = 44)
target_bc[, 44] <- bc_disease[, 44]
bc_ToN <- FillBCData(target_data = target_bc, target_col = "disease_state", target_col_position = 45)# DN pattern in data
target_bc[, 45] <- bc_ToN[, 45]
bc_msi <- FillBCData(target_data = target_bc, target_col = "msi", target_col_position = 46)# no record
bc_cin <- FillBCData(target_data = target_bc, target_col = "neoplasia", target_col_position = 47)# no record #try cin, carvical, intraepithelial, neoplasia as pattern
bc_prs <- FillBCData(target_data = target_bc, target_col = "progesteron", target_col_position = 48)
target_bc[, 48] <- bc_prs[, 48]
bc_ers <- FillBCData(target_data = target_bc, target_col = "estrogen", target_col_position = 49)
target_bc[, 49] <- bc_ers[, 49]
bc_her2 <- FillBCData(target_data = target_bc, target_col = "her2", target_col_position = 50)
target_bc[, 50] <- bc_her2[, 50]
bc_pr <- FillBCData(target_data = target_bc, target_col = "pr", target_col_position = 51, fix = T)
target_bc[, 51] <- bc_pr[, 51]
bc_er <- FillBCData(target_data = target_bc, target_col = "er", target_col_position = 52, fix = T)
target_bc[, 52] <- bc_er[, 52]
bc_her <- FillBCData(target_data = target_bc, target_col = "her", target_col_position = 53, fix = T)
target_bc[, 53] <- bc_her[, 53]
bc_ms <- FillBCData(target_data = target_bc, target_col = "menopausal", target_col_position = 54, fix = T)
target_bc[, 54] <- bc_ms[, 54]


colnames(target_bc) <- colnames(target_bc_pre)
bc_final <- rbind(target_bc, target_bc_pre)
bc_final <- bc_final[bc_final$GPL %in% bc_remain$GPL, ]
bc_final[is.na(bc_final)] <- ""

write.csv(bc_final, file = "breast_cancer_20161213.csv", row.names = F)






#[584:length(unique(target_bc$GSE_Summary))]
