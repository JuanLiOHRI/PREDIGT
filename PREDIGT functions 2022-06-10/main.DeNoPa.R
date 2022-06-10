# Main funtion for PREDIGT analysis on DeNoPa
# Created by Juan Li, 2022-03-21

library(dplyr)
library(pROC)
library(ggplot2)
library(ggpubr)
library(rms)
library(pracma)
library(boot)

# source functions
wd <- "/Users/juan/Desktop/OHRI/work/Cohorts/! PREDIGT functions"
source(paste(wd, '/plotBar.R', sep=""), echo=TRUE)
source(paste(wd, '/factor.transfer.R', sep=""), echo=TRUE)
source(paste(wd, '/data.impute.R', sep=""), echo=TRUE)
source(paste(wd, '/varScore.R', sep=""), echo=TRUE)
source(paste(wd, '/uniVar.R', sep=""), echo=TRUE)
source(paste(wd, '/PREDIGTScr.R', sep=""), echo=TRUE)
source(paste(wd, '/plotROC.R', sep=""), echo=TRUE)
source(paste(wd, '/plotDistBox.R', sep=""), echo=TRUE)
source(paste(wd, '/model2.R', sep=""), echo=TRUE)
source(paste(wd, '/model2point.R', sep=""), echo=TRUE)
source(paste(wd, '/model2ROC.R', sep=""), echo=TRUE)
source(paste(wd, '/model2box.R', sep=""), echo=TRUE)
source(paste(wd, '/model2boot.R', sep=""), echo=TRUE)

# set up working path
data_path  <- "/Users/juan/Desktop/OHRI/work/Cohorts/!DeNoPa PPMI Temporal/DeNoPa"
setwd(data_path)
plot_path  <- paste(data_path, "/Plots", sep="")
group_name <- c("Healthy Control", "Parkinson's Disease")
# ------------- Read data -----------
# read in all data
data_BL  <- read.csv(paste(data_path, "/", "DeNoPa_Baseline_curated.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
data_BL2 <- read.csv(paste(data_path, "/", "DeNoPa data cleaning.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
data_24  <- read.csv(paste(data_path, "/", "DeNoPa_24months_curated.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
data_48  <- read.csv(paste(data_path, "/", "DeNoPa_48months_curated.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
data_V   <- read.csv(paste(data_path, "/", "variable list.csv",sep=""), header = F, stringsAsFactors=FALSE, fileEncoding="latin1", na.strings=c("",".","NA"))
vfct <- data_V[which(data_V[,1]==1),2] # names of factor variables

genetic <- read.csv(paste(data_path, "/", "DeNoPa Genetic 2022-04-22.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))

# ------------ data prepatation -------------------
# select variables that will be used in PREDIGT model
df_BL1 <- data_BL %>% select(names(data_BL[which(names(data_BL) %in% data_V[,3])]))
df_BL2 <- data_BL2 %>% select(names(data_BL2[which(names(data_BL2) %in% data_V[,3])]))
df_BL  <- df_BL2 %>% inner_join(df_BL1, by="Patient ID")
df_BL  <- df_BL[,data_V[which(data_V[,3] %in% names(df_BL)),3]] # reorder by column name
names(df_BL) <- data_V[which(data_V[,3] %in% names(df_BL)),2]# update variable names
df_BL  <- df_BL %>% mutate(subID = as.character(subID))

df_24  <- data_24 %>% select(names(data_24[which(names(data_24) %in% data_V[,4])]))
df_24  <- df_24[,data_V[which(data_V[,4] %in% names(df_24)),4]] # reorder by column name
names(df_24) <- data_V[which(data_V[,4] %in% names(df_24)),2]# update variable names
df_24  <- df_24 %>% mutate(subID = as.character(subID))

df_48  <- data_48 %>% select(names(data_48[which(names(data_48) %in% data_V[,5])]))
df_48  <- df_48[,data_V[which(data_V[,5] %in% names(df_48)),5]] # reorder by column name
names(df_48) <- data_V[which(data_V[,5] %in% names(df_48)),2]# update variable names
df_48  <- df_48 %>% mutate(subID = as.character(subID))

# Metal - data only available in baseline
df_BL <- df_BL %>% 
  mutate(Metal.work1 = ifelse(is.na(Metal.work1),0,ifelse(Metal.work1==1, 1, 0))) %>% 
  mutate(Metal.work2 = ifelse(is.na(Metal.work2),0,ifelse(Metal.work2==1, 1, 0))) %>% 
  mutate(Metal.work3 = ifelse(is.na(Metal.work3),0,ifelse(Metal.work3==1, 1, 0))) %>% 
  mutate(Metal.poinsoning = ifelse(Metal.poinsoning=="Yes", 1, 0)) %>% 
  mutate(Metal.PREDIGT    = ifelse(Metal.work1|Metal.work2|Metal.work3|Metal.poinsoning, 1, 0))

# Pesticide and farm life - data only available in baseline
df_BL <- df_BL %>% 
  mutate(pesticide = ifelse(is.na(pesticide),NA,ifelse(pesticide==1, 1, 0))) %>% 
  mutate(Agriculture.work = ifelse(is.na(Agriculture.work),0,ifelse(Agriculture.work==1, 1, 0))) %>% 
  mutate(Farm.PREDIGT     = ifelse(is.na(pesticide),0,ifelse(pesticide|Agriculture.work, 1, 0)))

# Head trauma - data only available in baseline
df_BL <- df_BL %>% 
  mutate(Head.PREDIGT = case_when(
    Head.trauma == "No" ~ 0,
    Head.trauma == "Yes" & Lost.Consciousness == "Yes" ~ 2,
    TRUE ~ 1
  ))

# baseline - Sex and age group & olfaction classification
df_BL <- df_BL %>% 
  mutate(GTgroup = case_when(
    sex == "Female" & age <= 55 ~ 1,
    sex == "Female" & age > 55 ~ 2,
    sex == "Male" & age <= 55 ~ 3,
    sex == "Male" & age > 55 ~ 4
  )) %>% 
  mutate(Smell.PREDIGT = case_when(
    Sniffin.Sticks.TDI <= 16 ~ 2,
    GTgroup == 1 & Sniffin.Sticks.ID <= 12 ~ 1,
    GTgroup == 2 & Sniffin.Sticks.ID <= 9 ~ 1,
    GTgroup == 3 & Sniffin.Sticks.ID <= 11 ~ 1,
    GTgroup == 4 & Sniffin.Sticks.ID <= 9 ~ 1,
    TRUE ~ 0
  ))

# Baseline - Family History
df_BL <- df_BL %>% 
  mutate(FH.PREDIGT = case_when(
    FH1==1        ~ 4,
    FH2==1        ~ 2,
    FH3==1        ~ 1,
    TRUE          ~ 0
  )) %>% 
  mutate(FH.PREDIGT = ifelse(is.na(FH.any), NA, FH.PREDIGT)) 

# visit24 - update Family History
df_24 <- df_24 %>% 
  left_join(df_BL %>% select(subID, FH1, FH2, FH3), by="subID") %>% 
  mutate(FH1 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & (Fhchange2 == "Parent" | Fhchange2 == "Sibling"),1,FH1)) %>% 
  mutate(FH2 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & (Fhchange2 == "Grandparent" | Fhchange2 == "Uncle"),1,FH2)) %>% 
  mutate(FH3 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & Fhchange2 == "Cousin ",1,FH3)) %>% 
  mutate(FH.any = ifelse(is.na(FH1)&is.na(FH2)&is.na(FH3), NA, ifelse(FH1==1|FH2==1|FH3==1,1,0))) %>% 
  mutate(FH.PREDIGT = case_when(
    FH1==1        ~ 4,
    FH2==1        ~ 2,
    FH3==1        ~ 1,
    TRUE          ~ 0
  )) %>% 
  mutate(FH.PREDIGT = ifelse(is.na(FH.any), NA, FH.PREDIGT)) 

# visit48 - update Family History
df_48 <- df_48 %>% 
  left_join(df_24 %>% select(subID, FH1, FH2, FH3)) %>% 
  mutate(FH1 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & (Fhchange2 == "Parent" | Fhchange2 == "Sibling"),1,FH1)) %>% 
  mutate(FH2 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & (Fhchange2 == "Aunt" | Fhchange2 == "Uncle"),1,FH2)) %>% 
  mutate(FH3 = ifelse(!is.na(Fhchange1) & !is.na(Fhchange2) & Fhchange1=="Parkinson's Disease IPS/PD" & Fhchange2 == "Cousin ",1,FH3)) %>% 
  mutate(FH.any = ifelse(is.na(FH1)&is.na(FH2)&is.na(FH3), NA, ifelse(FH1==1|FH2==1|FH3==1,1,0))) %>% 
  mutate(FH.PREDIGT = case_when(
    FH1==1        ~ 4,
    FH2==1        ~ 2,
    FH3==1        ~ 1,
    TRUE          ~ 0
  )) %>% 
  mutate(FH.PREDIGT = ifelse(is.na(FH.any), NA, FH.PREDIGT))

# visit24 - update subjective sense of smell
vec <- c("PD.NMS.02", "PD.NMSS.28.Severity", "PD.NMSS.28.Frequency",
         "PD.NMS.02.old", "PD.NMSS.28.Severity.old", "PD.NMSS.28.Frequency.old")
names(df_24)[16:18] <- vec[4:6]
df_24 <- df_24 %>% 
  left_join(df_BL %>% select(subID, PD.NMS.02, PD.NMSS.28.Severity, PD.NMSS.28.Frequency), by="subID")
df_24[names(df_24)[which(names(df_24) %in% vec)]] <- 
  lapply(df_24[names(df_24)[which(names(df_24) %in% vec)]], as.character)

df_24 <- df_24 %>% 
  mutate(PD.NMS.02 = factor(ifelse((is.na(PD.NMS.02) & is.na(PD.NMS.02.old)), NA,
                                   ifelse((PD.NMS.02 == "Yes" | PD.NMS.02.old == "Yes"), "Yes", "No")))) %>% 
  mutate(PD.NMSS.28.Severity = factor(ifelse((is.na(PD.NMSS.28.Severity) & is.na(PD.NMSS.28.Severity.old)), NA,
                                             ifelse((is.na(PD.NMSS.28.Severity) |PD.NMSS.28.Severity == "0 - None"), 
                                                    PD.NMSS.28.Severity.old,PD.NMSS.28.Severity)))) %>% 
  mutate(PD.NMSS.28.Frequency = factor(ifelse((is.na(PD.NMSS.28.Frequency) & is.na(PD.NMSS.28.Frequency.old)), NA,
                                              ifelse((is.na(PD.NMSS.28.Frequency) |PD.NMSS.28.Frequency == "0 - None"),
                                                     PD.NMSS.28.Frequency.old,PD.NMSS.28.Frequency))))
df_24[names(df_24)[which(names(df_24) %in% vec)]] <- 
  lapply(df_24[names(df_24)[which(names(df_24) %in% vec)]], as.factor)

# visit48 - update subjective sense of smell
names(df_48)[16:18] <- vec[4:6]
df_48 <- df_48 %>% 
  left_join(df_24 %>% select(subID, PD.NMS.02, PD.NMSS.28.Severity, PD.NMSS.28.Frequency), by="subID") 
df_48[names(df_48)[which(names(df_48) %in% vec)]] <- 
  lapply(df_48[names(df_48)[which(names(df_48) %in% vec)]], as.character)  
df_48 <- df_48 %>% 
  mutate(PD.NMS.02 = factor(ifelse((is.na(PD.NMS.02) & is.na(PD.NMS.02.old)), NA,
                                   ifelse((PD.NMS.02 == "Yes" | PD.NMS.02.old == "Yes"), "Yes", "No")))) %>% 
  mutate(PD.NMSS.28.Severity = factor(ifelse((is.na(PD.NMSS.28.Severity) & is.na(PD.NMSS.28.Severity.old)), NA,
                                             ifelse((is.na(PD.NMSS.28.Severity) |PD.NMSS.28.Severity == "0 - None"), 
                                                    PD.NMSS.28.Severity.old,PD.NMSS.28.Severity)))) %>% 
  mutate(PD.NMSS.28.Frequency = factor(ifelse((is.na(PD.NMSS.28.Frequency) & is.na(PD.NMSS.28.Frequency.old)), NA,
                                              ifelse((is.na(PD.NMSS.28.Frequency) |PD.NMSS.28.Frequency == "0 - None"),
                                                     PD.NMSS.28.Frequency.old,PD.NMSS.28.Frequency))))
df_48[names(df_48)[which(names(df_48) %in% vec)]] <- 
  lapply(df_48[names(df_48)[which(names(df_48) %in% vec)]], as.factor)

# visit48 - Sex and age group & olfaction classification
df_48 <- df_48 %>% 
  mutate(GTgroup = case_when(
    sex == "Female" & age <= 55 ~ 1,
    sex == "Female" & age > 55 ~ 2,
    sex == "Male" & age <= 55 ~ 3,
    sex == "Male" & age > 55 ~ 4
  )) %>% 
  mutate(Smell.PREDIGT = case_when(
    #Sniffin.Sticks.TDI <= 16 ~ 2, # No TDI score in the curated form
    GTgroup == 1 & Sniffin.Sticks.ID <= 12 ~ 1,
    GTgroup == 2 & Sniffin.Sticks.ID <= 9 ~ 1,
    GTgroup == 3 & Sniffin.Sticks.ID <= 11 ~ 1,
    GTgroup == 4 & Sniffin.Sticks.ID <= 9 ~ 1,
    TRUE ~ 0
  ))

# transform categorical variables into numberic
vfct <- vfct[!vfct %in% c("groupA","groupB","groupC","sex")]

df_BL[names(df_BL)[which(names(df_BL) %in% vfct)]] <- 
  lapply(df_BL[names(df_BL)[which(names(df_BL) %in% vfct)]], as.numeric.factor2)
df_24[names(df_24)[which(names(df_24) %in% vfct)]] <- 
  lapply(df_24[names(df_24)[which(names(df_24) %in% vfct)]], as.numeric.factor2)
df_48[names(df_48)[which(names(df_48) %in% vfct)]] <- 
  lapply(df_48[names(df_48)[which(names(df_48) %in% vfct)]], as.numeric.factor2)

# ------------ combine data frames, only select PD and HC ---------------
df_BL <- df_BL %>% mutate(EVENT_ID = "baseline") %>% mutate(subID = as.character(subID))
df_24 <- df_24 %>% mutate(EVENT_ID = "24 months") %>% mutate(subID = as.character(subID))
df_48 <- df_48 %>% mutate(EVENT_ID = "48 months") %>% mutate(subID = as.character(subID))
data_join  <- df_BL %>% bind_rows(df_24) %>% bind_rows(df_48)
data_join[grep('group', names(data_join))] <- lapply(data_join[grep('group', names(data_join))],as.character)
data_join[grep("OND", data_join$groupC), c("groupA","groupB")] = "Other Neurological Patient OND"

data_join  <-  data_join %>% 
  mutate(group = factor(case_when(
    EVENT_ID == "baseline"  ~ groupA,
    EVENT_ID == "24 months" ~ groupB,
    EVENT_ID == "48 months" ~ groupC
  ))) %>% 
  filter(group == "Healthy Control HC" | group == "Parkinson's Disease Patient PD") %>% 
  mutate(EVENT_ID = factor(EVENT_ID)) %>% 
  droplevels() %>% 
  mutate(EVENT_ID = factor(EVENT_ID, levels = c("baseline", "24 months", "48 months"))) %>% 
  arrange(EVENT_ID) %>% 
  arrange(subID)

levels(data_join$group) <- group_name
names(data_join)[1]     <- "ID"

# ======= Tabel 1: Demographic summary ======
# Hoehn And Yahr Stage
data_BL3 <- read.csv(paste(data_path, "/", "Denopa - BL data.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))

df_HY <- data_BL3 %>% select("ï»¿Probanden_Nr", "a_H_Y")
names(df_HY) <- c("ID", "HYstage")

data_T1 <- data_join %>% 
  filter(EVENT_ID == "baseline") %>% 
  droplevels() %>% 
  select(ID, group, sex, age) %>% 
  left_join(df_HY, by = "ID") %>% 
  mutate(age_fct = factor(case_when(
    age < 40 ~ "<40",
    age < 50  ~ "40-50",
    age < 60  ~ "50-60",
    age < 70  ~ "60-70",
    age < 80  ~ "70-80",
    TRUE  ~ "80+"
  ))) 

data_T1 %>% 
  group_by(group, HYstage) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n/sum(n)*100)

summary(data_T1)
data_T1 %>% group_by(group,sex) %>% summarise(n())
kruskal.test(sex ~ group, data = data_T1)

data_T1 %>% group_by(group) %>% 
  summarise(min=min(age),
            median=median(age),
            max=max(age))
kruskal.test(age ~ group, data = data_T1)

data_T1 %>% group_by(group,sex) %>% 
  summarise(min=min(age),
            median=median(age),
            max=max(age))
kruskal.test(age ~ group, data = data_T1 %>% filter(sex=="Male"))
kruskal.test(age ~ group, data = data_T1 %>% filter(sex=="Female"))

nHC  <- sum(data_T1$group == group_name[1])
nPD  <- sum(data_T1$group == group_name[2])
data_T1 %>% group_by(group,age_fct) %>% summarise(count=n(), 
                                                  prop1=round(count/nHC*100,2), 
                                                  prop2=round(count/nPD*100,2))

# ========== imputate missing data using previous or next non-NA value ========== 
use_imp <- TRUE
if (use_imp)
{
  data <- data.impute(data_join)
} else
{
  data <- data_join
}

# ========== integrated score of each variable ========== 
varList <- list(
  c("MDS.UPDRS.1.11","SCOPA.AUT.05","SCOPA.AUT.06","PD.NMS.05","PD.NMSS.21.Severity","PD.NMSS.21.Frequency"),#
  c("MDS.UPDRS.1.03", "PD.NMS.16","PD.NMSS.10.Severity","PD.NMSS.10.Frequency","PDQ39.17", "BDI", "MADRS"),
  c("MDS.UPDRS.1.04", "PD.NMS.17", "PD.NMSS.09.Severity", "PD.NMSS.09.Frequency", "PDQ39.21"),
  c("PD.NMS.25", "RBD.SQ","RBD.diagnosis"),
  c("Alpha.Synuclein"),
  c("tTau")
)
varStr2 <- c("Constipation", "Depression", "Anxiety", "RBD","Alpha.Synuclein","tTau")
EVENT_vec <- c("baseline", "24 months", "48 months")

df_LR <- varScore(data, varList, varStr2, EVENT_vec, group_name)

# ========== univariate results =========
df_LR <- df_LR %>% 
  mutate(RBD.SQ.06.1to4 = select(., c(RBD.SQ.06.1:RBD.SQ.06.4)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(RBD.SQ.06.1to2 = select(., c(RBD.SQ.06.1,RBD.SQ.06.2)) %>% apply(1, sum, na.rm=TRUE)) 
varStr <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT",
            "Constipation.PREDIGT","MDS.UPDRS.1.11","SCOPA.AUT.05","SCOPA.AUT.06","PD.NMS.05",
            "PD.NMSS.21.Severity","PD.NMSS.21.Frequency",
            "Smell.PREDIGT","Sniffin.Sticks.ID","Sniffin.Sticks.TDI","PD.NMS.02","PD.NMSS.28.Severity", "PD.NMSS.28.Frequency",
            "Smoking.PREDIGT", "Caffeine.PREDIGT",
            "FH.PREDIGT", "FH1","FH2","FH.any",
            "Alpha.Synuclein","tTau",
            "Depression.PREDIGT", "MDS.UPDRS.1.03","PD.NMS.16","PD.NMSS.10.Severity","PD.NMSS.10.Frequency","PDQ39.17",
            "GDS","BDI","MADRS",
            "Anxiety.PREDIGT", "MDS.UPDRS.1.04","PD.NMS.17","PD.NMSS.09.Severity","PD.NMSS.09.Frequency","PDQ39.21",
            "RBD.PREDIGT","RBD.SQ","PD.NMS.25","RBD.diagnosis","sex","age",
            "RBD.SQ.01","RBD.SQ.02.","RBD.SQ.03.","RBD.SQ.04.","RBD.SQ.05.","RBD.SQ.06.1","RBD.SQ.06.2",         
            "RBD.SQ.06.3","RBD.SQ.06.4","RBD.SQ.07.","RBD.SQ.08.","RBD.SQ.09.","RBD.SQ.10.",
            "RBD.SQ.06.1to4","RBD.SQ.06.1to2")

#results <- uniVar(df_LR, varStr, EVENT_vec, group_name)

pltStr  <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT",
             "Smell.PREDIGT", "Smoking.PREDIGT", "Caffeine.PREDIGT",
             "FH.PREDIGT", "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
pltStr2  <- c("E: Metal", "E: Pesticide", "E: Head Trauma", "E: Constipation",
             "E: Hyposmia", "E: Smoking", "E: Caffeine",
             "D: Family History", "I: Depression", "I: Anxiety", "I: RBD")
results <- uniVar(df_LR, pltStr, EVENT_vec, group_name, pltStr = pltStr, pltStr2 = pltStr2)
results[[1]]

# ==========  Add Genetic Info ========== 
df_LR2 <- df_LR %>% 
  left_join(genetic, by = c("ID" = "Code"))
df_LR2_BL <- df_LR2 %>% filter(EVENT_ID == "baseline")

plotBar(df_LR2_BL, "MLPA", "MLPA")

names(df_LR2_BL)[114] <- "GBA"
plotBar(df_LR2_BL, "GBA", "GBA exons 8-11")
df_LR2_BL <- df_LR2_BL %>% 
  mutate(GBA_fct = ifelse(stringr::str_starts(GBA, "No"),"No mutation", "Mutation")) %>% 
  mutate(GBA_fct = factor(GBA_fct, levels = c("No mutation", "Mutation")))
plotBar(df_LR2_BL, "GBA_fct", "GBA", or_cal = TRUE)

names(df_LR2_BL)[115] <- "rs11931074"
df_LR2_BL <- df_LR2_BL %>% 
  mutate(rs11931074 = factor(rs11931074, levels = c("G/G", "G/T", "T/T")))
plotBar(df_LR2_BL, "rs11931074", "rs11931074 SNP", or_cal = TRUE, or_combine = TRUE)

names(df_LR2_BL)[116] <- "TAU"
df_LR2_BL <- df_LR2_BL %>% 
  mutate(TAU = factor(TAU, levels = c("H1/H1", "H1/H2", "H2/H2", "not enough gDNA")))
plotBar(df_LR2_BL, "TAU", "H0/H1 in TAU", or_cal = TRUE, or_combine = TRUE)

names(df_LR2_BL)[117] <- "REP1"
plotBar(df_LR2_BL, "REP1", "REP1 in 5'UTR of SNCA", or_cal = TRUE, or_combine = TRUE)

plotBar(df_LR2_BL, "APOE", "APOE", or_cal = TRUE, or_combine = TRUE)

# ========== save results and combine with PPMI ==========
write.csv(df_LR,paste(data_path,"/DeNoPa dfLR.csv",sep=""), row.names = FALSE)
# ========== calculate PREDIGT Score =========
FG  <- TRUE # including factor G?
FT  <- TRUE # including factor T?
df_LR <- df_LR %>% mutate(score=0)

theme_set(
  theme_classic()
)
My_Theme = theme(
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  axis.text.y = element_text(size = 25))

BL_only <- TRUE
AUC_only <- FALSE
print_result <- FALSE
plt <- TRUE
plot_x <- "cut1" # "n", "prop"
box_only <- TRUE
label_size <- 16
font_size  <- 30 # legend in density plot
font_size2  <- 10 # legend in ROC plot

# A colorblind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# PPV and NPV
prevalence <- 1/100 # prevalence
plotPPV <- TRUE

# compare AUC in Female vs Male
sex_ind <- 0 # 0: Both; 1: Female; 2: Male

# ======= select variables =========
# ------- 1. DeNoPa original -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT", "Smell.PREDIGT",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT","Alpha.Synuclein.PREDIGT","tTau.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa original"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 10: alpha-synuclein: 0.25,0.5,1
           0.25,    # 11: T-tau: 0.25,0.5,1
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:7)
ind_D <- 8
ind_I <- c(9:13)

# ------- 2. DeNoPa no CSF -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT", "Smell.PREDIGT",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa no CSF"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:7)
ind_D <- 8
ind_I <- c(9:11)

# ------- 2.5 DeNoPa no CSF, no E -------
varStr2 <- c("Constipation.PREDIGT", "Smell.PREDIGT",
             "FH.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa no CSF, no E"
# Coeffcients
coefV <- c(0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1,2)
ind_D <- 3
ind_I <- c(4:6)
# ------- 3. self-report -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "SCOPA.AUT.06", "PD.NMS.02",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT",
             "PDQ39.17", "PDQ39.21", "PD.NMS.25")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa self-report"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- c(4,5,9:11)
ind_E <- c(1:7)
ind_D <- 8
ind_I <- c(9:11)
# ------- 4. self-report without E vars -----------
varStr2 <- c("SCOPA.AUT.06", "PD.NMS.02",
             "FH.PREDIGT",
             "PDQ39.17", "PDQ39.21", "PD.NMS.25")

varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa self-report, no E"
# Coeffcients
coefV <- c(0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- c(1,2,4:6)
ind_E <- c(1:7)
ind_D <- 8
ind_I <- c(9:11)

# ------- 5. self-report, PPMI items -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "SCOPA.AUT.06",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT",
             "MDS.UPDRS.1.03", "MDS.UPDRS.1.04","RBD.SQ.06.1to2")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa self-report, PPMI"
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- c(4,8:10)
ind_E <- c(1:6)
ind_D <- 7
ind_I <- c(8:10)
# ------- 6. self-report without E vars, PPMI items -----------
varStr2 <- c("SCOPA.AUT.06",
             "FH.PREDIGT",
             "MDS.UPDRS.1.03", "MDS.UPDRS.1.04","RBD.SQ.06.1to2")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "DeNoPa self-report, PPMI, no E"
coefV <- c(0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)  
lenV <- length(coefV)
ind_question <- c(1,3:5)
ind_E <- 1
ind_D <- 2
ind_I <- c(3:5)

# ====== Calculate PREDIGT Score ===========
df_LR <- PREDIGTScr(df_LR, varStr1, varStr2, coefV, ind_question, EVENT_vec, FG, FT, method = 1)
# method = 1: Score with the original coefficients, w age-associated coefficients
# method = 2: score using simple additive LR *G*T

# ====== Plot: ROC ==============
save_plot <- TRUE
ROC_result <- plotROC(df_LR, EVENT_vec, group_name, plot_path, tableStr, My_Theme, save_plot)
cut_BL <- ROC_result$cut_BL
print(ROC_result$result)

df      <- df_LR %>% filter(!is.na(score)) %>% filter(EVENT_ID=="baseline")
scoreHC <- df %>% filter(group==group_name[1]) %>% select(score) 
scorePD <- df %>% filter(group==group_name[2]) %>% select(score) 
print(paste("HC: ", round(mean(scoreHC$score), 2), " (+/- ", round(sd(scoreHC$score), 2), ") ",
            "PD: ", round(mean(scorePD$score), 2), " (+/- ", round(sd(scorePD$score), 2), ") ", sep=""))

# ====== Plot: distribution and box plot ==========
cutline <- TRUE
violin <- FALSE
adj <- 1.2
xmax <- 300
plotDistBox(df_LR, cut_BL, EVENT_vec, group_name, plot_path, tableStr, My_Theme, save_plot, xmax = xmax) 

# ====== 2-step ===============
step2 <- c("Sniffin.Sticks.TDI",
           "GDS", "BDI", "MADRS",
           "RBD.SQ",
           "Alpha.Synuclein",
           "tTau")
Model2 <- model2(df_LR, step2, group_name)

# Step 1 only AUC
Model2$step1List[[1]]
Model2$step1List[[3]]
Model2$step1List[[6]]

# Step 2 only AUC
Model2$step2List[[1]][[1]]
Model2$step2List[[3]][[1]]
Model2$step2List[[6]][[1]]

# combined AUC
aucV <- Model2$resultList[[1]]$auc_combine
ind1 <- which(aucV == max(aucV, na.rm=TRUE))

# ====== 2-step: scatter plot ========
indTest  <- c(1, 3, 6)
model2point(Model2, df_LR, step2, indTest, plot_path, tableStr, My_Theme, save_plot)

# ====== 2-step: ROC plot =======
indTest  <- c(1, 3, 6)
model2ROC(Model2, df_LR, step2, indTest, plot_path, tableStr, My_Theme, save_plot)

# ====== 2-step: boxplot =======
indTest  <- 1:length(step2) # index of step 2 tests that will be include into the boxplot
model2box(Model2, step2, indTest, plot_path, tableStr, My_Theme, save_plot)

# ====== 2-step: bootstrapping =========
indTest  <- c(1, 3, 6)
set.seed(100)
reps  <- 1000
boot_results <- model2boot(df_LR, step2, indTest, reps)
# save(boot_results,file="DeNoPa bootstrapping.Rdata")
# 
# load("DeNoPa bootstrapping.Rdata")
for (istep in 1:length(indTest))
{
  print(boot.ci(boot_results[[istep]], type="perc")) #c("norm","basic", "stud", "perc", "bca")
}




# ====== correlation with motor ==========
motorPlt <- function(df_plt, ttlStr, pHC = FALSE)
{
  df <- df_plt %>% 
    group_by(HYstage) %>% 
    summarise(n=n())
  stage <- as.character(df$HYstage[df$n > 3])
  
  p1 <- ggplot(df_plt %>% filter(HYstage %in% stage) %>% droplevels(), aes(HYstage, score, fill = HYstage))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    stat_compare_means()+
    labs(title = ttlStr,
         y = "PREDIGT Score")
  print(p1)
  
  p2 <- ggplot(df_plt, aes(score, MDS_UPDRS_3))+
    geom_point()+
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = ttlStr,
         x = "PREDIGT Score")
  ymaxd <- layer_scales(p2)$y$range$range[2]  
  p2 <- p2 + stat_regline_equation(label.y = ymaxd+5, aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = ymaxd, aes(label = ..rr.label..))
  if (pHC) p2 <- p2+ facet_wrap(~group)
  print(p2)
  
  p3 <- ggplot(df_plt, aes(score, MDS_UPDRS_3, color = HYstage))+ #
    geom_point(alpha=0.5)+
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = ttlStr,
         x = "PREDIGT Score") 
  if (pHC) p3 <- p3+ facet_wrap(~group)
  print(p3)
}
# Hoehn And Yahr Stage
data_BL3 <- read.csv(paste(data_path, "/", "Denopa - BL data.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_motor <- data_BL3 %>% select("ï»¿Probanden_Nr", "a_H_Y")
names(df_motor) <- c("ID", "HYstage")

df_motor <- df_motor %>% 
  left_join(data_BL %>% select("Patient ID", "MDS UPDRS 2008: 3 Total Score: Sum Of: 3.1 - 3.18", 
                               "MDS UPDRS 2008: 4 Total Score: Sum Of: 4.1 - 4.6"), by = c("ID" = "Patient ID"))
names(df_motor)[c(3,4)] <- c("MDS_UPDRS_3", "MDS_UPDRS_4")

df_motor <- df_LR %>% 
  filter(EVENT_ID == "baseline") %>% 
  select(ID, group, score) %>% 
  left_join(df_motor, by = "ID") %>% 
  filter(!is.na(score)) %>% 
  mutate(HYstage = factor(HYstage))

df_HY <- df_motor %>% 
  group_by(HYstage) %>% 
  summarise(n=n())
levels(df_motor$HYstage) <- paste(df_HY$HYstage, " (n=", df_HY$n, ")", sep="")

ttlStr <- "DeNoPa Model1'"
df_plt <- df_motor
pHC    <- TRUE

df_plt <- df_motor %>% filter(group == group_name[2]) %>% droplevels()
pHC    <- FALSE

motorPlt(df_plt, ttlStr, pHC)
