# Main funtion for PREDIGT analysis on PPMI
# Created by Juan Li, 2022-03-22

library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)
library(ggpubr)
library(rms)
library(pracma)
library(boot)

FOUND_only <- TRUE
clean <- TRUE # only PD and HC

combine_DeNoPa <- FALSE

# source functions
wd <- "/Users/juan/Desktop/OHRI/work/Cohorts/! PREDIGT functions"
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
data_path  <- "/Users/juan/Desktop/OHRI/work/Cohorts/!DeNoPa PPMI Temporal/PPMI"
setwd(data_path)
plot_path  <- paste(data_path, "/Plots", sep="")
group_name <- c("Healthy Control", "Parkinson's Disease")
# ------------- Read data -----------
# read in all data
data_V   <- read.csv("variable list.csv", header = F, stringsAsFactors=FALSE, fileEncoding="latin1", na.strings=c("",".","NA"))

df_Demographics  <- read.csv("Screening___Demographics.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_FH            <- read.csv("Family_History__PD.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_DNA           <- read.csv("DNA.csv",header = T, fileEncoding="latin1", check.names=TRUE, na.strings=c("",".","NA"))
df_GDS           <- read.csv("Geriatric_Depression_Scale__Short.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_MDSUPDRSI1    <- read.csv("MDS_UPDRS_Part_I__Patient_Questionnaire-2.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_MDSUPDRSI2    <- read.csv("MDS_UPDRS_Part_I-2.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_MDSUPDRSI     <- df_MDSUPDRSI1 %>% full_join(df_MDSUPDRSI2, by=c("PATNO","EVENT_ID"))
df_RBDQ          <- read.csv("REM_Sleep_Disorder_Questionnaire.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_SCOPA_AUT     <- read.csv("SCOPA-AUT.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_STAI          <- read.csv("State-Trait_Anxiety_Inventory.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_UPSIT         <- read.csv("University_of_Pennsylvania_Smell_ID_Test.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_CSF           <- read.csv("labs02.csv",header = T, fileEncoding="latin1", check.names=TRUE, na.strings=c("",".","NA"))
df_motor         <- read.csv("motor04.csv",header = T, fileEncoding="latin1", check.names=TRUE, na.strings=c("",".","NA"))
# FOUND
df_Caffeine      <- read.csv("FOUND_RFQ_Caffeine.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Head          <- read.csv("FOUND_RFQ_Head_Injury.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Pesticide1    <- read.csv("FOUND_RFQ_Pesticides_at_Work.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Pesticide2    <- read.csv("FOUND_RFQ_Pesticides_Non-Work.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Residential   <- read.csv("FOUND_RFQ_Residential_History.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Pesticide     <- df_Pesticide1 %>% full_join(df_Pesticide2, by="patno")  %>% full_join(df_Residential, by="patno") 
df_Smoking       <- read.csv("FOUND_RFQ_Smoking_History.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
df_Toxicant      <- read.csv("FOUND_RFQ_Toxicant_History.csv",header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))

# ------------ data prepatation -------------------
# ---------- PPMI --------------
Demographics <- df_Demographics %>% 
  filter(is.na(DECLINED) & is.na(EXCLUDED)) %>% 
  select(PATNO, APPRDX, CURRENT_APPRDX, BIRTHDT, GENDER) %>% 
  mutate(gender = factor(case_when(
    GENDER == 2   ~ "Male",
    TRUE          ~ "Female"
  ))) %>% 
  mutate(group1 = factor(case_when(
    APPRDX == 1   ~ "Parkinson disease",
    APPRDX == 2   ~ "Healthy Control",
    APPRDX == 3   ~ "SWEDD",
    APPRDX == 4   ~ "Prodromal",
    APPRDX == 5   ~ "Genetic Cohort - PD",
    APPRDX == 6   ~ "Genetic Cohort - Unaffected",
    APPRDX == 7   ~ "Genetic Registry - PD",
    APPRDX == 8   ~ "Genetic Registry - Unaffected"
  ))) %>% 
  mutate(group2 = factor(case_when(
    CURRENT_APPRDX == 1   ~ "Parkinson disease",
    CURRENT_APPRDX == 2   ~ "Healthy Control",
    CURRENT_APPRDX == 3   ~ "SWEDD",
    CURRENT_APPRDX == 4   ~ "Prodromal",
    CURRENT_APPRDX == 5   ~ "Genetic Cohort - PD",
    CURRENT_APPRDX == 6   ~ "Genetic Cohort - Unaffected",
    CURRENT_APPRDX == 7   ~ "Genetic Registry - PD",
    CURRENT_APPRDX == 8   ~ "Genetic Registry - Unaffected"
  ))) %>% 
  mutate(group3 = factor(case_when(
    APPRDX == 1 | APPRDX == 5 | APPRDX == 7  ~ "Parkinson disease",
    APPRDX == 2 | APPRDX == 6 | APPRDX == 8  ~ "Healthy Control",
    APPRDX == 3   ~ "SWEDD",
    APPRDX == 4   ~ "Prodromal"
  ))) %>% 
  mutate(group4 = factor(case_when(
    CURRENT_APPRDX == 1 | CURRENT_APPRDX == 5 | CURRENT_APPRDX == 7   ~ "Parkinson disease",
    CURRENT_APPRDX == 2 | CURRENT_APPRDX == 6 | CURRENT_APPRDX == 8   ~ "Healthy Control",
    CURRENT_APPRDX == 3   ~ "SWEDD",
    CURRENT_APPRDX == 4   ~ "Prodromal"
  )))

FH <- df_FH %>% 
  select(PATNO, EVENT_ID, INFODT, BIOMOM:KIDSPD) %>% 
  separate(INFODT, c("infomonth", "FH_infoyear"),remove = FALSE) %>% 
  mutate(KIDSPD   = ifelse(KIDSNUM==0, 0, ifelse(is.na(KIDSPD),NA,KIDSPD))) %>% 
  mutate(FULSIBPD = ifelse(FULSIB==0, 0, ifelse(is.na(FULSIBPD),NA,FULSIBPD))) %>% 
  mutate(HAFSIBPD = ifelse(HAFSIB==0, 0, ifelse(is.na(HAFSIBPD),NA,HAFSIBPD))) %>% 
  mutate(MAGPARPD = ifelse(MAGPAR==0, 0, ifelse(is.na(MAGPARPD),NA,MAGPARPD))) %>% 
  mutate(PAGPARPD = ifelse(PAGPAR==0, 0, ifelse(is.na(PAGPARPD),NA,PAGPARPD))) %>% 
  mutate(MATAUPD  = ifelse(MATAU==0, 0, ifelse(is.na(MATAUPD),NA,MATAUPD))) %>% 
  mutate(PATAUPD  = ifelse(PATAU==0, 0, ifelse(is.na(PATAUPD),NA,PATAUPD))) %>% 
  mutate(FH1      = case_when(
    BIOMOMPD==1 | BIODADPD==1 | KIDSPD==1 | FULSIBPD==1 ~ 1,
    TRUE  ~ 0
  )) %>%
  mutate(FH2      = case_when(
    HAFSIBPD==1 | MAGPARPD==1 | PAGPARPD==1 | MATAUPD==1 | PATAUPD==1 ~ 1,
    TRUE  ~ 0
  )) %>% 
  mutate(FH.any = ifelse(FH1==1|FH2==1, 1, 0)) %>% 
  mutate(FH.PREDIGT = ifelse(FH1==1, 4, ifelse(FH2==1, 2, 0))) %>% 
  mutate(FH1_fct = factor(FH1)) %>% 
  mutate(FH2_fct = factor(FH2)) %>% 
  mutate(FH.any_fct = factor(FH.any)) %>% 
  mutate(FH.PREDIGT_fct = factor(FH.PREDIGT))
FH2 <- FH %>% 
  select(PATNO, FH1, FH2, FH.any, FH.PREDIGT) %>% 
  group_by(PATNO) %>% 
  summarise(FH1=last(FH1),
            FH2=last(FH2),
            FH.any=last(FH.any),
            FH.PREDIGT=last(FH.PREDIGT)) 
data_unchange <- Demographics %>% left_join(FH2) # data that do not change much in each visit

GDS <- df_GDS %>% 
  select(PATNO, EVENT_ID, INFODT, 
         GDSSATIS:GDSBETER) %>% 
  separate(INFODT, c("infomonth", "GDS_infoyear"),remove = FALSE) %>% 
  mutate(GDSSATIS = ifelse(is.na(GDSSATIS), NA, 1-GDSSATIS)) %>% 
  mutate(GDSGSPIR = ifelse(is.na(GDSGSPIR), NA, 1-GDSGSPIR)) %>% 
  mutate(GDSHAPPY = ifelse(is.na(GDSHAPPY), NA, 1-GDSHAPPY)) %>% 
  mutate(GDSALIVE = ifelse(is.na(GDSALIVE), NA, 1-GDSALIVE)) %>% 
  mutate(GDSENRGY = ifelse(is.na(GDSENRGY), NA, 1-GDSENRGY)) %>% 
  mutate(GDS = select(., GDSSATIS:GDSBETER) %>% apply(1, sum, na.rm=TRUE))
# GDS$na_count <- rowSums( is.na( GDS[,5:19]))
# nrow(GDS %>% filter(na_count>0))
# in total 21 with incomplete responses, mostly 12-15, only one 6 response. Assume NA=0

MDS_UPDRS_I <- df_MDSUPDRSI %>% 
  select(PATNO, EVENT_ID, INFODT.x, NP1SLPN:NP1FATG, NP1COG:NP1DDS) %>% 
  rename(INFODT=INFODT.x) %>% 
  separate(INFODT, c("infomonth", "UPDRS_infoyear"),remove = FALSE)
data_temporal <- GDS %>% 
  select(PATNO, EVENT_ID, GDS_infoyear, GDS) %>% 
  full_join(MDS_UPDRS_I %>% select(PATNO, EVENT_ID, UPDRS_infoyear, NP1CNST, NP1DPRS, NP1ANXS), by=c("PATNO", "EVENT_ID"))

RBD.SQ <- df_RBDQ %>% 
  select(PATNO, EVENT_ID, INFODT, DRMVIVID:CNSOTH) %>% 
  separate(INFODT, c("infomonth", "RBDQ_infoyear"),remove = FALSE) %>% 
  mutate(NSDISEASE = select(., STROKE:CNSOTH) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(NSDISEASE = ifelse(NSDISEASE>0, 1, 0)) %>% 
  mutate(RBD.SQ = select(., c(DRMVIVID:SLPDSTRB,NSDISEASE)) %>% apply(1, sum, na.rm=TRUE))
# RBD.SQ$na_count <- rowSums( is.na( RBD.SQ[,c(6:17,27)]))
# nrow(RBD.SQ %>% filter(na_count>0))
# in total 22 with incomplete responses. Assume NA=0
data_temporal <- data_temporal %>% 
  full_join(RBD.SQ %>% select(PATNO, EVENT_ID, RBDQ_infoyear, RBD.SQ, DRMVIVID:CNSOTH), by=c("PATNO", "EVENT_ID"))

SCOPA_AUT <- df_SCOPA_AUT %>% 
  select(PATNO, EVENT_ID, INFODT, SCAU1:SCAU23A, SCAU24:SCAU26A, SCAU26B, SCAU26C, SCAU26D) %>% 
  separate(INFODT, c("infomonth", "SCOPAAUT_infoyear"),remove = FALSE)
data_temporal <- data_temporal %>% 
  full_join(SCOPA_AUT %>% select(PATNO, EVENT_ID, SCOPAAUT_infoyear, SCAU5, SCAU6), by=c("PATNO", "EVENT_ID"))

STAI <- df_STAI %>% 
  select(PATNO, EVENT_ID, INFODT, STAIAD1:STAIAD40) %>% 
  separate(INFODT, c("infomonth", "STAI_infoyear"),remove = FALSE) 
STAI[,c("STAIAD1", "STAIAD2", "STAIAD5", "STAIAD8", "STAIAD10", "STAIAD11", "STAIAD15", "STAIAD16", "STAIAD19", "STAIAD20",
        "STAIAD21", "STAIAD23", "STAIAD26", "STAIAD27", "STAIAD30", "STAIAD33", "STAIAD34", "STAIAD36", "STAIAD39")] <-
  5-STAI[,c("STAIAD1", "STAIAD2", "STAIAD5", "STAIAD8", "STAIAD10", "STAIAD11", "STAIAD15", "STAIAD16", "STAIAD19", "STAIAD20",
            "STAIAD21", "STAIAD23", "STAIAD26", "STAIAD27", "STAIAD30", "STAIAD33", "STAIAD34", "STAIAD36", "STAIAD39")]
STAI <- STAI %>%   
  mutate(STAIS = select(., STAIAD1:STAIAD20) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(STAIT = select(., STAIAD21:STAIAD40) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(STAI  = select(., STAIAD1:STAIAD40) %>% apply(1, sum, na.rm=TRUE))  
STAI$na_countS <- rowSums( is.na( STAI[,c(6:25)]))
STAI$na_countT <- rowSums( is.na( STAI[,c(26:45)]))
STAI$na_count <- rowSums( is.na( STAI[,c(6:45)]))
STAI <- STAI %>%   
  mutate(STAIS = STAIS + na_countS) %>% 
  mutate(STAIT = STAIT + na_countT) %>% 
  mutate(STAI  = STAI  + na_count) 
# STAI$na_count <- rowSums( is.na( STAI[,c(6:45)])): 
# in total 53 with incomplete responses. Assume NA=1
data_temporal <- data_temporal %>% 
  full_join(STAI %>% select(PATNO, EVENT_ID, STAI_infoyear, STAIS, STAIT, STAI), by=c("PATNO", "EVENT_ID"))

UPSIT <- df_UPSIT %>% 
  select(PATNO, EVENT_ID, INFODT, UPSITBK1:NORMATIVE_SCORE) %>% 
  separate(INFODT, c("infomonth", "UPSIT_infoyear"),remove = FALSE) %>% 
  mutate(UPSIT = select(., UPSITBK1:UPSITBK4) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(Smell.PREDIGT = case_when(
    NORMATIVE_SCORE == "Anosmia" ~ 2,
    NORMATIVE_SCORE == "Hyposmia" ~ 1,
    NORMATIVE_SCORE == "Normosmia" ~ 0,
    TRUE ~ -1
  )) 
UPSIT$Smell.PREDIGT[UPSIT$Smell.PREDIGT==-1]  <- NA
UPSIT <- UPSIT %>% 
  mutate(Smell.PREDIGT.fct = factor(Smell.PREDIGT))
data_temporal <- data_temporal %>% 
  full_join(UPSIT %>% select(PATNO, EVENT_ID, UPSIT_infoyear, UPSITBK1:UPSITBK4, UPSIT, Smell.PREDIGT), by=c("PATNO", "EVENT_ID"))

CSF <- df_CSF %>% 
  select(SUBJ_ID, VISIT_NAME, Test.Type..TESTNAME., Test.Value..TESTVALUE.) %>% 
  rename(PATNO=SUBJ_ID, Test.Type=Test.Type..TESTNAME., Test.Value = Test.Value..TESTVALUE.) %>% 
  filter(Test.Type == "CSF Alpha-synuclein (pg/ml)" |
           Test.Type == "p-Tau181P (pg/ml)" | 
           Test.Type == "Total tau (pg/ml)") %>% 
  pivot_wider(names_from = Test.Type, values_from = Test.Value) %>% 
  rename(Alpha.Synuclein="CSF Alpha-synuclein (pg/ml)", 
         pTau="p-Tau181P (pg/ml)", 
         tTau = "Total tau (pg/ml)") %>% 
  mutate(Alpha.Synuclein=as.numeric(Alpha.Synuclein), 
         pTau=as.numeric(pTau), 
         tTau = as.numeric(tTau)) %>% 
  mutate(EVENT_ID = factor(case_when(
    VISIT_NAME == "Baseline" ~ "BL",
    VISIT_NAME == "Baseline (2nd)" ~ "BL",
    VISIT_NAME == "Baseline (3d)" ~ "BL",
    VISIT_NAME == "Month 03" ~ "V01",
    VISIT_NAME == "Month 06" ~ "V02",
    VISIT_NAME == "Month 06 (2nd)" ~ "V02",
    VISIT_NAME == "Month 12" ~ "V04",
    VISIT_NAME == "Month 12 (2nd)" ~ "V04",
    VISIT_NAME == "Month 24" ~ "V06",
    TRUE ~"other"
  )))
data_temporal <- data_temporal %>% 
  full_join(CSF %>% select(PATNO, EVENT_ID, Alpha.Synuclein, pTau, tTau), by=c("PATNO", "EVENT_ID"))

DNA <- df_DNA %>% 
  select(SUBJ_ID, VISIT_NAME, Test.Type..TESTNAME., Test.Value..TESTVALUE.) %>% 
  rename(PATNO=SUBJ_ID, Test.Type=Test.Type..TESTNAME., Test.Value = Test.Value..TESTVALUE.) %>% 
  filter(Test.Type == "Genetic Risk Score") %>% 
  rename(GRS=Test.Value) %>% 
  mutate(GRS=as.numeric(GRS)) %>% 
  mutate(EVENT_ID = factor(case_when(
    VISIT_NAME == "Screening Visit" ~ "SC",
    VISIT_NAME == "Baseline" ~ "BL",
    VISIT_NAME == "Month 03" ~ "V01",
    VISIT_NAME == "Month 06" ~ "V02",
    VISIT_NAME == "Month 09" ~ "V03",
    VISIT_NAME == "Month 12" ~ "V04",
    TRUE ~"other"
  )))
data_temporal <- data_temporal %>% 
  full_join(DNA %>% select(PATNO, EVENT_ID, GRS), by=c("PATNO", "EVENT_ID"))


# -------- FOUND-------- 
Caffeine <- df_Caffeine %>% 
  select(patno, caffeine_timestamp, 
         cfqa1, cfa3, cfqa5day, cfqa5week,
         cfqb1, cfb3, cfqb5day, cfqb5week,
         cfqc1, cfc3, cfqc5day, cfqc5week) %>% 
  rename(PATNO=patno) %>% 
  separate(caffeine_timestamp,c("infomonth", "infoday","Caffeine_infoyear","infohour","infomin"),remove = FALSE) 
Caffeine[Caffeine==9999] <- NA  
Caffeine[Caffeine==7777] <- NA  
Caffeine <- Caffeine %>%  
  mutate(cfa3      = ifelse(cfqa1==0, 0, ifelse(is.na(cfa3), NA, cfa3))) %>% 
  mutate(cfqa5week = ifelse(cfqa1==0, 0, ifelse(is.na(cfqa5week), NA, cfqa5week))) %>% 
  mutate(cfqa5day  = ifelse(cfqa1==0, 0, ifelse(is.na(cfqa5day), ifelse(is.na(cfqa5week), NA, cfqa5week/7), cfqa5day)))%>% 
  mutate(cfb3      = ifelse(cfqb1==0, 0, ifelse(is.na(cfb3), NA, cfb3))) %>% 
  mutate(cfqb5week = ifelse(cfqb1==0, 0, ifelse(is.na(cfqb5week), NA, cfqb5week))) %>% 
  mutate(cfqb5day  = ifelse(cfqb1==0, 0, ifelse(is.na(cfqb5day), ifelse(is.na(cfqb5week), NA, cfqb5week/7), cfqb5day)))%>% 
  mutate(cfc3      = ifelse(cfqc1==0, 0, ifelse(is.na(cfc3), NA, cfc3))) %>% 
  mutate(cfqc5week = ifelse(cfqc1==0, 0, ifelse(is.na(cfqc5week), NA, cfqc5week))) %>% 
  mutate(cfqc5day  = ifelse(cfqc1==0, 0, ifelse(is.na(cfqc5day), ifelse(is.na(cfqc5week), NA, cfqc5week/7), cfqc5day)))%>% 
  mutate(ncups = select(., c(cfqa5day, cfqb5day, cfqc5day)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(Caffeine.PREDIGT = ifelse(is.na(ncups),NA, ifelse(ncups>=2, 2, ifelse(ncups>=1, 1, 0)))) %>% 
  mutate(Caffeine.PREDIGT.fct = factor(Caffeine.PREDIGT))
data_unchange <- data_unchange %>% 
  left_join(Caffeine %>% select(PATNO, Caffeine_infoyear, Caffeine.PREDIGT)) 

Head <- df_Head %>% 
  select(patno, head_injury_timestamp, 
         hiq1, hiq2, hiqa2, hiqb2,hiqc3unc, hiqc4unc, hiqc5unc) %>% 
  rename(PATNO=patno) %>% 
  separate(head_injury_timestamp,c("infomonth", "infoday","Head_infoyear","infohour","infomin"),remove = FALSE) 
Head[Head==9999] <- NA
Head <- Head %>% 
  mutate(concussion = select(., c(hiqa2:hiqc5unc)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(concussion = ifelse(is.na(hiq1), NA, concussion)) %>% 
  mutate(Head.PREDIGT = ifelse(is.na(hiq1),NA, ifelse(concussion>0, 2, ifelse(hiq1>0, 1, 0)))) %>% 
  mutate(Head.PREDIGT.fct = factor(Head.PREDIGT))
data_unchange <- data_unchange %>% 
  left_join(Head %>% select(PATNO, Head_infoyear, Head.PREDIGT)) 

Pstcdfarm <- df_Pesticide %>% 
  select(patno, pesticides_at_work_timestamp, 
         pwlabelintro1, phintro, rsqa5, rsqb5, rsqc5, rsqd5, rsqe5, rsqf5, rsqg5, 
         rsqa4, rsqb4) %>% 
  rename(PATNO=patno, pstcdfarm_timestamp=pesticides_at_work_timestamp) %>% 
  separate(pstcdfarm_timestamp,c("infomonth", "infoday","Pstcdfarm_infoyear","infohour","infomin"),remove = FALSE) 
Pstcdfarm[Pstcdfarm==9999] <- NA  
Pstcdfarm[Pstcdfarm==7777] <- NA  
Pstcdfarm <- Pstcdfarm %>% 
  mutate(pesticide = select(., c(pwlabelintro1:rsqg5)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(farm = select(., c(rsqa4:rsqb4)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(pest.farm = select(., c(pesticide:farm)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(Farm.PREDIGT = ifelse(pest.farm>0,1,0)) %>% 
  mutate(Farm.PREDIGT.fct = factor(Farm.PREDIGT))
data_unchange <- data_unchange %>% 
  left_join(Pstcdfarm %>% select(PATNO, Pstcdfarm_infoyear, Farm.PREDIGT)) 

Smoking <- df_Smoking %>% 
  select(patno, smoking_timestamp, 
         smq1, smq2, smq4, smq5, smq5a1_more, smq5a2_more, 
         smq3_year, smq4_year, smq3_age, smq4_age, 
         smq5astart1, smq5astop1, smq5astart2, smq5astop2, smq5astart3, smq5astop3) %>% 
  rename(PATNO=patno) %>% 
  separate(smoking_timestamp,c("infomonth", "infoday","Smoking_infoyear","infohour","infomin"),remove = FALSE) 
Smoking[Smoking==9999] <- NA  
Smoking[Smoking==7777] <- NA  
Smoking <- Smoking %>% 
  left_join(Demographics %>% select(PATNO, BIRTHDT)) 
Smoking <- Smoking %>% 
  mutate(age      = as.numeric(Smoking_infoyear) - BIRTHDT) %>% 
  mutate(smq3_age = ifelse(is.na(smq3_age),ifelse(is.na(smq3_year),NA, smq3_year-BIRTHDT),smq3_age)) %>% 
  mutate(smq4_age = ifelse(is.na(smq4_age),ifelse(is.na(smq4_year),NA, smq4_year-BIRTHDT),smq4_age)) %>% 
  mutate(smq4_age = ifelse(smq1+smq1>0, ifelse(is.na(smq4_age),age,smq4_age),NA)) %>% 
  mutate(smq2     = ifelse(smq1==0, 0, ifelse(is.na(smq2), NA, smq2))) %>% 
  mutate(smq4     = ifelse(smq1==0, 0, ifelse(is.na(smq4), NA, smq4))) %>% 
  mutate(smq5     = ifelse(smq1==0, 0, ifelse(is.na(smq5), NA, smq5))) %>% 
  mutate(smq3_age = -smq3_age) %>% 
  mutate(smq5astop1 = -smq5astop1) %>% 
  mutate(smq5astop2 = -smq5astop2) %>%
  mutate(smq5astop3 = -smq5astop3) %>%
  mutate(smk.dur = select(., c(smq3_age:smq5astop3)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(Smoking.PREDIGT = case_when(
    !is.na(smq4) & smq4 == 1 & smk.dur >= 20 ~ 12,
    !is.na(smq4) & smq4 == 1 & smk.dur >= 11 ~ 8,
    !is.na(smq4) & smq4 == 0 & smk.dur >= 20 ~ 4,
    !is.na(smq4) & smq4 == 0 & smk.dur >= 11 ~ 2,
    !is.na(smq1) & smq1 == 1                 ~ 1,
    TRUE ~ 0
  ))
Smoking[which(is.na(Smoking$smq1)),"Smoking.PREDIGT"] <- NA
data_unchange <- data_unchange %>% 
  left_join(Smoking %>% select(PATNO, Smoking_infoyear, Smoking.PREDIGT)) 

Toxicant <- df_Toxicant %>% 
  select(patno, toxicant_timestamp, 
         tx3, tx4, tx5) %>% 
  rename(PATNO=patno) %>% 
  separate(toxicant_timestamp,c("infomonth", "infoday","Toxicant_infoyear","infohour","infomin"),remove = FALSE) 
Toxicant[Toxicant==9999] <- NA  
Toxicant[Toxicant==7777] <- NA  
Toxicant <- Toxicant %>% 
  mutate(Metal.PREDIGT = select(., c(tx3:tx5)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(Metal.PREDIGT = ifelse(Metal.PREDIGT>0,1,0)) %>% 
  mutate(Metal.PREDIGT.fct = factor(Metal.PREDIGT))
data_unchange <- data_unchange %>% 
  left_join(Toxicant %>% select(PATNO, Toxicant_infoyear, Metal.PREDIGT)) 

# ------------ combine data frames, only select PD and HC ---------------
data_temporal %>% group_by(EVENT_ID) %>% summarise(count=n()) %>% arrange(desc(count)) %>% print(n = Inf)

data_join <- data_temporal %>% 
  filter(PATNO %in% data_unchange$PATNO) %>% 
  left_join(data_unchange)
data_join <- data_join %>% 
  rename(MDS.UPDRS.1.11=NP1CNST, MDS.UPDRS.1.03=NP1DPRS, MDS.UPDRS.1.04=NP1ANXS,
         SCOPA.AUT.05=SCAU5, SCOPA.AUT.06=SCAU6,
         group=group4)
if (clean)
{
  data_join <- data_join %>% filter(group2 == "Healthy Control" | group2 == "Parkinson disease") %>% droplevels()
} else
{
  data_join <- data_join %>% filter(group == "Healthy Control" | group == "Parkinson disease") %>% droplevels()
}
data_join[grep('_infoyear', names(data_join))] <- lapply(data_join[grep('_infoyear', names(data_join))],as.numeric)
data_join <- data_join %>% 
  mutate(infoyear = select(., c(GDS_infoyear,UPDRS_infoyear,RBDQ_infoyear,SCOPAAUT_infoyear,
                                STAI_infoyear,UPSIT_infoyear)) %>% apply(1, Mode)) %>% 
  mutate(age      = infoyear - BIRTHDT) %>% 
  mutate(EVENT_ID = factor(EVENT_ID)) %>% 
  select(-c(GDS_infoyear,UPDRS_infoyear,RBDQ_infoyear,SCOPAAUT_infoyear,STAI_infoyear,UPSIT_infoyear,
            Caffeine_infoyear,Head_infoyear,Pstcdfarm_infoyear,Smoking_infoyear,Toxicant_infoyear,
            group1,group2,group3))

EVENT_vec <- c("SC","BL",
               "V01","V02","V03","V04","V05","V06","V07","V08",
               "V09","V10","V11","V12","V13","V14","V15","V16",
               "PW","RS1","ST", "U01","other")
data_join <- data_join %>% 
  mutate(EVENT_ID = factor(EVENT_ID, levels=EVENT_vec)) %>% 
  arrange(EVENT_ID) %>% 
  arrange(PATNO)

EVENT_vec[2] <- "baseline"
levels(data_join$EVENT_ID)[2] <- "baseline"

levels(data_join$group) <- group_name
names(data_join)[1]     <- "ID"
names(data_join)[48]    <- "sex"

# ======= Tabel 1: Demographic summary ======
# Hoehn And Yahr Stage
df_HY <- df_motor %>% select("SUBJ_ID", "X3.21.MDS.UPDRS...Hoehn.And.Yahr.Stage..UPD2HY.", "VISIT_NAME")
names(df_HY)[1:2] <- c("ID", "HYstage")
df_HY <- df_HY %>% filter(VISIT_NAME == "Baseline")

found <- data_join %>% filter(ID %in% Caffeine$PATNO)
table1 <- data_join 
table1 <- found

data_T1 <- table1 %>% 
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
  c("MDS.UPDRS.1.11","SCOPA.AUT.05","SCOPA.AUT.06"),
  c("MDS.UPDRS.1.03", "GDS"),
  c("MDS.UPDRS.1.04", "STAIS"),
  c("RBD.SQ"),
  c("Alpha.Synuclein"),
  c("tTau"),
  c("GRS")
)
varStr2 <- c("Constipation", "Depression", "Anxiety", "RBD","Alpha.Synuclein","tTau","GRS")
EVENT_vec <- c("baseline","V04","V06","V08","V10","V12","V14")

df_LR <- varScore(data, varList, varStr2, EVENT_vec, group_name)

# ========== univariate results =========
df_LR <- df_LR %>% 
  mutate(RBD.SQ.06.1to4 = select(., c(DRMVERBL:DRMOBJFL)) %>% apply(1, sum, na.rm=TRUE)) %>% 
  mutate(RBD.SQ.06.1to2 = select(., c(DRMVERBL,DRMFIGHT)) %>% apply(1, sum, na.rm=TRUE)) 
varStr <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT",
            "Constipation.PREDIGT","MDS.UPDRS.1.11","SCOPA.AUT.05","SCOPA.AUT.06",
            "Smell.PREDIGT","UPSIT",
            "Smoking.PREDIGT", "Caffeine.PREDIGT",
            "FH.PREDIGT", "FH1","FH2","FH.any","GRS.PREDIGT","GRS",
            "Alpha.Synuclein","tTau",
            "Depression.PREDIGT", "MDS.UPDRS.1.03","GDS",
            "Anxiety.PREDIGT", "MDS.UPDRS.1.04", "STAIS","STAIT","STAI",
            "RBD.SQ","sex","age",
            "RBD.SQ.06.1to4","RBD.SQ.06.1to2")

#results <- uniVar(df_LR, varStr, EVENT_vec, group_name)

pltStr  <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT",
             "Smell.PREDIGT", "Smoking.PREDIGT", "Caffeine.PREDIGT",
             "FH.PREDIGT", "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.SQ")
pltStr2  <- c("E: Metal", "E: Pesticide", "E: Head Trauma", "E: Constipation",
              "E: Hyposmia", "E: Smoking", "E: Caffeine",
              "D: Family History", "I: Depression", "I: Anxiety", "I: RBD")
results <- uniVar(df_LR, pltStr, EVENT_vec, group_name, pltStr = pltStr, pltStr2 = pltStr2)
results[[1]]

# three cohort combined
resPlt <- read.csv(paste(wd,"/3-cohorts figure 2.csv", sep = ""), header = T, stringsAsFactors=FALSE, fileEncoding="latin1", na.strings=c("",".","NA"))


# ========== read in saved results of DeNoPa and combine results ===========
if (combine_DeNoPa)
{
  df_LR0  <- df_LR
  df_LR_P <- df_LR
  df_LR_D <- read.csv("/Users/juan/Desktop/OHRI/work/Cohorts/!DeNoPa PPMI Temporal/DeNoPa/DeNoPa dfLR.csv", 
                      header = T, stringsAsFactors=FALSE, fileEncoding="latin1", na.strings=c("",".","NA"))
  
  df_LR_P <- df_LR_P %>% 
    droplevels() %>% 
    mutate(
      EVENT_ID = as.character(EVENT_ID),
      cohort = "PPMI",
      ID = as.character(ID))
  df_LR_D <- df_LR_D %>% mutate(cohort = "DeNoPa")
  
  df_LR <- bind_rows(df_LR_P, df_LR_D) %>% 
    mutate(group = factor(group))
}
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
# ------- 1. PPMI original -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT", "Smell.PREDIGT",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", 
             "FH.PREDIGT","GRS.PREDIGT","Alpha.Synuclein.PREDIGT","tTau.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
tableStr <- "PPMI original"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 9:  genetic risk score: 0.25,0.5,0.75,1
           0.25,    # 10: alpha-synuclein: 0.25,0.5,1
           0.25,    # 11: T-tau: 0.25,0.5,1
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)    
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:7)
ind_D <- c(8,9)
ind_I <- c(10:14)

# ------- 2. PPMI original without FOUND-------
varStr2 <- c("Constipation.PREDIGT", "Smell.PREDIGT",
             "FH.PREDIGT","GRS.PREDIGT","Alpha.Synuclein.PREDIGT","tTau.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
tableStr <- "PPMI original without FOUND"
# Coeffcients
coefV <- c(0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           0.125,   # 8:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 9:  genetic risk score: 0.25,0.5,0.75,1
           0.25,    # 10: alpha-synuclein: 0.25,0.5,1
           0.25,    # 11: T-tau: 0.25,0.5,1
           0.25,    # 12: depression
           0.25,    # 13: anxiety
           0.25     # 14: RBD
)    
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1,2)
ind_D <- c(3,4)
ind_I <- c(5:9)

# ------- 3. PPMI original without CRF and GRS-------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "Constipation.PREDIGT", "Smell.PREDIGT",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
if (combine_DeNoPa) varStr1 <- append(varStr1, "cohort")
tableStr <- "PPMI original without CRF and GRS"
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

# ------- 4. PPMI original without CRF, GRS, FOUND-------
varStr2 <- c("Constipation.PREDIGT", "Smell.PREDIGT",
             "FH.PREDIGT",
             "Depression.PREDIGT", "Anxiety.PREDIGT", "RBD.PREDIGT")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
tableStr <- "PPMI original without CRF, GRS, FOUND"
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

# ------- 5. self-report -------
varStr2 <- c("Metal.PREDIGT", "Farm.PREDIGT", "Head.PREDIGT", "SCOPA.AUT.06",
             "Smoking.PREDIGT", "Caffeine.PREDIGT", "FH.PREDIGT",
             "MDS.UPDRS.1.03", "MDS.UPDRS.1.04","RBD.SQ.06.1to2")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
tableStr <- "PPMI self-report"
# Coeffcients
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

# ------- 6. self-report, no E -----------
varStr2 <- c("SCOPA.AUT.06",
             "FH.PREDIGT",
             "MDS.UPDRS.1.03", "MDS.UPDRS.1.04","RBD.SQ.06.1to2")
varStr1 <- append(c("ID", "group", "sex","age","EVENT_ID"), varStr2)
tableStr2 <- "PPMI self-report, no E"
# Coeffcients
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
if (combine_DeNoPa)
{
  scoreHC_D <- df %>% filter(cohort == "DeNoPa", group==group_name[1]) %>% select(score) 
  scorePD_D <- df %>% filter(cohort == "DeNoPa", group==group_name[2]) %>% select(score) 
  scoreHC_P <- df %>% filter(cohort == "PPMI", group==group_name[1]) %>% select(score) 
  scorePD_P <- df %>% filter(cohort == "PPMI", group==group_name[2]) %>% select(score) 
  print(paste("DeNoPa - HC: ", round(mean(scoreHC_D$score), 2), " (+/- ", round(sd(scoreHC_D$score), 2), ") ",
              "DeNoPa - PD: ", round(mean(scorePD_D$score), 2), " (+/- ", round(sd(scorePD_D$score), 2), ") ", sep=""))
  print(paste("PPMI - HC: ", round(mean(scoreHC_P$score), 2), " (+/- ", round(sd(scoreHC_P$score), 2), ") ",
              "PPMI - PD: ", round(mean(scorePD_P$score), 2), " (+/- ", round(sd(scorePD_P$score), 2), ") ", sep=""))
} 

# ====== Plot: distribution and box plot ==========
cutline <- TRUE
violin <- FALSE
adj <- 1.2
xmax <- 300
plotDistBox(df_LR, cut_BL, EVENT_vec, group_name, plot_path, tableStr, My_Theme, save_plot, xmax=xmax) 

# ====== 2-step ===============
step2 <- c("UPSIT", 
           "GRS",
           "GDS", 
           "STAIS",
           "RBD.SQ",
           "Alpha.Synuclein",
           "tTau")
Model2 <- model2(data_LR, step2, group_name)

# ====== 2-step: scatter plot ========
indTest  <- c(1, 3, 6)
model2point(Model2, df_LR, step2, indTest, plot_path, tableStr, My_Theme,
            save_plot)

# ====== 2-step: ROC plot =======
indTest  <- c(1, 3, 6)
model2ROC(Model2, df_LR, step2, indTest, plot_path, tableStr, My_Theme,
          save_plot=FALSE)

# ====== 2-step: boxplot =======
indTest  <- 1:length(step2) # index of step 2 tests that will be include into the boxplot
model2box(Model2, step2, indTest, plot_path, tableStr, My_Theme,
          save_plot=FALSE)

# ====== 2-step: bootstrapping =========
indTest  <- c(1, 3, 6)
reps  <- 1000
boot_results <- model2boot(df_LR, step2, indTest, reps)
save(boot_results,file="DeNoPa bootstrapping.Rdata")

load("DeNoPa bootstrapping.Rdata")
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
df_HY <- df_motor %>% select("SUBJ_ID", "X3.21.MDS.UPDRS...Hoehn.And.Yahr.Stage..UPD2HY.", "Summary.Score", "VISIT_NAME")
names(df_HY)[1:3] <- c("ID", "HYstage", "MDS_UPDRS_3")
df_HY <- df_HY %>% filter(VISIT_NAME == "Baseline")

df_HY <- df_LR %>% 
  filter(EVENT_ID == "baseline") %>% 
  select(ID, group, score) %>% 
  left_join(df_HY, by = "ID") %>% 
  filter(!is.na(score), !is.na(HYstage)) %>% 
  mutate(HYstage = factor(HYstage))

df_HY2 <- df_HY %>% 
  group_by(HYstage) %>% 
  summarise(n=n())
levels(df_HY$HYstage) <- paste(df_HY2$HYstage, " (n=", df_HY2$n, ")", sep="")

ttlStr <- "PPMI Model1'"
df_plt <- df_HY
pHC    <- TRUE

df_plt <- df_HY %>% filter(group == group_name[2]) %>% droplevels()
pHC    <- FALSE

motorPlt(df_plt, ttlStr, pHC)


