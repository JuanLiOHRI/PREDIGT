# Main funtion for PREDIGT analysis on DeNoPa
# Created by Juan Li, 2022-04-01

library(dplyr)
library(pROC)
library(ggplot2)
library(ggpubr)
library(rms)
library(pracma)
library(boot)

plt <- FALSE

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
data_path  <- "/Users/juan/Desktop/OHRI/work/Cohorts/C-OPN/data"
plot_path  <- "/Users/juan/Desktop/OHRI/work/Cohorts/C-OPN/Plots"
group_name <- c("Healthy Control", "Parkinson's Disease")
# ------------- Read data -----------
# read in all data
COPN  <- read.csv(paste(data_path, "/", "COPN_DATA_2021-10-20_1226.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
CaPRI <- read.csv(paste(data_path, "/", "CaPRI_DATA_2021-10-20_1234.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))
QPN   <- read.csv(paste(data_path, "/", "QPN Healthy Control Data_PREDIGT_February 10 2022_Final.csv",sep=""),header = T, fileEncoding="latin1", check.names=FALSE, na.strings=c("",".","NA"))

# ------------ data prepatation -------------------
# ----ID----
names(COPN)[1]  <- "ID"
names(CaPRI)[1] <- "ID"
names(QPN)[1]   <- "ID"
COPN   <- COPN %>% 
  filter(ID != "PD01483" & ID != "PD01738") # remove two participants that are also in QPN controls
QPN    <- QPN %>% 
  filter(ID != "PD01483" & ID != "PD01738") # remove two participants that are also in C-OPN

COPN <- COPN %>% 
  mutate(site = case_when(
    grepl('CA', ID)   ~ "Calgary",
    grepl('MNI', ID)  ~ "MNI",
    grepl('PD', ID)   ~ "QPN",
    grepl('TOH', ID)  ~ "Ottawa",
    grepl('UBC', ID)  ~ "Vancouver",
    grepl('UDM', ID)  ~ "UDM",
    grepl('UOA', ID)  ~ "Edmonton",
    grepl('UOT', ID)  ~ "Toronto"
  ))

ggplot(COPN, aes(site)) +
  geom_bar(stat = "count") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  ylim(0,250)
  

# ----Diagnosis----
group1 <- COPN %>% select("ID", contains("diagnosis"), contains("duration"), contains("first_symptoms")) %>% 
  mutate(cohort = "C-OPN",
         group  = group_name[2])
group2 <- CaPRI %>% select("ID") %>% 
  mutate(cohort = "CaPRI",
         group  = group_name[1])
group3 <- QPN %>% select("ID") %>% 
  mutate(cohort = "QPN",
         group  = group_name[1])
group  <- bind_rows(group1, group2, group3) %>% 
  mutate(diagnosis_pd = factor(diagnosis_pd),
         diagnosis_probability = factor(diagnosis_probability),
         diagnosis_confirmation = factor(diagnosis_confirmation),
         diagnosis_pdplus = factor(diagnosis_pdplus))
levels(group$diagnosis_pd) <- c("0, No", "1, Yes", "2, Uncertain")
levels(group$diagnosis_probability) <- c("1, >90% (high certainty)",
                                         "2, 50-89% (likely)",
                                         "4, Unknown", 
                                         "5, Not applicable")
levels(group$diagnosis_confirmation) <- c("1, Patient chart",
                                          "2, General neurologist",
                                          "3, Movement disorders neurologist",
                                          "4, Other")
levels(group$diagnosis_pdplus) <- c("1, Progressive Supranuclear Palsy (PSP)",
                                    "2, Multiple System Atrophy (MSA)",
                                    "3, Corticobasal Syndrome (CBS)",
                                    "8, Not Determined")

ggplot(group, aes(diagnosis_pd)) +
  geom_bar()+
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  ylim(0,650)
  
ggplot(group %>% filter(!is.na(diagnosis_pdplus)), aes(diagnosis_pdplus)) +
  geom_bar()+
  geom_text(stat='count', aes(label=..count..), hjust=-1)+
  coord_flip()+
  ylim(0,8)

group <- group %>% 
  mutate(group = ifelse((!is.na(diagnosis_pd) & (diagnosis_pd == "0, No" | diagnosis_pd == "2, Uncertain")),
                        as.character(diagnosis_pdplus), group)) %>% 
  mutate(group = ifelse(group==group_name[2] & is.na(diagnosis_pd), NA, group)) %>% 
  mutate(group = factor(group)) %>% 
  mutate(cohort = factor(cohort))
colSums(group[,15:20],na.rm=TRUE) # first symptoms
summary(group$diagnosis_pd)
summary(group %>% select(diagnosis_pdplus))
summary(group %>% filter(diagnosis_pd=="0, No") %>% select(diagnosis_pdplus))
summary(group %>% filter(diagnosis_pd=="2, Uncertain") %>% select(diagnosis_pdplus))

# ----Age----
names(QPN)[2] <- "study_visit_age"
age1 <- COPN %>% select("ID", contains("age")) 
age2 <- CaPRI %>% select("ID", contains("age")) 
age3 <- QPN %>% select("ID", "study_visit_age") 
age  <- bind_rows(age1, age2, age3)
age  <- age %>% 
  mutate(study_visit_age = ifelse(is.na(study_visit_age),ifelse(!is.na(age_ottawa), age_ottawa,study_visit_age),study_visit_age),
         age_onset = ifelse(is.na(age_onset),ifelse(!is.na(age_onset_4), age_onset_4,age_onset),age_onset),
         age_onset_2 = ifelse(is.na(age_onset_2),ifelse(!is.na(age_onset_3), age_onset_3,age_onset_2),age_onset_2)) %>% 
  select(ID, study_visit_age, age_onset, age_onset_2)

data <- group %>% select(ID, group, cohort, duration_disease, duration_disease_2, duration_first_symptoms) 
data <- data %>% 
  left_join(age, by="ID") 
# %>% 
#   filter(group %in% group_name) %>% droplevels() 
sum(!is.na(data %>% filter(cohort=="C-OPN", group %in% group_name) %>% select(study_visit_age)))
sum(!is.na(data %>% filter(cohort=="CaPRI") %>% select(study_visit_age)))
sum(!is.na(data %>% filter(cohort=="QPN") %>% select(study_visit_age)))

data <- data %>% 
  filter(!is.na(study_visit_age))
  
data <- data %>% 
  mutate(age_onset_cal = ifelse(!is.na(age_onset), age_onset, 
                                ifelse(!is.na(duration_disease_2), round(study_visit_age - duration_disease_2),
                                       ifelse(!is.na(duration_disease), round(study_visit_age - duration_disease),NA))),
         age_onset_2_cal = ifelse(!is.na(age_onset_2), age_onset_2, 
                                ifelse(!is.na(duration_first_symptoms), round(study_visit_age - duration_first_symptoms),NA))) %>% 
  mutate(age_onset_cal = ifelse(group == group_name[1], study_visit_age, age_onset_cal),
         age_onset_2_cal = ifelse(group == group_name[1], study_visit_age, age_onset_2_cal))

ggplot(data %>% filter(group %in% group_name), aes(duration_disease_2))+
  geom_histogram(breaks=c(0,5,10,15,20,40))+
  stat_bin(breaks=c(0,5,10,15,20,40),aes(y=..count.., label=..count..), geom="text", vjust=-.5)+
  xlab("Disease duration (years)")+
  ylim(0,310)

ggplot(data %>% filter(group %in% group_name), aes(study_visit_age,fill=group))+
  geom_histogram(alpha=0.5)+
  scale_fill_manual(values=c("Healthy Control"="blue",
                             "Parkinson's Disease"="red"))+
  theme(legend.position = "top")
data %>% filter(study_visit_age<=40) %>% select(ID, cohort, group, study_visit_age)

age_min <- 40
age_max <- 100
p1 <- ggplot(data %>% filter(group %in% group_name), aes(group, study_visit_age,fill = group))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  ylab("Age at study visit (years)")+ 
  scale_fill_manual(values=c("Healthy Control"="blue",
                             "Parkinson's Disease"="red"))+
  ylim(age_min,age_max)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     label.x = 1.35, label.y = age_max) # p-value
p2 <- ggplot(data %>% filter(group %in% group_name), aes(group, age_onset_cal,fill = group))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  ylab("Age at diagnosis of Parkinson's disease (years)")+ 
  scale_fill_manual(values=c("Healthy Control"="blue",
                             "Parkinson's Disease"="red"))+
  ylim(age_min,age_max)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     label.x = 1.35, label.y = age_max) # p-value
my_comparisons <- list( c("C-OPN", "CaPRI"), c("C-OPN", "QPN"), c("CaPRI", "QPN") )
p3 <- ggplot(data %>% filter(study_visit_age>=40, group %in% group_name), aes(cohort, study_visit_age,fill = cohort))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  ylab("Age (years)")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif")

ggarrange(p1,p3,p2, ncol=3,common.legend = TRUE)

data <- data %>% filter(study_visit_age >= 40)
nrow(data %>% filter(cohort=="C-OPN", group %in% group_name))
nrow(data %>% filter(cohort=="CaPRI"))
nrow(data %>% filter(cohort=="QPN"))

# ----Gender----
names(QPN)[3] <- "gender"
gender1 <- COPN %>% select(ID, gender) 
gender2 <- CaPRI %>% select(ID, gender) 
gender3 <- QPN %>% select(ID, gender) %>% 
  mutate(gender = ifelse(gender=="Femme",2,1))
gender  <- bind_rows(gender1, gender2, gender3) %>% 
  mutate(gender_fct = factor(gender))
levels(gender$gender_fct) <- c("Male", "Female")

data <- data %>% 
  left_join(gender, by="ID") %>% 
  filter(!is.na(gender))
nrow(data %>% filter(cohort=="C-OPN", group %in% group_name))
nrow(data %>% filter(cohort=="CaPRI"))
nrow(data %>% filter(cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "gender_fct", "Gender", or_cal = TRUE)
  plotBar(data %>% filter(cohort!="QPN",group %in% group_name), "gender_fct", "Gender", or_cal = TRUE)
  plotBar(data %>% filter(cohort!="CaPRI",group %in% group_name), "gender_fct", "Gender", or_cal = TRUE)
  plotBar(data %>% filter(cohort!="C-OPN") %>% droplevels(), "gender_fct", varlab = "Gender", col2 = "cohort")
}

# site
ind  <- data %>% filter(group == group_name[2]) %>% select(ID) %>% droplevels()
COPN_temp <- COPN %>% filter(ID %in% levels(ind$ID))
ggplot(COPN_temp, aes(site)) +
  geom_bar(stat = "count") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1)+
  ylim(0,250)
# ---- Table 1 ---------
d1 <- data %>% 
  group_by(group,gender_fct) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n))

d2 <- data %>% 
  filter(group %in% group_name) %>% 
  group_by(cohort,gender_fct) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n))

PD <- data %>% filter(group == group_name[2])
HC <- data %>% filter(group == group_name[1])
HC_C <- data %>% filter(cohort == "CaPRI")
HC_Q <- data %>% filter(cohort == "QPN")
quantile(PD$study_visit_age,na.rm = TRUE)
quantile(HC$study_visit_age,na.rm = TRUE)
quantile(HC_C$study_visit_age,na.rm = TRUE)
quantile(HC_Q$study_visit_age,na.rm = TRUE)

# ----Ethnicity----
race1 <- COPN %>% select("ID", contains("ethnicity")) 
race2 <- CaPRI %>% select("ID", contains("ethnicity")) 
race  <- bind_rows(race1, race2) 
race   <- race %>% 
  mutate(ethnicity___1  = ifelse(ethnicity_father___1 == 1 | ethnicity_mother___1 == 1, 1, 0),
         ethnicity___2  = ifelse(ethnicity_father___2 == 1 | ethnicity_mother___2 == 1, 1, 0),
         ethnicity___3  = ifelse(ethnicity_father___3 == 1 | ethnicity_mother___3 == 1, 1, 0),
         ethnicity___4  = ifelse(ethnicity_father___4 == 1 | ethnicity_mother___4 == 1, 1, 0),
         ethnicity___5  = ifelse(ethnicity_father___5 == 1 | ethnicity_mother___5 == 1, 1, 0),
         ethnicity___6  = ifelse(ethnicity_father___6 == 1 | ethnicity_mother___6 == 1, 1, 0),
         ethnicity___7  = ifelse(ethnicity_father___7 == 1 | ethnicity_mother___7 == 1, 1, 0),
         ethnicity___8  = ifelse(ethnicity_father___8 == 1 | ethnicity_mother___8 == 1, 1, 0),
         ethnicity___9  = ifelse(ethnicity_father___9 == 1 | ethnicity_mother___9 == 1, 1, 0),
         ethnicity___10 = ifelse(ethnicity_father___10 == 1 | ethnicity_mother___10 == 1, 1, 0),
         ethnicity___11 = ifelse(ethnicity_father___11 == 1 | ethnicity_mother___11 == 1, 1, 0),
         ethnicity___12 = ifelse(ethnicity_father___12 == 1 | ethnicity_mother___12 == 1, 1, 0),
         ethnicity___13 = ifelse(ethnicity_father___13 == 1 | ethnicity_mother___13 == 1, 1, 0),
         ethnicity___14 = ifelse(ethnicity_father___14 == 1 | ethnicity_mother___14 == 1, 1, 0),
         ethnicity___15 = ifelse(ethnicity_father___15 == 1 | ethnicity_mother___15 == 1, 1, 0),
         ethnicity___16 = ifelse(ethnicity_father___16 == 1 | ethnicity_mother___16 == 1, 1, 0))

multi1 <- rowSums(race[,2:17],na.rm=TRUE)  # father multiracial
multi2 <- rowSums(race[,19:34],na.rm=TRUE) # mother multiracial
multi  <- rowSums(race[,36:51],na.rm=TRUE) # child multiracial
race   <- race %>% 
  mutate(ethnicity_missing = ifelse(multi1 == 0 | multi2 == 0, 1, 0),
         ethnicity_multi = ifelse(multi > 1, 1, 0),
         ethnicity_count = multi) %>% 
  mutate(ethnicity = factor(case_when(
    ethnicity_multi == 1   ~ "17, Multiracial",
    ethnicity___1   == 1   ~ "01, Caucasian",
    ethnicity___2   == 1   ~ "02, French Canadian",
    ethnicity___3   == 1   ~ "03, First Nations",
    ethnicity___4   == 1   ~ "04, Hispanic or Latino",
    ethnicity___5   == 1   ~ "05, African or Caribbean or Afro American",
    ethnicity___6   == 1   ~ "06, North African",
    ethnicity___7   == 1   ~ "07, Middle Eastern",
    ethnicity___8   == 1   ~ "08, Indian",
    ethnicity___9   == 1   ~ "09, Chinese",
    ethnicity___10  == 1   ~ "10, Southeast Asian",
    ethnicity___11  == 1   ~ "11, Pacific Islander",
    ethnicity___12  == 1   ~ "12, Sephardic Jewish",
    ethnicity___13  == 1   ~ "13, Ashkenazi Jewish",
    ethnicity___14  == 1   ~ "14, Uncertain",
    ethnicity___15  == 1   ~ "15, Other",
    ethnicity___16  == 1   ~ "16, Unknown",
    TRUE ~ "18, Missing"
  ))) %>% 
  select(ID, ethnicity)

data <- data %>% 
  left_join(race, by="ID") 
nrow(data %>% filter(cohort=="C-OPN" & !is.na(ethnicity) & group %in% group_name))
nrow(data %>% filter(cohort=="CaPRI" & !is.na(ethnicity)))
nrow(data %>% filter(cohort=="QPN" & !is.na(ethnicity)))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "ethnicity", "Ethnicity", level_rvs = TRUE)
}

# ----Handedness----
hand1 <- COPN %>% select("ID", contains("handedness")) %>% 
  mutate(handedness_fct = factor(case_when(
    handedness == 1 ~ "Right hand",
    handedness == 2 ~ "Left hand",
    handedness == 3 ~ "Ambidextrous"
  ))) %>% 
  select(ID, handedness_fct)
hand2 <- CaPRI %>% select(ID, handedness, ehi_hand_class) %>% 
  mutate(handedness_fct = factor(case_when(
    ehi_hand_class == 3 ~ "Right hand",
    ehi_hand_class == 1 ~ "Left hand",
    ehi_hand_class == 2 ~ "Ambidextrous"
  )))%>% 
  select(ID, handedness_fct)
hand3 <- QPN %>% select("ID", "Lat\u008eralit\u008e") 
names(hand3)[2] <- "handedness"
hand3 <- hand3 %>% 
  mutate(handedness_fct = factor(case_when(
    handedness == "Droitier" ~ "Right hand",
    handedness == "Gaucher" ~ "Left hand",
    handedness == "Gaucher contrari\u008e" ~ "Left hand"
  )))%>% 
  select(ID, handedness_fct)

hand <- bind_rows(hand1, hand2, hand3)
data <- data %>% 
  left_join(hand, by="ID") 
nrow(data %>% filter(cohort=="C-OPN" & !is.na(handedness_fct)))
nrow(data %>% filter(cohort=="CaPRI" & !is.na(handedness_fct)))
nrow(data %>% filter(cohort=="QPN" & !is.na(handedness_fct)))

# Plot
if (plt)
{
  plotBar(data, "handedness_fct", "Handedness", or_cal = TRUE)
}

# ----Education----
edu1 <- COPN %>% select(ID, education) %>% 
  mutate(education_fct = factor(case_when(
    education == 1 ~ "1, None",
    education == 2 ~ "2, Grade School",
    education == 3 ~ "3, High School",
    education == 4 ~ "4, Trade Certificate or Diploma",
    education == 5 ~ "5, Non-University Diploma",
    education == 6 ~ "6, Bachelor's Degree",
    education == 7 ~ "7, Postgraduate Degree"
  ))) %>% 
  select(ID, education_fct)
edu2 <- CaPRI %>% select(ID, education) %>% 
  mutate(education_fct = factor(case_when(
    education == 1 ~ "1, None",
    education == 2 ~ "2, Grade School",
    education == 3 ~ "3, High School",
    education == 4 ~ "4, Trade Certificate or Diploma",
    education == 5 ~ "5, Non-University Diploma",
    education == 6 ~ "6, Bachelor's Degree",
    education == 7 ~ "7, Postgraduate Degree"
  ))) %>% 
  select(ID, education_fct)
edu3 <- QPN %>% select(ID, "Scolarit\u008e")
names(edu3)[2] <- "education"
edu3 <- edu3 %>% 
  mutate(education_fct2 = factor(case_when(
    education == "Primaire" ~ "1, Primaire",
    education == "Universitaire" ~ "4, Universitaire",
    education == "Secondaire" ~ "2, Secondaire",
    education == "Coll\u008egiale" ~ "3, Collégiale"
  ))) %>% 
  select(ID, education_fct2)
edu <- bind_rows(edu1,edu2,edu3)

data <- data %>% left_join(edu, by="ID")
nrow(data %>% filter(cohort=="C-OPN" & !is.na(education_fct)))
nrow(data %>% filter(cohort=="CaPRI" & !is.na(education_fct)))
nrow(data %>% filter(cohort=="QPN" & !is.na(education_fct2)))

# Plot
if (plt)
{
  plotBar(data, "education_fct", "Education", level_rvs = TRUE)
  plotBar(data, "education_fct2", "Education", level_rvs = TRUE)
  plotBar(data %>% filter(gender==1), "education_fct", "Education", level_rvs = TRUE) #male
  plotBar(data %>% filter(gender==2), "education_fct", "Education", level_rvs = TRUE) #female
  plotBar(data %>% filter(gender==1), "education_fct2", "Education", level_rvs = TRUE)
  plotBar(data %>% filter(gender==2), "education_fct2", "Education", level_rvs = TRUE)
}

# ----Comorbidities----
com1 <- COPN %>% select(ID, contains("comorbidities"))
com2 <- CaPRI %>% select(ID, contains("comorbidities"))
com3 <- QPN %>% select(ID, "Hypertension", "Hypotension", "Diab\u008fte")
names(com3)[2:4] <- c("comorbidities___8", "comorbidities___9",
                        "comorbidities___10")
com3 <- com3 %>% 
  mutate(comorbidities___8 = ifelse(!is.na(comorbidities___8),
                                    ifelse(comorbidities___8 == "Non",0,1),NA),
         comorbidities___9 = ifelse(!is.na(comorbidities___9),
                                    ifelse(comorbidities___9 == "Non",0,1),NA),
         comorbidities___10 = ifelse(!is.na(comorbidities___10),
                                     ifelse(comorbidities___10 == "Non",0,1),NA))
com  <- bind_rows(com1, com2, com3)

data <- data %>% left_join(com, by="ID") 

levelStr <- c("01, Cancer",
              "02, Myocardial infarction or angina (heart attack)",
              "03, Transient ischemic attack (TIA)",
              "04, Respiratory problems",
              "05, Lung disease",
              "06, Cerebrovascular accident (stroke)",
              "07, Anemia or other blood disease",
              "08, Hypertension (high blood pressure)",
              "09, Hypotension (low blood pressure)",
              "10, Diabetes",
              "11, Hypercholesterolemia (high cholesterol)",
              "12, Ulcers or stomach disease",
              "13, Inflammatory bowel disease (IBD)",
              "14, Liver disease",
              "15, Kidney disease",
              "16, Osteoarthritis (degenerative)",
              "17, Rheumatoid Arthritis (autoimmune)",
              "18, Other",
              "19, None of these")

d1 <- data %>% select(cohort, contains("comorbidities")) %>% filter(cohort=="C-OPN")
d2 <- data %>% select(cohort, contains("comorbidities")) %>% filter(cohort=="CaPRI")
d3 <- data %>% select(cohort, contains("comorbidities")) %>% filter(cohort=="QPN")

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "comorbidities", "Comorbidities", or_cal = TRUE, level_rvs = TRUE, 
          sel_contains = TRUE, levelStr = levelStr, ylim = c(-0.5,0.5))
}

# ----Metal Pesticide Farm----
exp1 <- COPN %>% select(ID, welding, pesticides, living_environment, 
                        work_environment, contains("urban"), contains("rural"))
exp2 <- CaPRI %>% select(ID, welding, pesticides, living_environment, 
                         work_environment, contains("urban"), contains("rural"))
exp3 <- QPN %>% select(ID, "Milieu Vie", "Milieu Travail")
names(exp3)[2:3] <- c("living_environment", "work_environment")
exp3 <- exp3 %>% 
  mutate(living_environment = ifelse(!is.na(living_environment), 
                                     ifelse(living_environment == "Urbain", 1, 2), NA),
         work_environment   = ifelse(!is.na(work_environment), 
                                     ifelse(work_environment == "Urbain", 1, 2), NA))
exp  <- bind_rows(exp1, exp2, exp3)
data <- data %>% left_join(exp, by="ID") 
names(data)[42:51] <- paste("env_",names(data)[42:51],sep="")

nrow(data %>% filter(!is.na(welding) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(welding) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(pesticides) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(pesticides) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(living_environment) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(living_environment) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(living_environment) & cohort=="QPN"))
nrow(data %>% filter(!is.na(work_environment) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(work_environment) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(work_environment) & cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "welding", "Welding", or_cal = TRUE, level_rvs = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(group %in% group_name), "pesticides", "Pesticides", or_cal = TRUE, level_rvs = TRUE, levelStr = c("No", "Yes"))
  
  plotBar(data %>% filter(group %in% group_name), "living_environment", "Living environment", or_cal = TRUE, level_rvs = TRUE, 
          levelStr = c("Urban" , "Rural"))
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "living_environment", "Living environment", 
          or_cal = TRUE, level_rvs = TRUE, levelStr = c("Urban" , "Rural"))
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "living_environment", "Living environment", 
          or_cal = TRUE, level_rvs = TRUE, levelStr = c("Urban" , "Rural"))
  
  plotBar(data %>% filter(group %in% group_name), "work_environment", "Work environment", or_cal = TRUE, level_rvs = TRUE, 
          levelStr = c("Urban" , "Rural"))
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "work_environment", "Work environment", 
          or_cal = TRUE, level_rvs = TRUE, levelStr = c("Urban" , "Rural"))
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "work_environment", "Work environment", 
          or_cal = TRUE, level_rvs = TRUE, levelStr = c("Urban" , "Rural"))
  
  
  plotBar(data %>% filter(cohort != "C-OPN") %>% droplevels(), "living_environment", "Living environment", 
          col2 = "cohort", level_rvs = TRUE, levelStr = c("Urban" , "Rural"))
  plotBar(data %>% filter(cohort != "C-OPN") %>% droplevels(), "work_environment", "Work environment", 
          col2 = "cohort", level_rvs = TRUE, levelStr = c("Urban" , "Rural"))

# Plot (set up manually)
levelStr <- c("1.1. living: Urban - near factory",
              "1.2. living: Urban - near farm",
              "1.3. living: Rural - near factory",
              "1.4. living: Rural - near farm",
              "1.5. living: Rural - on farm",
              "2.1. working: Urban - near factory",
              "2.2. working: Urban - near farm",
              "2.3. working: Rural - near factory",
              "2.4. working: Rural - near farm",
              "2.5. working: Rural - on farm")
plotBar(data %>% filter(group %in% group_name), "env_", "Environment", or_cal = TRUE, level_rvs = TRUE, 
        sel_contains = TRUE, levelStr = levelStr, ylim = c(-0.5,0.5))
}

# ----Head trauma----
ht1 <- COPN %>% select(ID, head_blow, head_blow_number, contact_sports) %>% 
  mutate(head_blow0 = case_when(
    head_blow == 0  ~ 0,
    head_blow_number == 1  ~ 1,
    head_blow_number > 1  ~ 2
  ))
ht2 <- CaPRI %>% select(ID, head_blow, head_blow_number, contact_sports)%>% 
  mutate(head_blow0 = case_when(
    head_blow == 0  ~ 0,
    head_blow_number == 1  ~ 1,
    head_blow_number > 1  ~ 2
  ))
ht3 <- QPN %>% select(ID, "Choc \u0088 la t\u0090te")
names(ht3)[2] <- c("head_blow0")
ht3 <- ht3 %>% 
  mutate(head_blow0 = case_when(
    head_blow0 == "Non"  ~ 0,
    head_blow0 == "Une seule fois"  ~ 1,
    head_blow0 == "Ë quelques reprises"  ~ 2
  ),
  head_blow = ifelse(is.na(head_blow0), NA,
                     ifelse(head_blow0>0, 1, 0)))
ht  <- bind_rows(ht1, ht2, ht3)
data <- data %>% left_join(ht, by="ID") 
nrow(data %>% filter(!is.na(head_blow) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(head_blow) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(contact_sports) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(contact_sports) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(head_blow) & cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "head_blow", "Head blow", or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE,
          levelStr = c("No", "Yes", "Uncertain"))
  plotBar(data %>% filter(group %in% group_name), "contact_sports", "Contact sports", or_cal = TRUE, level_rvs = TRUE,
          levelStr = c("No", "Yes"))
}

# group "Uncertain" as "No"
data <- data %>% 
  mutate(head_blow = ifelse(is.na(head_blow), NA,
                            ifelse(head_blow == 1, 1, 0)))
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "head_blow", "Head blow", or_cal = TRUE, level_rvs = TRUE,
          levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort!="QPN", group %in% group_name), "head_blow", "Head blow", or_cal = TRUE, level_rvs = TRUE,
          levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort!="CaPRI", group %in% group_name), "head_blow", "Head blow", or_cal = TRUE, level_rvs = TRUE,
          levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "C-OPN"), "head_blow", "Head blow", col2 = "cohort", level_rvs = TRUE,
          levelStr = c("No", "Yes"))
}

# ----Constipation----
cstp1 <- COPN %>% select(ID, contains("constipation")) 
cstp2 <- CaPRI %>% select(ID, contains("constipation"))
cstp3 <- QPN %>% select(ID, Constipation) %>% 
  mutate(Constipation = ifelse(is.na(Constipation), NA,
                               ifelse(Constipation == "Non", 0, 1)))
names(cstp3)[2] <- "constipation"
cstp  <- bind_rows(cstp1, cstp2, cstp3)
data <- data %>% left_join(cstp, by="ID") 
nrow(data %>% filter(!is.na(constipation) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(constipation) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(constipation) & cohort=="QPN"))
nrow(data %>% filter(!is.na(constipation_yrs) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(constipation_yrs) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(constipation_yrs) & cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "constipation", "Constipation", or_cal = TRUE, level_rvs = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "constipation", "Constipation", or_cal = TRUE, level_rvs = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "constipation", "Constipation", or_cal = TRUE, level_rvs = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "C-OPN"), "constipation", "Constipation", col2 = "cohort", level_rvs = TRUE, levelStr = c("No", "Yes"))
  
  ggplot(data, aes(group, constipation_yrs, fill=group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("Duration of constipation (years)")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"))+
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = max(data$constipation_yrs, na.rm=TRUE)) 
}

# ----Hyposmia----
sml1 <- COPN %>% select(ID, smell) %>% 
  mutate(hyposmia = case_when(
    smell == 0 ~ 0,
    smell == 1 | smell == 3 ~ 1,
    smell == 2 ~ 2
  ))
sml2 <- CaPRI %>% select(ID, smell) %>% 
  mutate(hyposmia = case_when(
    smell == 0 ~ 0,
    smell == 1 | smell == 3 ~ 1,
    smell == 2 ~ 2
  ))
sml <- bind_rows(sml1, sml2)
data <- data %>% left_join(sml, by="ID") 
nrow(data %>% filter(!is.na(hyposmia) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(hyposmia) & cohort=="CaPRI"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "hyposmia", "Hyposmia", or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE, 
          levelStr = c("No", "Yes/Already lost", "Uncertain"))
  plotBar(data %>% filter(group %in% group_name), "smell", "Hyposmia", or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE, 
          levelStr = c("No", "Yes", "Uncertain", "Already lost"))
}

# group "Uncertain" as "No"
data <- data %>% 
  mutate(smell = ifelse(is.na(smell), NA,
                        ifelse(smell == 2, 0, smell)))
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "smell", "Hyposmia", or_cal = TRUE, level_rvs = TRUE,or_combine = TRUE, 
          levelStr = c("No", "Yes", "Already lost"))
}

# ----Drugs----
drug1 <- COPN %>% select(ID, contains("drugs"))
drug2 <- CaPRI %>% select(ID, contains("drugs"))
drug  <- bind_rows(drug1, drug2) %>% select(-recreational_drugs)

data <- data %>% left_join(drug, by="ID") 

levelStr <- c("1, Alcohol",
              "2, Cigarettes",
              "3, Other tobacco products",
              "4, Recreational drugs",
              "5, Cannabis",
              "6, None of the above")

nrow(data %>% filter(!is.na(drugs___1) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(drugs___1) & cohort=="CaPRI"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "drugs__", "Drugs", or_cal = TRUE, level_rvs = TRUE, 
          sel_contains = TRUE, levelStr = levelStr)
}

# ----Smoking----
smk1 <- COPN %>% select(ID, drugs___2, drugs___3, smoking, smoking_yrs, cigarettes, 
                        past_smoking_years, past_cigarettes) %>% 
  mutate(smoke = case_when(
    drugs___2 == 1 | drugs___3 == 1 ~ 2,
    smoking == 1 ~ 1,
    drugs___2 == 0 & drugs___3 == 0 & smoking == 0 ~ 0
  ),
  smoke_yrs = ifelse(!is.na(smoking_yrs), smoking_yrs,
                     ifelse(!is.na(past_smoking_years), past_smoking_years, NA)),
  smoke_n = ifelse(!is.na(cigarettes), cigarettes,
                   ifelse(!is.na(past_cigarettes), past_cigarettes, NA)))
smk2 <- CaPRI %>% select(ID, drugs___2, drugs___3, smoking, smoking_yrs, cigarettes, 
                        past_smoking_years, past_cigarettes) %>% 
  mutate(smoke = case_when(
    drugs___2 == 1 | drugs___3 == 1 ~ 2,
    smoking == 1 ~ 1,
    drugs___2 == 0 & drugs___3 == 0 & smoking == 0 ~ 0
  ),
  smoke_yrs = ifelse(!is.na(smoking_yrs), smoking_yrs,
                     ifelse(!is.na(past_smoking_years), past_smoking_years, NA)),
  smoke_n = ifelse(!is.na(cigarettes), cigarettes,
                   ifelse(!is.na(past_cigarettes), past_cigarettes, NA)))
smk3 <- QPN %>% select(ID, "Fumeur") %>% 
  mutate(smoke = case_when(
    Fumeur == "Oui" ~ 2,
    Fumeur == "Jadis" ~ 1,
    Fumeur == "Non" ~ 0
  ))
smk <- bind_rows(smk1, smk2, smk3) %>% 
  mutate(smoke_PREDIGT = case_when(
    !is.na(smoke) & smoke == 2 & smoke_yrs >= 20 ~ 12,
    !is.na(smoke) & smoke == 2 & smoke_yrs >= 11 ~ 8,
    !is.na(smoke) & smoke == 1 & smoke_yrs >= 20 ~ 4,
    !is.na(smoke) & smoke == 1 & smoke_yrs >= 11 ~ 2,
    !is.na(smoke) & smoke > 0                  ~ 1,
    !is.na(smoke) & smoke == 0                 ~ 0
  ))
data <- data %>% left_join(smk, by="ID") 
nrow(data %>% filter(!is.na(smoke_PREDIGT) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(smoke_PREDIGT) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(smoke_PREDIGT) & cohort=="QPN"))
nrow(data %>% filter(!is.na(smoke_yrs) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(smoke_yrs) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(smoke_yrs) & cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "smoke", "Smoking", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "smoke", "Smoking", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "smoke", "Smoking", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "C-OPN"), "smoke", "Smoking", levelStr = c("No", "Past", "Current"), 
          col2 = "cohort", level_rvs = TRUE)
  
  plotBar(data %>% filter(cohort!="QPN", group %in% group_name), "smoke_PREDIGT", "Smoking", 
          levelStr = c("No", "Any smoking <=10 yrs", "Past 11-19 yrs", "Past >=20 yrs",
                       "Current 11-19 yrs", "Current >=20 yrs"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  
  ggplot(data, aes(group, smoke_yrs, fill=group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("Duration of smoking (years)")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"))+
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = max(data$smoke_yrs, na.rm = TRUE))
}

# ----Caffeine----
cafe1 <- COPN %>% select(ID, contains("coffee")) %>% 
  mutate(cafe = case_when(
    coffee      == 1 ~ 2,
    coffee_past == 1 ~ 1,
    coffee      == 0 ~ 0
  ),
  cafe_yrs = ifelse(!is.na(coffee_yrs), coffee_yrs,
                     ifelse(!is.na(coffee_past_yrs), coffee_past_yrs, NA)),
  cafe_n = ifelse(!is.na(coffee_cups), coffee_cups,
                   ifelse(!is.na(coffee_past_cups), coffee_past_cups, NA))) %>% 
  select(ID, cafe, cafe_yrs, cafe_n)
cafe2 <- CaPRI %>% select(ID, contains("coffee")) %>% 
  mutate(cafe = case_when(
    coffee      == 1 ~ 2,
    coffee_past == 1 ~ 1,
    coffee      == 0 ~ 0
  ),
  cafe_yrs = ifelse(!is.na(coffee_yrs), coffee_yrs,
                    ifelse(!is.na(coffee_past_yrs), coffee_past_yrs, NA)),
  cafe_n = ifelse(!is.na(coffee_cups), coffee_cups,
                  ifelse(!is.na(coffee_past_cups), coffee_past_cups, NA))) %>% 
  select(ID, cafe, cafe_yrs, cafe_n)
cafe3 <- QPN %>% select(ID, "Caf\u008e")
names(cafe3)[2] <- "coffee"
cafe3 <- cafe3 %>% 
  mutate(cafe = case_when(
    coffee == "Oui" ~ 2,
    coffee == "Jadis" ~ 1,
    coffee == "Non" ~ 0
  )) %>% 
  select(ID, cafe)
cafe <- bind_rows(cafe1, cafe2, cafe3) %>% 
  mutate(cafe_yrs = ifelse(is.na(cafe_yrs),NA,
                           ifelse(cafe_yrs>100, 2022-cafe_yrs, cafe_yrs)),
         cafe_n   = ifelse(is.na(cafe_n),NA,
                           ifelse(cafe_n>100, cafe_n/365, cafe_n)))%>% 
  mutate(cafe_PREDIGT = case_when(
    !is.na(cafe) & cafe == 2 & cafe_n >= 2 ~ 2,
    !is.na(cafe) & cafe == 2 & cafe_n >= 1 ~ 1,
    !is.na(cafe)                           ~ 0
  ))
data <- data %>% left_join(cafe, by="ID") 
nrow(data %>% filter(!is.na(cafe_PREDIGT) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(cafe_PREDIGT) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(cafe_PREDIGT) & cohort=="QPN"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "cafe", "Coffee", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort!="QPN", group %in% group_name), "cafe", "Coffee", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort!="CaPRI", group %in% group_name), "cafe", "Coffee", levelStr = c("No", "Past", "Current"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  
  plotBar(data %>% filter(cohort!="C-OPN"), "cafe", "Coffee", col2 = "cohort", levelStr = c("No", "Past", "Current"), 
          level_rvs = TRUE)
  
  plotBar(data %>% filter(cohort!="QPN", group %in% group_name), "cafe_PREDIGT", "Coffee", levelStr = c("Other", "Current >=1 cup/day", "Current >= 2 cups/day"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  
  p1 <- ggplot(data %>% filter(group %in% group_name), aes(group, cafe_yrs, fill=group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("Duration of drinking coffee (years)")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"))+
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = max(data$cafe_yrs, na.rm = TRUE))
  
  p2 <- ggplot(data %>% filter(group %in% group_name), aes(group, cafe_n, fill=group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("Cups of coffee per day")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"))+
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = max(data$cafe_n, na.rm = TRUE))
  ggarrange(p1, p2, ncol=2, common.legend = TRUE)
}

# ----Excercise----
exc1 <- COPN %>% select(ID, contains("exercise")) %>% 
  select(-c(qpn_exercise_intensity, qpn_exercise_intensity_2, qpn_exercise_intensity_3))
exc2 <- CaPRI %>% select(ID, contains("exercise")) 
exc  <- bind_rows(exc1, exc2)
data <- data %>% left_join(exc, by="ID") %>% select(-exercise_other)
data <- data %>% 
  mutate(exercise_duration = ifelse(!is.na(exercise_duration), exercise_duration,
                                    ifelse(exercise == 0, 0, NA)))
nrow(data %>% filter(!is.na(exercise_duration) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(exercise_duration) & cohort=="CaPRI"))

# Plot
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "exercise", "Exercise", levelStr = c("No", "Yes"), 
          or_cal = TRUE, level_rvs = TRUE)
  
  plotBar(data %>% filter(group %in% group_name), "exercise_duration", "Exercise duration", 
          levelStr = c("0, No","1, < 30 min", "2, 30-60 min", "3, 1-2 hours", "4, >2 hours"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  
  levelStr <- c("01, Walking or hiking",
                "02, Running",
                "03, Cycling",
                "04, Swimming",
                "05, Yoga",
                "06, Pilates",
                "07, Dancing",
                "08, Weight Training",
                "09, Golf",
                "10, Skiing or X-Country Skiing",
                "11, Hockey",
                "12, Parkinson Wellness Recovery (PWR!Moves)",
                "13, Boxing",
                "99, Other")
  plotBar(data %>% filter(group %in% group_name), "exercise_type", "Exercise type", or_cal = TRUE, level_rvs = TRUE, 
          sel_contains = TRUE, levelStr = levelStr)
}

# ----Family History----
fh1 <- COPN %>% select(ID, father, mother, grandparents, contains("_pd")) %>% select(-c("diagnosis_pd", "diagnosis_pdplus"))
fh2 <- CaPRI %>% select(ID, father, mother, grandparents, contains("_pd")) %>% select(-c("diagnosis_pd", "diagnosis_pdplus"))
fhC <- bind_rows(fh1, fh2) %>% 
  mutate(
    father_fh    = father,
    mother_fh    = mother,
    parent_fh    = ifelse(father_fh == 1 | mother_fh == 1, 1, 
                          ifelse(is.na(father_fh) | is.na(mother_fh), NA, 0)),
    grandparents_fh = grandparents,
    siblings_fh  = ifelse(siblings_pd___1 == 1 | siblings_pd___2 == 1, 1, 0),
    brothers_fh = ifelse(siblings_pd___1 == 1, 1, 0),
    sisters_fh  = ifelse(siblings_pd___2 == 1, 1, 0),
    children_fh  = ifelse(children_pd___1 == 1 | children_pd___2 == 1, 1, 0),
    sons_fh = ifelse(children_pd___1 == 1, 1, 0),
    daughters_fh  = ifelse(children_pd___2 == 1, 1, 0),
    aunts_uncles_fh  = ifelse(aunts_uncles_pd___1 == 1 | aunts_uncles_pd___2 == 1, 1, 0),
    maternal_auntuncle_fh = ifelse(aunts_uncles_pd___1 == 1, 1, 0),
    paternal_auntuncle_fh  = ifelse(aunts_uncles_pd___2 == 1, 1, 0),
    cousins_fh  = ifelse(cousins_pd___1 == 1 | cousins_pd___2 == 1, 1, 0),
    maternal_cousins_fh = ifelse(cousins_pd___1 == 1, 1, 0),
    paternal_cousins_fh  = ifelse(cousins_pd___2 == 1, 1, 0),
    nephews_nieces_fh  = ifelse(nephews_nieces_pd___1 == 1 | nephews_nieces_pd___2 == 1, 1, 0),
    nephews_fh = ifelse(nephews_nieces_pd___1 == 1, 1, 0),
    nieces_fh  = ifelse(nephews_nieces_pd___2 == 1, 1, 0)
  )
names(fhC)[19:22] <- c("grandparents_fh___1", "grandparents_fh___2", 
                       "grandparents_fh___3", "grandparents_fh___4")
fhC <- fhC %>% select(ID, contains("_fh"))
fh3 <- QPN %>% select(ID, "Histoire familiale") 
names(fh3)[2] <- "FH"
fh3 <- fh3 %>% 
  mutate(father_fh = ifelse(is.na(FH), 0, ifelse(FH == "P\u008fre MP", 1, 0)),
         mother_fh = ifelse(is.na(FH), 0, ifelse(FH == "M\u008fre MP", 1, 0)),
         parent_fh = ifelse(is.na(FH), 0, ifelse(FH == "P\u008fre MP" | FH == "M\u008fre MP" | FH == "Parent\u008e MP", 1, 0)),
         siblings_fh = ifelse(is.na(FH), 0, ifelse(FH == "Fratrie MP", 1, 0))) %>% 
  select(-FH)
fh <- bind_rows(fhC, fh3)

data <- data %>% left_join(fh, by="ID") %>% 
  mutate(
    fh1 = case_when(
      parent_fh == 1 | children_fh == 1 | siblings_fh == 1 ~ 1,
      parent_fh == 0 & siblings_fh == 0 ~ 0
    ),
    fh2 = ifelse(grandparents_fh == 1 | aunts_uncles_fh == 1 | nephews_nieces_fh == 1, 1, 0),
    fh3 = cousins_fh,
    fh.PREDIGT  = case_when(
      fh1 == 1 ~ 4,
      fh2 == 1 ~ 2,
      fh3 == 1 ~ 1,
      TRUE     ~ 0
    )
  ) %>% 
  mutate(fh.PREDIGT = ifelse(fh.PREDIGT>0, fh.PREDIGT,
                             ifelse(is.na(fh1)&is.na(fh2),NA,0)))
nrow(data %>% filter(!is.na(fh1) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(fh1) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(fh1) & cohort=="QPN"))
nrow(data %>% filter(!is.na(fh2) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(fh2) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(fh2) & cohort=="QPN"))
nrow(data %>% filter(!is.na(fh3) & cohort=="C-OPN" & group %in% group_name))
nrow(data %>% filter(!is.na(fh3) & cohort=="CaPRI"))
nrow(data %>% filter(!is.na(fh3) & cohort=="QPN"))

if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "fh.PREDIGT", "Family History", levelStr = c("No", "3rd", "2nd", "1st"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "fh.PREDIGT", "Family History", levelStr = c("No", "3rd", "2nd", "1st"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "fh.PREDIGT", "Family History", levelStr = c("No", "3rd", "2nd", "1st"), 
          or_cal = TRUE, or_combine = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "C-OPN"), "fh.PREDIGT", "Family History", levelStr = c("No", "3rd", "2nd", "1st"), 
          col2 = "cohort", level_rvs = TRUE)
}

# FH
df1 <- data %>% 
  select(group, fh1, fh2, fh3) 
varlab   <- "Family history"
levelStr <- c("1st", "2nd", "3rd")
ylim  <- c(-0.5, 0.5)
or_cal <- TRUE
sel_contains <- TRUE

# Parent (1st)
df1 <- data %>% 
  select(group, parent_fh, father_fh, mother_fh) 
varlab   <- "Parent"
levelStr <- c("Parent", "Father", "Mother")
ylim  <- c(-0.5, 0.5)

# Children (1st)
df1 <- data %>% 
  select(group, children_fh, sons_fh, daughters_fh) 
varlab   <- "Children"
levelStr <- c("Children", "Son", "Daughter")
ylim  <- c(-0.5, 0.5)

# Siblings (1st)
df1 <- data %>% 
  select(group, siblings_fh, brothers_fh, sisters_fh) 
varlab   <- "Siblings"
levelStr <- c("Siblings", "Brother", "Sister")
ylim  <- c(-0.5, 0.5)

# Grandparent (2nd)
df1 <- data %>% 
  select(group, contains("grandparents_fh")) 
varlab   <- "Grandparent"
levelStr <- c("Paternal grandfather", "Paternal grandmother", 
              "Maternal grandfather", "Maternal grandmother",
              "Grandparent")
ylim  <- c(-0.5, 0.5)

# Aunt & Uncle (2nd)
df1 <- data %>% 
  select(group, aunts_uncles_fh, maternal_auntuncle_fh, paternal_auntuncle_fh) 
varlab   <- "Aunt & Uncle"
levelStr <- c("Aunt & Uncle", "Maternal", "Paternal")
ylim  <- c(-0.5, 0.5)

# Nephews & Nieces (2nd)
df1 <- data %>% 
  select(group, nephews_nieces_fh, nephews_fh, nieces_fh) 
varlab   <- "Nephews & Nieces"
levelStr <- c("Nephews & Nieces", "Nephews", "Nieces")
ylim  <- c(-0.5, 0.5)

# Cousins (3rd)
df1 <- data %>% 
  select(group, cousins_fh, maternal_cousins_fh, paternal_cousins_fh) 
varlab   <- "Cousins"
levelStr <- c("Cousins", "Maternal", "Paternal")
ylim  <- c(-0.5, 0.5)

# ----Psychiatric----
psych1 <- COPN %>% select(ID, contains("psychiatric"))
psych2 <- CaPRI %>% select(ID, contains("psychiatric"))
psych3 <- QPN %>% select(ID, "Anxi\u008et\u008e", "D\u008epression", "Maladie bipolaire",
                         "Hallucinations", "Impulsions/compulsions",
                         "D\u008esordres psychiatriques Autre")
names(psych3)[2:7] <- c("psychiatric_disorder___1", "psychiatric_disorder___2",
                   "psychiatric_disorder___6", "psychiatric_disorder___5",
                   "Impulsions", "psychiatric_disorder___8")
psych3 <- psych3 %>% 
  mutate(psychiatric_disorder___1 = ifelse(!is.na(psychiatric_disorder___1),
                                           ifelse(psychiatric_disorder___1 == "Non",0,1),NA),
         psychiatric_disorder___2 = ifelse(!is.na(psychiatric_disorder___2),
                                           ifelse(psychiatric_disorder___2 == "Non",0,1),NA),
         psychiatric_disorder___6 = ifelse(!is.na(psychiatric_disorder___6),
                                           ifelse(psychiatric_disorder___6 == "Non",0,1),NA),
         psychiatric_disorder___5 = ifelse(!is.na(psychiatric_disorder___5),
                                           ifelse(psychiatric_disorder___5 == "Non",0,1),NA),
         psychiatric_disorder___8 = ifelse(!is.na(psychiatric_disorder___8),
                                           ifelse(psychiatric_disorder___8 == "Non",0,1),NA))
psychi  <- bind_rows(psych1, psych2, psych3) %>% select(-Impulsions)

data <- data %>% left_join(psychi, by="ID") 

if (plt)
{
  levelStr <- c("0, None",
                "1, Anxiety",
                "2, Depression",
                "3, Parkinson's Disease Psychosis",
                "4, Apathy",
                "5, Hallucinations",
                "6, Bipolar Disorder",
                "7, Primary Psychosis (e.g. Schizophrenia)",
                "8, Other")
  d1 <- data %>% select(cohort, contains("psychiatric_disorder")) %>% filter(cohort=="C-OPN")
  d2 <- data %>% select(cohort, contains("psychiatric_disorder")) %>% filter(cohort=="CaPRI")
  d3 <- data %>% select(cohort, contains("psychiatric_disorder")) %>% filter(cohort=="QPN")
  
  # Plot
  plotBar(data %>% filter(group %in% group_name), "psychiatric_disorder", "Psychiatric disorder", level_rvs = TRUE, or_cal = TRUE, 
          sel_contains = TRUE, levelStr = levelStr)
}

# ----Depression----
dps1 <- COPN %>% select(ID, psychiatric_disorder___2, bdsii_total_score) %>% 
  mutate(depression = psychiatric_disorder___2) %>% 
  mutate(bdsii_class = case_when(
    !is.na(bdsii_total_score) & bdsii_total_score <=13   ~ 1,
    !is.na(bdsii_total_score) & bdsii_total_score <=19   ~ 2,
    !is.na(bdsii_total_score) & bdsii_total_score <=28   ~ 3,
    !is.na(bdsii_total_score) & bdsii_total_score <=63   ~ 4
  ))
# Plot
if (plt)
{
  p <- plotBar(dps1, "bdsii_class", "Depression", col2 = "psychiatric_disorder___2", 
               levelStr = c("Minimal", "Mild", "Moderate", "Severe"))
  p <- p + 
    scale_fill_hue(name = "Depression",
                   labels = c("No","Yes"))
  print(p)
}

# logistic to find a better cut: bdsii_class? bdsii_total_score?
# datai <- dps1 %>% select(psychiatric_disorder___2, bdsii_total_score)
# roc_1 <- roc(datai$psychiatric_disorder___2, datai$bdsii_total_score, levels=c(0,1), na.rm=TRUE,
#              ci=TRUE)
# roc_1$auc # 0.6561
# 
# datai <- dps1 %>% select(psychiatric_disorder___2, bdsii_class)
# roc_1 <- roc(datai$psychiatric_disorder___2, datai$bdsii_class, levels=c(0,1), na.rm=TRUE,
#              ci=TRUE)
# roc_1$auc # 0.6407

dps2 <- CaPRI %>% select(ID, bdsii_total_score) %>% 
  mutate(bdsii_class = case_when(
    !is.na(bdsii_total_score) & bdsii_total_score <=13   ~ 1,
    !is.na(bdsii_total_score) & bdsii_total_score <=19   ~ 2,
    !is.na(bdsii_total_score) & bdsii_total_score <=28   ~ 3,
    !is.na(bdsii_total_score) & bdsii_total_score <=63   ~ 4
  )) %>% 
  mutate(depression = ifelse(is.na(bdsii_class),NA,ifelse(bdsii_class>1,1,0)))
dps3 <- QPN %>% select(ID, "D\u008epression") %>% 
  rename("Depression" = "D\u008epression") %>% 
  mutate(depression = ifelse(!is.na(Depression),ifelse(Depression == "Oui", 1, 0), NA)) %>% 
  select(ID, depression)
dps  <- bind_rows(dps1, dps2, dps3) %>% select(-psychiatric_disorder___2)
data <- data %>% left_join(dps, by="ID") %>% 
  mutate(depression_fct = factor(depression))
levels(data$depression_fct) <- c("No", "Yes")
nrow(dps1 %>% filter(!is.na(bdsii_class)))
nrow(dps2 %>% filter(!is.na(bdsii_class)))

# violin plot
if (plt)
{
  min <- min(data$bdsii_total_score, na.rm=TRUE)
  max <- max(data$bdsii_total_score, na.rm=TRUE)
  nHC <- nrow(data %>% filter(group==group_name[1] & !is.na(bdsii_total_score)))
  nPD <- nrow(data %>% filter(group==group_name[2] & !is.na(bdsii_total_score)))
  ggplot(data %>% filter(group %in% group_name), aes(group, bdsii_total_score,fill = group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("BDSII total score")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"),
                      label=c(paste("Healthy Control (n=",nHC,")",sep=""),
                              paste("Parkinson's Disease (n=",nPD,")",sep="")))+
    ylim(min,max)+
    geom_hline(yintercept = c(13, 19,28), color="red")+
    annotate("text", x = 1.46, y = c(9, 15, 24, 33), label = c("Minimal", "Mild", "Moderate", "Severe")) +
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = max) # p-value
  
  # Plot
  plotBar(data %>% filter(group %in% group_name), "bdsii_class", "Depression", or_cal = TRUE, or_combine = TRUE,level_rvs = TRUE, 
          levelStr = c("Minimal", "Mild", "Moderate", "Severe"))
  plotBar(data %>% filter(group %in% group_name), "depression_fct", "Depression", level_rvs = TRUE,or_cal = TRUE)
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "depression_fct", "Depression", level_rvs = TRUE,or_cal = TRUE)
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "depression_fct", "Depression", level_rvs = TRUE,or_cal = TRUE)
  plotBar(data %>% filter(cohort != "C-OPN"), "depression_fct", "Depression", level_rvs = TRUE,col2 = "cohort")
}

# ----Anxiety----
axt1 <- COPN %>% select(ID, psychiatric_disorder___1, bai_total_score, updrs_1_4) %>% 
  mutate(anxiety = psychiatric_disorder___1) %>% 
  mutate(bai_class = case_when(
    !is.na(bai_total_score) & bai_total_score <=7   ~ 1,
    !is.na(bai_total_score) & bai_total_score <=15   ~ 2,
    !is.na(bai_total_score) & bai_total_score <=25   ~ 3,
    !is.na(bai_total_score) & bai_total_score <=63   ~ 4
  ))
axt2 <- CaPRI %>% select(ID, bai_total_score)%>% 
  mutate(bai_class = case_when(
    !is.na(bai_total_score) & bai_total_score <=7   ~ 1,
    !is.na(bai_total_score) & bai_total_score <=15   ~ 2,
    !is.na(bai_total_score) & bai_total_score <=25   ~ 3,
    !is.na(bai_total_score) & bai_total_score <=63   ~ 4
  )) %>% 
  mutate(anxiety = ifelse(is.na(bai_class),NA,ifelse(bai_class>1,1,0)))
axt3 <- QPN %>% select(ID, "Anxi\u008et\u008e") %>% 
  rename("Anxiety" = "Anxi\u008et\u008e") %>% 
  mutate(anxiety = ifelse(!is.na(Anxiety),ifelse(Anxiety == "Oui", 1, 0), NA))%>% 
  select(ID, anxiety)
axt  <- bind_rows(axt1, axt2, axt3) %>% select(-psychiatric_disorder___1)
data <- data %>% left_join(axt, by="ID") %>% 
  mutate(anxiety_fct = factor(anxiety))
levels(data$anxiety_fct) <- c("No", "Yes")
nrow(axt1 %>% filter(!is.na(bai_class)))
nrow(axt2 %>% filter(!is.na(bai_class)))

# violin plot
if (plt)
{
  min <- min(data$bai_total_score, na.rm=TRUE)
  max <- max(data$bai_total_score, na.rm=TRUE)
  nHC <- nrow(data %>% filter(group==group_name[1] & !is.na(bai_total_score)))
  nPD <- nrow(data %>% filter(group==group_name[2] & !is.na(bai_total_score)))
  ggplot(data %>% filter(group %in% group_name), aes(group, bai_total_score,fill = group))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
    ylab("BAI total score")+ 
    scale_fill_manual(values=c("Healthy Control"="blue",
                               "Parkinson's Disease"="red"),
                      label=c(paste("Healthy Control (n=",nHC,")",sep=""),
                              paste("Parkinson's Disease (n=",nPD,")",sep="")))+
    ylim(min,20)+
    geom_hline(yintercept = c(7, 15, 25), color="red")+
    annotate("text", x = 1.46, y = c(6, 14, 17, 26), label = c("Minimal", "Mild", "Moderate", "Severe")) +
    stat_compare_means(label = "p.signif", method = "t.test",
                       label.x = 1.35, label.y = 20) # p-value
  
  # Plot
  plotBar(data %>% filter(group %in% group_name), "bai_class", "Anxiety")
  plotBar(data %>% filter(group %in% group_name), "anxiety_fct", "Anxiety", or_cal = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "anxiety_fct", "Anxiety", or_cal = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "anxiety_fct", "Anxiety", or_cal = TRUE, level_rvs = TRUE)
  plotBar(data %>% filter(cohort != "C-OPN"), "anxiety_fct", "Anxiety", col2 = "cohort", level_rvs = TRUE)
}

# ----RBD----
rbd1 <- COPN %>% select(ID, dreams, qpn_dreams) %>% 
  mutate(qpn_dreams = ifelse(is.na(qpn_dreams),NA,ifelse(qpn_dreams>1000, 2021-qpn_dreams, qpn_dreams)))
rbd2 <- CaPRI %>% select(ID, dreams, dreams_yrs) 
rbdC <- bind_rows(rbd1,rbd2) %>% 
  mutate(rbd_fct = factor(dreams))
levels(rbdC$rbd_fct) <- c("0, No", "1, Yes", "2, Uncertain")
rbd3 <- QPN %>% select(ID, "RBD") %>% 
  mutate(rbd_fct = factor(RBD))
levels(rbd3$rbd_fct) <- c("0, No", "1, Yes")
rbd <- bind_rows(rbdC, rbd3) %>% select(ID, rbd_fct)
data <- data %>% left_join(rbd, by="ID") %>% 
  mutate(rbd = case_when(
    rbd_fct == "0, No" ~ 0,
    rbd_fct == "1, Yes" ~ 1,
    rbd_fct == "2, Uncertain" ~ 2
  ))
nrow(data %>% filter(cohort=="C-OPN" & !is.na(rbd_fct) & group %in% group_name))
nrow(data %>% filter(cohort=="CaPRI" & !is.na(rbd_fct)))
nrow(data %>% filter(cohort=="QPN" & !is.na(rbd_fct)))

if (plt)
{
  ggplot(rbd1, aes(qpn_dreams))+
    geom_histogram()+
    xlab("RBD (years)")
  
  # Plot
  plotBar(data %>% filter(group %in% group_name), "rbd_fct", "RBD", level_rvs = TRUE, or_cal = TRUE, or_combine = TRUE,
          levelStr = c("No", "Yes", "Uncertain"))
}
# group "Uncertain" as "No"
data <- data %>% 
  mutate(rbd = ifelse(is.na(rbd), NA,
                      ifelse(rbd == 1, 1, 0)))
if (plt)
{
  plotBar(data %>% filter(group %in% group_name), "rbd", "RBD", level_rvs = TRUE, or_cal = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "QPN", group %in% group_name), "rbd", "RBD", level_rvs = TRUE, or_cal = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "CaPRI", group %in% group_name), "rbd", "RBD", level_rvs = TRUE, or_cal = TRUE, levelStr = c("No", "Yes"))
  plotBar(data %>% filter(cohort != "C-OPN"), "rbd", "RBD", level_rvs = TRUE, col2 = "cohort", levelStr = c("No", "Yes"))
  
  data <- data %>% left_join(COPN %>% select(ID, site), by = "ID")
  d <- data %>% 
    filter(group == group_name[2]) %>% 
    select(group, cohort, site, rbd) %>% 
    group_by(site) %>% 
    summarise(n1 = n(),
              n2 = sum(rbd, na.rm = TRUE),
              prop = n2/n1)
  
  ggplot(d, aes(x=site,y=prop))+
    geom_bar(stat='identity',position = "dodge", colour = "black")+
    geom_text(aes(x = site, y = prop, label = n2), 
              position = position_dodge(width = 1), vjust = -1) +
    ylim(0,0.5)
}

# ========== slice and prepare data ==========
df_0 <- data %>%
  select(ID, group, cohort, duration_disease, study_visit_age, age_onset_cal,
         gender_fct, welding, pesticides, head_blow,
         constipation, smell, smoke, smoke_PREDIGT, cafe, cafe_PREDIGT,
         exercise, fh1, fh.PREDIGT, depression, anxiety, rbd) %>%
  mutate(smell.PREDIGT = case_when(
           smell == 0             ~ 0,
           smell == 1| smell == 2 ~ 1,
           smell == 3             ~ 2
         ),
         smoke = ifelse(is.na(smoke), NA,
                        ifelse(smoke == 0, 0, 
                               ifelse(smoke == 1, 2, 8))), # similar scale as smoke_PREDIGT
         cafe  = ifelse(is.na(cafe), NA,
                        ifelse(cafe==2, 1, 0)), # combine "Past" with "No"
         fh1   = fh1 * 4
         )

names(df_0)[c(5,7)] <- c("age", "sex")
#df_0 <- df_0 %>% mutate(sex = ifelse(sex == "Female", 1, 2))

# ========== calculate PREDIGT Score =========
FG  <- TRUE # including factor G?
FT  <- TRUE # including factor T?
EVENT_vec <- c("BL")
df_LR <- df_0 %>% 
  mutate(
    score=0,
    EVENT_ID = EVENT_vec[1])

pltStr  <- c("welding", "pesticides", "head_blow", "constipation",
             "smell", "smoke", "cafe",
             "fh1", "depression", "anxiety", "rbd")
pltStr2  <- c("E: Metal", "E: Pesticide", "E: Head Trauma", "E: Constipation",
              "E: Hyposmia", "E: Smoking", "E: Caffeine",
              "D: Family History", "I: Depression", "I: Anxiety", "I: RBD")
results <- uniVar(df_LR %>% filter(group %in% group_name) %>% droplevels(), pltStr, EVENT_vec, group_name, pltStr = pltStr, pltStr2 = pltStr2)
results[[1]]

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
# ------- 1. C-OPN vs CaPRI+QPN -------
varStr2 <- c("head_blow", "constipation",
             "smoke", "cafe", "fh1",
             "depression", "anxiety", "rbd")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "Both"
# Coeffcients
coefV <- c(0.5,     # 1:  head trauma: no concussion*1, concussion*2
           0.5,     # 2:  constipation: 0.25, 0.5, 1
           -0.0625, # 3:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 4:  caffeine: >=1*1, >=2*2
           0.125,   # 5:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 6: depression
           0.25,    # 7: anxiety
           0.25     # 8: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:4)
ind_D <- 5
ind_I <- c(6:8)

# ------- 2. C-OPN vs CaPRI -------
varStr2 <- c("welding", "pesticides", "head_blow", "constipation", "smell.PREDIGT",
             "smoke_PREDIGT", "cafe_PREDIGT", "exercise", "fh.PREDIGT",
             "depression", "anxiety", "rbd")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "CaPRI"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  head trauma: no concussion*1, concussion*2
           0.5,     # 4:  constipation: 0.25, 0.5, 1
           0.5,     # 5:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.0625, # 6:  smoking: any<10y*1, Ex11-19*2, Ex20*4, Current11-19*8, Current20*12
           -0.125,  # 7:  caffeine: >=1*1, >=2*2
           -0.125,  # 8:  exercise
           0.125,   # 9:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 10: depression
           0.25,    # 11: anxiety
           0.25     # 12: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:8)
ind_D <- 9
ind_I <- c(10:12)

# ------- 3. C-OPN vs CaPRI+QPN No Smoking and Anxiety-------
varStr2 <- c("head_blow", "constipation",
             "cafe", "fh.PREDIGT",
             "depression", "rbd")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "CaPRI QPN-3"
# Coeffcients
coefV <- c(0.5,     # 1:  head trauma: no concussion*1, concussion*2
           0.5,     # 2:  constipation: 0.25, 0.5, 1
           -0.125,  # 3:  caffeine: >=1*1, >=2*2
           0.125,   # 4:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 5: depression
           0.25     # 6: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:3)
ind_D <- 4
ind_I <- c(5:6)

# ------- 4. C-OPN vs CaPRI No Smoking, Coffee, head-------
varStr2 <- c("welding", "pesticides", "constipation", "smell.PREDIGT",
             "exercise", "fh.PREDIGT",
             "depression", "anxiety", "rbd")
varStr1 <- append(c("ID", "group", "sex","age"), varStr2)
tableStr <- "CaPRI"
# Coeffcients
coefV <- c(0.5,     # 1:  Metal
           0.25,    # 2:  Pesiticide
           0.5,     # 3:  constipation: 0.25, 0.5, 1
           0.5,     # 4:  olfaction: hyposmia*1, anosmia*2. self-report*1
           -0.125,  # 5:  exercise
           0.125,   # 6:  family history: 3rd*1, 2nd*2, 1st*4
           0.25,    # 7: depression
           0.25,    # 8: anxiety
           0.25     # 9: RBD
)  
lenV <- length(coefV)
ind_question <- NULL # Only in self-report
ind_E <- c(1:5)
ind_D <- 6
ind_I <- c(7:9)

# ====== Calculate PREDIGT Score ===========
df_LR <- PREDIGTScr(df_LR, varStr1, varStr2, coefV, ind_question, EVENT_vec, FG, FT, method = 1)

# ====== Plot: ROC ==============
save_plot <- TRUE

df_ROC    <- cbind(df_LR,df_0$cohort) %>% filter(group %in% group_name)
df_CaPRI  <- df_ROC %>% filter(cohort!="QPN")
df_QPN    <- df_ROC %>% filter(cohort!="CaPRI")

dfPlot    <- df_ROC

ROC_result <- plotROC(dfPlot, EVENT_vec, group_name, plot_path, tableStr, My_Theme, save_plot)
cut_BL <- ROC_result$cut_BL
print(ROC_result$result)

df      <- dfPlot %>% filter(!is.na(score)) %>% filter(EVENT_ID=="BL")
scoreHC <- df %>% filter(group==group_name[1]) %>% select(score) 
scorePD <- df %>% filter(group==group_name[2]) %>% select(score) 
print(paste("HC: ", round(mean(scoreHC$score), 2), " (+/- ", round(sd(scoreHC$score), 2), ") ",
            "PD: ", round(mean(scorePD$score), 2), " (+/- ", round(sd(scorePD$score), 2), ") ", sep=""))

nCOPN <- nrow(df_ROC %>% filter(cohort == "C-OPN", !is.na(score)))
nCaPRI <- nrow(df_ROC %>% filter(cohort == "CaPRI", !is.na(score)))
nQPN <- nrow(df_ROC %>% filter(cohort == "QPN", !is.na(score)))

df_OND <- df_LR %>% filter(!(group %in% group_name), !is.na(group)) %>% droplevels()
nPSP <- nrow(df_OND %>% filter(group == "1, Progressive Supranuclear Palsy (PSP)"))
nMSA <- nrow(df_OND %>% filter(group == "2, Multiple System Atrophy (MSA)"))
nCBS <- nrow(df_OND %>% filter(group == "3, Corticobasal Syndrome (CBS)"))
nND  <- nrow(df_OND %>% filter(group == "8, Not Determined"))

dfOND_res <- df_OND %>% select(all_of(varStr1), score) %>% arrange(group)
  
  
ggplot(df_ROC, aes(score, fill = cohort)) + 
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = c(paste("C-OPN (n=", nCOPN, ")", sep=""),
                                 paste("CaPRI (n=", nCaPRI, ")", sep=""),
                                 paste("QPN (n=", nQPN, ")", sep=""))) 
  
my_comparisons <- list( c("C-OPN", "CaPRI"), c("C-OPN", "QPN"), c("CaPRI", "QPN") )
ggplot(df_ROC, aes(cohort, score,fill = cohort))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  ylab("PREDIGT Score")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test",label = "p.signif")

df_plt <- df_LR %>% filter(!is.na(group)) %>% select(group, score)
nHC    <- nrow(df_plt %>% filter(group == group_name[1], !is.na(score)))
nPD    <- nrow(df_plt %>% filter(group == group_name[2], !is.na(score)))
levels(df_plt$group) <- c(paste("PSP (n=", nPSP, ")", sep=""),
                          paste("MSA (n=", nMSA, ")", sep=""),
                          paste("CBS (n=", nCBS, ")", sep=""),
                          paste("Not determined (n=", nND, ")", sep=""),
                          paste("HC (n=", nHC, ")", sep=""),
                          paste("PD (n=", nPD, ")", sep=""))
ggplot(df_plt, aes(group, score,fill = group))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  ylab("PREDIGT Score")
  
# ====== Plot: distribution and box plot ==========
cutline <- FALSE
violin <- FALSE
adj <- 1.2
xmax <- 250
DB <- plotDistBox(dfPlot, cut_BL, EVENT_vec, group_name, plot_path, tableStr, My_Theme, 
                  save_plot, xmax = xmax, cutline= cutline) 

dp <- DB$dp
dp <- dp + 
  geom_vline(data = df_OND, aes(xintercept = score, colour = group)) +
  scale_color_discrete(labels = c(paste("PSP (n=", nPSP, ")", sep=""),
                                  paste("MSA (n=", nMSA, ")", sep=""),
                                  paste("CBS (n=", nCBS, ")", sep=""),
                                  paste("Not determined (n=", nND, ")", sep="")))
print(dp)






