# 2-step model
model2 <- function(data_LR, step2, group_name, sex_ind = 0, plot = FALSE) {
  ntest <- length(step2)
  
  if (sex_ind == 0) # Both
  {
    df0 <- data_LR 
  } else if (sex_ind == 1) # Female
  {
    df0 <- data_LR %>% filter(sex == "Female") 
  } else # Male
  {
    df0 <- data_LR %>% filter(sex == "Male") 
  }
  
  df0 <- df0 %>% 
    select(ID,EVENT_ID, group, score,names(data_LR[which(names(data_LR) %in% step2)])) %>% 
    filter(!is.na(score)) %>% 
    filter(EVENT_ID == "baseline") %>% 
    droplevels()
  nPD0 <- nrow(df0 %>% filter(group == group_name[2]))
  nHC0 <- nrow(df0 %>% filter(group == group_name[1]))
  
  step1List <- list() # length = ntest, roc objs of step 1 questionnaire of each df1 defined later
  step2List <- list() # list of list, roc objs of step 2 assessment with each C1
  combList  <- list() # list of list, length = ntest, results of combined ROC (gray zone)
  resultList <- list()
  
  for (itest in 1:ntest)
  {
    col <- step2[itest]
    df1 <- df0 %>% select(.dots = c("ID", "group", "score", all_of(col))) %>% 
      rename(ID=.dots1, group=.dots2, score = .dots3, col=.dots4) %>% 
      filter(!is.na(col))
    nPD1 <- nrow(df1 %>% filter(group == group_name[2]))
    nHC1 <- nrow(df1 %>% filter(group == group_name[1]))
    
    # step 1
    roc_1 <- roc(df1$group, df1$score, levels=group_name, na.rm=TRUE, ci=TRUE, quiet=TRUE)
    cut1  <- roc_1$thresholds
    ncut1 <- length(cut1)
    step1List[[itest]] <- roc_1
    
    # step 2
    rocList     <- list() # list of step 2 roc objects
    combineList <- list() # list of step 2 sensibility and specificity sections (gray zone)
    nPD2        <- rep(NA, ncut1)
    nHC2        <- rep(NA, ncut1)
    auc_combine <- rep(NA, ncut1) # combined AUC
    cut_opt_combine  <- rep(NA, ncut1) # optimal cut-value of step2 based on Youden Index
    sens_opt_combine <- rep(NA, ncut1) 
    spec_opt_combine <- rep(NA, ncut1)
    cut_opt_combine_PPV  <- rep(NA, ncut1) # optimal cut-value of step2 based on PPV
    sens_opt_combine_PPV <- rep(NA, ncut1)
    spec_opt_combine_PPV <- rep(NA, ncut1)
    PPV_opt_combine_PPV <- rep(NA, ncut1)
    NPV_opt_combine_PPV <- rep(NA, ncut1)
    
    for (icut1 in 1:ncut1)
    {
      # ROC of step 2
      df2 <- df1 %>% 
        filter(score > cut1[icut1]) %>% 
        droplevels()
      nPD2[icut1] <- nrow(df2 %>% filter(group==group_name[2]))
      nHC2[icut1] <- nrow(df2 %>% filter(group==group_name[1]))
      if (sum(df2$group==group_name[1])>=3)
      {
        roc_2 <- roc(df2$group, df2$col, levels=group_name,na.rm=TRUE,
                     ci=TRUE, quiet=TRUE)
        rocList[[icut1]] <- roc_2
        
        if (icut1 == 1) # step 2 test only
        {
          sens_combine <- roc_2$sensitivities
          spec_combine <- roc_2$specificities
          sensC <- roc_2$sensitivities
          specC <- roc_2$specificities
          
          # Combined ROC
          auc_combine[icut1]  <- roc_2$auc
        } else
        {
          sens_combine <- roc_2$sensitivities*roc_1$sensitivities[icut1]
          spec_combine <- 1-(1-roc_2$specificities)*(1-roc_1$specificities[icut1])
          sensC <- c(roc_1$sensitivities[1:icut1], sens_combine)
          specC <- c(roc_1$specificities[1:icut1], spec_combine)
          
          # Combined ROC
          auc_combine[icut1]  <- max(-cumtrapz(1-specC, sensC),na.rm=TRUE)
        }
        combineList[[icut1]] <- list(sens = sens_combine, spec = spec_combine)
        
        # find 'optimal' - Youden Index
        sum   <- sens_combine + spec_combine
        ind2  <- which(sum==max(sum,na.rm=TRUE))
        cut_opt  <- roc_2$thresholds[ind2]
        sens_opt <- sens_combine[ind2]
        spec_opt <- spec_combine[ind2]
        if (length(cut_opt)>1)
        {
          cut_opt  <- cut_opt[1]
          sens_opt <- sens_opt[1]
          spec_opt <- spec_opt[1]
        }
        cut_opt_combine[icut1]  <- cut_opt
        sens_opt_combine[icut1] <- sens_opt
        spec_opt_combine[icut1] <- spec_opt
        
        # find 'optimal' - Populational PPV and NPV
        PPV_combine <- sens_combine * prevalence     / 
          (sens_combine * prevalence + (1-spec_combine)*(1-prevalence)) #population PPV
        NPV_combine <- spec_combine * (1-prevalence) / 
          ((1 - sens_combine) * prevalence + spec_combine*(1-prevalence)) #population NPV
        
        sum_PPV      <- PPV_combine + NPV_combine
        ind2_PPV  <- which(sum_PPV==max(sum_PPV,na.rm=TRUE))
        cut_opt_PPV  <- roc_2$thresholds[ind2_PPV]
        sens_opt_PPV <- sens_combine[ind2_PPV]
        spec_opt_PPV <- spec_combine[ind2_PPV]
        PPV_opt_PPV <- PPV_combine[ind2_PPV]
        NPV_opt_PPV <- NPV_combine[ind2_PPV]
        if (length(cut_opt_PPV)>1)
        {
          cut_opt_PPV  <- cut_opt_PPV[1]
          sens_opt_PPV <- sens_opt_PPV[1]
          spec_opt_PPV <- spec_opt_PPV[1]
          PPV_opt_PPV  <- PPV_opt_PPV[1]
          NPV_opt_PPV  <- NPV_opt_PPV[1]
        }
        cut_opt_combine_PPV[icut1]  <- cut_opt_PPV
        sens_opt_combine_PPV[icut1] <- sens_opt_PPV
        spec_opt_combine_PPV[icut1] <- spec_opt_PPV
        PPV_opt_combine_PPV[icut1] <- PPV_opt_PPV
        NPV_opt_combine_PPV[icut1] <- NPV_opt_PPV
      }
    }
    
    step2List[[itest]] <- rocList
    combList[[itest]]  <- combineList
    resultList[[itest]] <- data.frame(nPD2  = nPD2,
                                      nHC2  = nHC2,
                                      cut1  = cut1,
                                      auc_combine = auc_combine,
                                      cut_opt_combine = cut_opt_combine,
                                      sens_opt_combine = sens_opt_combine,
                                      spec_opt_combine = spec_opt_combine,
                                      cut_opt_combine_PPV = cut_opt_combine_PPV,
                                      sens_opt_combine_PPV = sens_opt_combine_PPV,
                                      spec_opt_combine_PPV = spec_opt_combine_PPV,
                                      PPV_opt_combine_PPV = PPV_opt_combine_PPV,
                                      NPV_opt_combine_PPV = NPV_opt_combine_PPV)
    
    if (plot){
      plotdata <- resultList[[itest]]
      ncut     <- nrow(plotdata)
      plotdf <- data.frame(cut = rep(plotdata$cut1,5),
                           value = c(plotdata$auc_combine, 
                                     plotdata$sens_opt_combine, plotdata$spec_opt_combine, 
                                     plotdata$PPV_opt_combine_PPV, plotdata$NPV_opt_combine_PPV),
                           labels = c(rep("Combined AUC", ncut),
                                      rep("Sensitivity", ncut),
                                      rep("Specificity", ncut),
                                      rep("Pop PPV", ncut),
                                      rep("Pop NPV", ncut)))
      p <- ggplot(plotdf, aes(cut, value, color=factor(labels)))+
        geom_path()+
        labs(title = col)
      print(p)
    }
  }
  return(list(step1List  = step1List, step2List  = step2List,
              combList   = combList, resultList = resultList))
} 
