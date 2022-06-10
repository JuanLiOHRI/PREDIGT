# prepare for bootstrapping
model2boot <- function(data_LR, step2, indTest, reps = 1000) {
  ntest      <- length(indTest)
  step2vec   <- step2[indTest]
  
  df0 <- data_LR %>% 
    filter(EVENT_ID == "baseline") %>% 
    select(ID, group, score,names(data_LR[which(names(data_LR) %in% step2)])) %>% 
    filter(!is.na(score)) %>% 
    droplevels()
  
  boot_results <- list()
  for (istep in 1:ntest)
  {
    print(paste("istep:", istep))
    col <- step2vec[istep]
    df1 <- df0 %>% select(.dots = c("ID", "group", "score", all_of(col))) %>% 
      rename(ID=.dots1, group = .dots2, score=.dots3, col=.dots4) %>% 
      filter(!is.na(col))
    
    results <- boot(data=df1, statistic=rocauc, R=5, parallel="multicore")
    boot_results[[istep]] <- results
  }
  return(boot_results)
}


# ------- function to obtain R-Squared from the data (fast) -----
rocauc <- function(data, indices) {
  df1 <- data[indices,] # allows boot to select sample
  group_name <- levels(df1$group)
  nPD0 <- nrow(df1 %>% filter(group == group_name[2]))
  nHC0 <- nrow(df1 %>% filter(group == group_name[1]))
  
  # step 1
  roc_1 <- roc(df1$group, df1$score, levels=group_name,na.rm=TRUE,
               ci=TRUE, quiet=TRUE)
  cut1  <- roc_1$thresholds
  sens1 <- roc_1$sensitivities
  spec1 <- roc_1$specificities
  ncut1 <- length(cut1)
  
  # step 2
  nPD1       <- rep(NA, ncut1)
  nHC1       <- rep(NA, ncut1)
  auc_combine <- rep(NA, ncut1)
   
  for (icut1 in 1:ncut1)
  {
    # ROC of step 2
    df2 <- df1 %>% 
      filter(score >= cut1[icut1]) %>% 
      droplevels()
    nPD1[icut1] <- nrow(df2 %>% filter(group==group_name[2]))
    nHC1[icut1] <- nrow(df2 %>% filter(group==group_name[1]))
    
    if (sum(df2$group==group_name[1])>=3)
    {
      roc_2 <- roc(df2$group, df2$col, levels=group_name,na.rm=TRUE,
                   ci=TRUE, quiet=TRUE)
      
      if (icut1 == 1) # step 2 test only
      {
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
    }
  }
  ind_max <- which(auc_combine==max(auc_combine, na.rm=TRUE))
  if(length(ind_max) > 1) ind_max <- ind_max[1]
  auc_stepmax <- auc_combine[ind_max]
  return(auc_stepmax)
}