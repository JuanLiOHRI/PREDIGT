PREDIGTScr <- function(df_LR_Scr, varStr1, varStr2, coefV, ind_question, EVENT_vec, FG, FT, method = 1) {
  nvisit <- length(EVENT_vec)
  
  df_LR_Scr <- df_LR_Scr %>% 
    mutate(min_E = case_when(
      age <= 50 ~ 0,
      age < 60 ~ 0.005,
      age < 70 ~ 0.0075,
      age < 80 ~ 0.02,
      TRUE ~ 0.03
    )) %>% 
    mutate(min_I = case_when(
      age <= 50 ~ 0,
      age < 60 ~ 0.001,
      age < 70 ~ 0.002,
      age < 80 ~ 0.003,
      TRUE ~ 0.004
    ))
  
  if (method == 1){
    # 1. Score with the original coefficients, w age-associated coefficients
    for (ivisit in 1:nvisit)
    {
      dfScr <- df_LR_Scr %>% filter(EVENT_ID == EVENT_vec[ivisit]) %>% 
        select(all_of(varStr1), min_E, min_I)
      
      # check data
      # dfScr[,c(2,3,5:ncol(dfScr))] <- lapply(dfScr[,c(2,3,5:ncol(dfScr))],factor)
      # summary(dfScr)
      
      dfScr2 <- dfScr %>% select(all_of(varStr2))
      if (!is.null(ind_question))
      {
        dfScr2[names(dfScr2)[ind_question]] <- lapply(dfScr2[names(dfScr2)[ind_question]], question)
      }
      # dfScr2[,c(1:ncol(dfScr2))] <- lapply(dfScr2[,c(1:ncol(dfScr2))],factor)
      # summary(dfScr2)
      
      dfScr2 <- as.matrix(dfScr2)
      score_E <- t(coefV[ind_E] %*% t(dfScr2)[ind_E,])
      score_E <- ifelse(score_E<=0, dfScr$min_E,score_E)
      
      score_D <- t(coefV[ind_D] %*% t(dfScr2)[ind_D,])
      score_D <- ifelse(score_D<=0, 0.01,score_D)
      
      score_I <- t(coefV[ind_I] %*% t(dfScr2)[ind_I,])
      score_I <- ifelse(score_I<=0, dfScr$min_I,score_I)
      
      score   <- score_E + score_D + score_I
      
      if (FG) # factor G
      {
        score <- ifelse(dfScr$sex=="Female",score*0.8,score*1.2)
      }
      if (FT) # factor T
      {
        score <- score * dfScr$age
      }
      
      df_LR_Scr[which(df_LR_Scr$EVENT_ID == EVENT_vec[ivisit]),"score"] <- score
    }
  } else{
    # 2. score using simple additive LR *G*T
    for (ivisit in 1:nvisit)
    {
      dfScr <- df_LR_Scr %>% filter(EVENT_ID == EVENT_vec[ivisit]) %>% 
        select(all_of(varStr1))
      dfScr2 <- dfScr %>% select(all_of(varStr2)) %>% 
        mutate(group = dfScr$group)
      if (ivisit == 1)
      {
        vStr <- names(dfScr2)[names(dfScr2)!="group"]
        vStr <- paste(vStr, collapse="+")
        
        formula  <- as.formula(paste("group ~", vStr, sep=""))
        lr_score <- rms::lrm(formula, dfScr2, x=TRUE, y=TRUE)
        lr_score
        bsv <- rms::validate(lr_score, B=300) # bootstrapping
        cbind((bsv[1,1]+1)/2,(bsv[1,5]+1)/2) # Original AUC, Corrected AUC
        
        # lr prediction
        score <- predict(lr_score)
      } else 
      {
        score <- predict(lr_score, newdata=dfScr2)
      } 
      
      if (FG) # factor G
      {
        score <- ifelse(dfScr$sex=="Female",score*0.8,score*1.2)
      }
      if (FT) # factor T
      {
        score <- score * dfScr$age
      }
      
      df_LR_Scr[which(df_LR_Scr$EVENT_ID == EVENT_vec[ivisit]),"score"] <- score
    }
  }
  return(df_LR_Scr)
}
