# integrated score of each variable 
varScore <- function(data, varList, varStr2, EVENT_vec, group_name) {
  nvisit    <- length(EVENT_vec)
  nV     <- length(varStr2)
  cutVec <- rep(0,nV)
  out <- NULL
  
  # ----- baseline ----- 
  dflr <- data %>% filter(EVENT_ID == EVENT_vec[1])
  
  for (ivar in 1:nV)
  {
    vStr <- varList[[ivar]]
    vStr <- paste(vStr, collapse="+")
    
    formula  <- as.formula(paste("group ~", vStr, sep=""))
    lr  <- rms::lrm(formula, dflr, x=TRUE, y=TRUE)
    lr
    bsv <- rms::validate(lr, B=300) # bootstrapping
    cbind((bsv[1,1]+1)/2,(bsv[1,5]+1)/2) # Original AUC, Corrected AUC
    
    # lr prediction
    var <- paste(varStr2[ivar],".LR", sep="")
    dflr <- dflr %>% mutate(!!var:=predict(lr))
    
    datai <- dflr[,ncol(dflr)]
    roc_1 <- pROC::roc(dflr$group, datai, levels=group_name,na.rm=TRUE)
    # find finite thresholds
    ind1  <- which(!is.infinite(roc_1$thresholds))
    sum   <- roc_1$sensitivities[ind1]+roc_1$specificities[ind1]
    cuti  <- roc_1$thresholds[ind1]
    sensi <- roc_1$sensitivities[ind1]
    speci <- roc_1$specificities[ind1]
    # find 'optimal'
    ind2  <- which(sum==max(sum))
    cuti  <- cuti[ind2]
    sensi <- sensi[ind2]
    speci <- speci[ind2]
    if (length(cuti)>1)
    {
      cuti  <- cuti[1]
      sensi <- sensi[1]
      speci <- speci[1]
    }
    
    # determine if variable is true
    var <- paste(varStr2[ivar],".PREDIGT", sep="")
    dflr <- dflr %>% 
      mutate(!!var:=ifelse(is.na(datai), NA, ifelse(datai>cuti, 1, 0)))
    
    out[[ivar]] <- lr
    cutVec[ivar] <- cuti
  }
  df_LR <- dflr
  
  # ----- 24 and 48 using lr models at baseline-----
  for (ivisit in 2:nvisit)
  {
    dflr <- data %>% filter(EVENT_ID == EVENT_vec[ivisit])
    
    for (ivar in 1:nV)
    {
      # lr prediction
      var <- paste(varStr2[ivar],".LR", sep="")
      dflr <- dflr %>% mutate(!!var:=predict(out[[ivar]], newdata = dflr))
      
      # determin if variable is true
      datai <- dflr[,ncol(dflr)]
      var   <- paste(varStr2[ivar],".PREDIGT", sep="")
      dflr  <- dflr %>% 
        mutate(!!var:=ifelse(is.na(datai), NA, ifelse(datai>cutVec[ivar], 1, 0)))
    }
    
    # save complete data
    df_LR <- rbind(df_LR, dflr)
  }
  
  return(df_LR)
} 
