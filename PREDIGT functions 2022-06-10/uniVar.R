# univariate results
uniVar <- function(data_LR, varStr, EVENT_vec, group_name, pltStr = NULL, pltStr2 = NULL) {
  data_LR <- data_LR %>% mutate(sex = as.numeric(sex))
  nVar <- length(varStr)
  cutVec <- c(0,nVar)
  nvisit <- length(EVENT_vec)
  
  results <- list()
  resPlt  <- data.frame(variable = "",
                        nHC = 0,
                        pHC = 0,
                        nPD = 0,
                        pPD = 0,
                        nTotal = 0,
                        pTotal = 0,
                        OR = 0,
                        ORl = 0,
                        ORu = 0)
  for (ivisit in 1:nvisit)
  {
    df <- data_LR %>% filter(EVENT_ID == EVENT_vec[ivisit]) %>% droplevels()
    result <- rep(NA, nVar+1)
    result[1] <- "variable; nHC (%); nPD (%); OR (95% CI); AUC (95% CI); sensitivity; specificity; cut-value"
    
    for (ivar in 1:nVar)
    {
      col <- varStr[ivar]
      df1 <- df %>% select(.dots = c("ID", "group", all_of(col))) %>% 
        rename(ID=.dots1, group=.dots2, col=.dots3) %>% 
        filter(!is.na(col)) %>% 
        droplevels()
      
      if (nrow(df1) >6 & length(levels(df1$group)) == 2)
      {
        datai <- df1$col
        df1 <- df1 %>% mutate(col.fct = as.factor(col))
        roc_1 <- pROC::roc(df1$group, datai, levels=group_name,na.rm=TRUE,ci=TRUE)
        # find finite thresholds
        ind1  <- which(!is.infinite(roc_1$thresholds))
        sum   <- roc_1$sensitivities[ind1]+roc_1$specificities[ind1]
        cuti  <- roc_1$thresholds[ind1]
        sensi <- roc_1$sensitivities[ind1]
        speci <- roc_1$specificities[ind1]
        if (ivisit == 1)
        {
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
          cutVec[ivar] <- cuti
        } else
        {
          if (length(cuti)>1)
          {
            ind   <- which.min(abs(cuti - cutVec[ivar]))
            cuti  <- cuti[ind]
            sensi <- sensi[ind]
            speci <- speci[ind]
          }
        }
        
        aucStr1 <- paste(round(roc_1$auc,2), " (",round(roc_1$ci[1],2),"-",round(roc_1$ci[3],2),")",
                         "; ", round(speci,2), 
                         "; ", round(sensi,2), 
                         "; ", roc_1$direction,round(cuti, 2),sep="")
        
        # OR
        scorei    <- rep(0,length(datai))
        scorei[is.na(datai)] <- NA
        if (roc_1$direction == "<")
        {
          scorei[datai >= cuti] <- 1
        } else
        {
          scorei[datai <= cuti] <- 1
        }
        
        M    <- table(df1$group, scorei)
        OR   <- fisher.test(M)$estimate
        ORl  <- fisher.test(M)$conf.int[1]
        ORu  <- fisher.test(M)$conf.int[2]
        orStr1 <- paste(round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
        
        df1 <- df1 %>% mutate(score = factor(scorei))
        total1 <- df1 %>% group_by(group) %>% summarise(count=n()) %>% pull(count)
        total2 <- df1 %>% group_by(group,score) %>% summarise(count=n()) %>% pull(count)
        
        nHC    <- total2[2]
        pHC    <- total2[2]/total1[1]*100
        nPD    <- total2[4]
        pPD    <- total2[4]/total1[2]*100
        nTotal <- nHC + nPD
        pTotal <- round(nTotal/(total1[1] + total1[2])*100,2)
        
        result[ivar+1] <- paste(col, ";", 
                                nHC, "(",round(pHC,2),"); ",
                                nPD, "(",round(pPD,2),"); ",
                                orStr1, "; ", 
                                aucStr1,sep="")
      }
      results[[ivisit]] <- result
      
      if (ivisit == 1 & col %in% pltStr)
      {
        resPlti <- data.frame(variable = col,
                              nHC = nHC,
                              pHC = pHC,
                              nPD = nPD,
                              pPD = pPD,
                              nTotal = nTotal,
                              pTotal = pTotal,
                              OR = OR,
                              ORl = ORl,
                              ORu = ORu)
        resPlt <- bind_rows(resPlt, resPlti)
      }
    }
  }
  
  if (!is.null(pltStr))
  {
    resPlt <- resPlt[-1,]
    p1 <- ggplot(resPlt, aes(pPD, OR, color = factor(variable), label = factor(pltStr2)))+
      geom_point()+
      geom_text(position = position_nudge(y = -0.3))+
      geom_hline(yintercept = 1, color = "red")+
      xlim(0,100)+
      scale_color_discrete(label = factor(pltStr2))+
      xlab("Freq. in PD cases (%)")+
      ylab("Odds Ratio")
    p2 <- ggplot(resPlt, aes(pTotal, OR, color = factor(variable), label = factor(pltStr2)))+
      geom_point()+
      geom_text(position = position_nudge(y = -0.3))+
      geom_hline(yintercept = 1, color = "red")+
      xlim(0,100)+
      scale_color_discrete(label = factor(pltStr2))+
      xlab("Freq. in the whole cohort (%)")+
      ylab("Odds Ratio")
    
    print(ggarrange(p1, p2, ncol = 2, common.legend = TRUE))
  }
  
  return(results)
} 
