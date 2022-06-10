plotROC <- function(dfPlt, EVENT_vec, group_name, path, tableStr, My_Theme,
                    save_plot = TRUE, cutBL = FALSE, BL_only = TRUE, plotPPV = FALSE,
                    prevalence = 1/100, PPVc = 0.1) {
  nvisit <- length(EVENT_vec)
  
  df <- dfPlt %>% 
    select(EVENT_ID, group, score) %>% 
    filter(!is.na(score))
  # remove EVENT_ID group with too few subjects.
  A <- df %>% 
    group_by(group, EVENT_ID) %>% 
    summarise(number=n()) %>% 
    filter(number<3) %>% 
    mutate(groupEVENT_ID = paste(group,EVENT_ID))
  df <- df %>% 
    mutate(groupEVENT_ID = paste(group,EVENT_ID)) %>% 
    filter(!(groupEVENT_ID %in% A$groupEVENT_ID)) %>% 
    droplevels()
  totals <- df %>% 
    group_by(group, EVENT_ID) %>% 
    summarise(number=n())
  
  result <- rep(NA,nvisit+1)
  result[1] <- "visit; nHC; nPD; AUC (95% CI); sensitivity (95% CI); specificity (95% CI); cut-value"
  rocList   <- list()
  cut_BL <- 0
  for (ivisit in 1:nvisit)
  {
    df1 <- df %>% filter(EVENT_ID == EVENT_vec[ivisit]) %>% droplevels()
    
    if (nrow(df1) >6 & length(levels(df1$group)) == 2)
    {
      datai <- df1$score
      roc_1 <- roc(df1$group, datai, levels=group_name, na.rm=TRUE,
                   ci=TRUE)
      # find finite thresholds
      ind1  <- which(!is.infinite(roc_1$thresholds))
      sum   <- roc_1$sensitivities[ind1]+roc_1$specificities[ind1]
      cuti  <- roc_1$thresholds[ind1]
      sensi <- roc_1$sensitivities[ind1]
      speci <- roc_1$specificities[ind1]
      if (ivisit == 1 | !cutBL)
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
        if (ivisit == 1)
        {
          cut_BL <- cuti
        }
        cut_v1 <- cuti
      } else
      {
        if (length(cuti)>1)
        {
          ind   <- which.min(abs(cuti - cut_v1))
          cuti  <- cuti[ind]
          sensi <- sensi[ind]
          speci <- speci[ind]
        }
      }
      
      # ci of sensitivity and specificity
      ci_cut <- ci(roc_1, of = "thresholds", thresholds = cuti)
      se_ci1 <- ci_cut$sensitivity[1]
      se_ci2 <- ci_cut$sensitivity[3]
      sp_ci1 <- ci_cut$specificity[1]
      sp_ci2 <- ci_cut$specificity[3]
      
      aucStr1 <- paste(round(roc_1$auc,2), "(",round(roc_1$ci[1],2),"-",round(roc_1$ci[3],2),")",
                       "; ", round(sensi,2), "(",round(se_ci1,2),"-",round(se_ci2,2),")",
                       "; ", round(speci,2), "(",round(sp_ci1,2),"-",round(sp_ci2,2),")",
                       "; ", round(cuti, 2),sep="")
      total <- df1 %>% group_by(group) %>% summarise(count=n()) %>% pull(count)
      result[ivisit+1] <- paste(EVENT_vec[ivisit], "; ",total[1], "; ",total[2], "; ",aucStr1,sep="")
    }
    rocList[[ivisit]] <- roc_1
  }
  
  # ROC plot
  if (BL_only)
  {
    if (save_plot){
      tiff(paste(path,"/Model 1/", tableStr," BL ROC.tiff",sep=""), width = 800, height = 1000)
    }
    
    roc_obj <- rocList[[1]]
    ciobj  <- ci(roc_obj, of = "thresholds", thresholds = roc_obj$thresholds)
    ci_cut <- ci(roc_obj, of = "thresholds", thresholds = cut_BL)
    se_ci1 <- ci_cut$sensitivity[1]
    sensi  <- ci_cut$sensitivity[2]
    se_ci2 <- ci_cut$sensitivity[3]
    sp_ci1 <- ci_cut$specificity[1]
    speci  <- ci_cut$specificity[2]
    sp_ci2 <- ci_cut$specificity[3]
    
    dat.ci <- data.frame(x = ciobj$specificity[,2],
                         lower = ciobj$sensitivity[,1],
                         upper = ciobj$sensitivity[,3])
    
    p <- ggroc(roc_obj, size=2) + 
      geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2)+
      geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), linetype = "dotted", color = "black") +
      annotate("point", x = speci, y = sensi,
               colour = "red", size=5) +
      annotate("segment", x = speci, xend = speci, y = se_ci1, yend = se_ci2,
               colour = "red", size=1) +
      annotate("segment", x = sp_ci2, xend = sp_ci1, y = sensi, yend = sensi,
               colour = "red", size=1) 
    if (plotPPV)
    {
      sens <- roc_obj$sensitivities
      spec <- roc_obj$specificities
      
      PPV <- sens * prevalence     / (sens * prevalence + (1-spec)*(1-prevalence)) #population PPV
      NPV <- spec * (1-prevalence) / ((1 - sens) * prevalence + spec*(1-prevalence)) #population NPV
      slopeP <- (1-prevalence)/prevalence*PPVc/(1-PPVc)
      
      nPD <- sum(df1$group == group_name[2])
      nHC <- sum(df1$group == group_name[1])
      TP  <- sens*nPD
      FN  <- nPD - TP
      TN  <- spec*nHC
      FP  <- nHC - TN
      PPVs <- TP/(TP+FP)
      NPVs <- TN/(TN+FN)
      
      dfPPV <- data.frame(sens=sens,
                          spec=spec,
                          PPV = PPV,
                          NPV = NPV,
                          PPVs = PPVs,
                          NPVs = NPVs)
      
      min(PPV, na.rm = TRUE)
      min(NPV, na.rm = TRUE)
      
      pp <- ggroc(roc_obj, size=2)+
        geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), linetype = "dotted", color = "black")
      
      pp <- pp + 
        geom_path(data=dfPPV, mapping = aes(spec, PPV), size =1, color="red")+
        geom_hline(yintercept = PPVc, color="black")+
        annotate("segment", x = 1, xend = 1-1/slopeP, y = 0, yend = 1,size =1, colour = "forestgreen")+
        #geom_path(data=dfPPV, mapping = aes(spec, NPV), size =1, color="blue") + 
        #geom_path(data=dfPPV, mapping = aes(spec, PPVs), size =1, color="red", linetype = "dashed") + 
        #geom_path(data=dfPPV, mapping = aes(spec, NPVs), size =1, color="blue", linetype = "dashed") + 
        coord_fixed()
    }
    
    p <- p + coord_fixed() + My_Theme
    p <- annotate_figure(p,
                         top = text_grob(paste(result[2],collapse="\n"), face = "bold", size = 14))
    print(p)
    if (save_plot){
      dev.off()
    }
  } else
  {
    if (save_plot){
      setwd(paste(path,"/ROC",sep=""))
      tiff(paste(tableStr,"ROC.tiff"), width = 800, height = 1000)
    }
    p <- ggroc(rocList, size=2) + 
      geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), linetype = "dotted", color = "black")
    p <- p + coord_fixed() + My_Theme
    p <- annotate_figure(p,
                         top = text_grob(paste(result,collapse="\n"), face = "bold", size = 14))
    print(p)
    if (save_plot){
      dev.off()
    }
  }
  
  return(list(cut_BL = cut_BL, result = result))
}