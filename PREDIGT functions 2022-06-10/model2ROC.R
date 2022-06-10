# scatter plot of the 2-step model
model2ROC <- function(Model2, step2, indTest, path, tableStr, My_Theme,
                      save_plot = TRUE, plotGray = FALSE, 
                      plotPPV = FALSE, prevalence = 0.01) {
  ntest      <- length(indTest)
  step2vec   <- step2[indTest]
  resultList <- Model2$resultList
  step1List  <- Model2$step1List
  combList   <- Model2$combList
  
  for (itest in 1:length(step2vec))
  {
    col <- step2vec[itest]
    
    # step 1
    roc_1 <- step1List[[indTest[itest]]]
    cut1  <- roc_1$thresholds
    
    # best combined
    auc_veci <- resultList[[indTest[itest]]]$auc_combine
    sensC    <- resultList[[indTest[itest]]]$sens_opt_combine
    specC    <- resultList[[indTest[itest]]]$spec_opt_combine
    combineList <- combList[[indTest[itest]]]
    indmax   <- which(auc_veci == max(auc_veci, na.rm = TRUE))
    
    if (save_plot)
    {
      setwd(paste(path,"/2-step model",sep=""))
      tiff(paste(tableStr, col, "ROC.tiff"), width = 600, height = 600)
    }
    
    rp <- ggroc(roc_1, size=1) + 
      #geom_point(data=data.frame(x=roc_1$specificities, y=roc_1$sensitivities),aes(x,y),size =2, shape=17)+ #
      geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), linetype = "dotted", color = "black")+
      coord_fixed()
    
    #gray zone
    if (plotGray)
    {
      for (icut1 in 2:length(cut1))
      {
        if (icut1 != indmax)
        {
          df_plt <- data.frame(sens=combineList[[icut1]]$sens, 
                               spec=combineList[[icut1]]$spec)
          rp <- rp + geom_path(data=df_plt, mapping=aes(spec, sens),size = 1, color=cbPalette[1],alpha=.1) +
            theme(legend.text = element_text(size = font_size))+My_Theme
        }
      }
    }
    
    # step2 only
    df_plt <- data.frame(sens=combineList[[1]]$sens, 
                         spec=combineList[[1]]$spec)
    rp <- rp + 
      geom_path(data=df_plt, mapping = aes(spec, sens),size =1, color="purple2") 
    
    # best combined
    df_plt <- data.frame(sens=combineList[[indmax]]$sens, 
                         spec=combineList[[indmax]]$spec)
    rp <- rp + 
      geom_path(data=df_plt, mapping = aes(spec, sens),size =1, color="springgreen4")+
      annotate("point", x = specC[indmax], y = sensC[indmax],colour = "darkorange", size=5, shape=15) 
    
    # equiPPV
    if (plotPPV)
    {
      sens <- c(roc_1$sensitivities[1:(indmax-1)], df_plt$sens)
      spec <- c(roc_1$specificities[1:(indmax-1)], df_plt$spec)
      
      PPV <- sens * prevalence     / (sens * prevalence + (1-spec)*(1-prevalence)) #population PPV
      NPV <- spec * (1-prevalence) / ((1 - sens) * prevalence + spec*(1-prevalence)) #population NPV
      slopeP <- (1-prevalence)/prevalence*PPVc/(1-PPVc)
      
      nPD <- resultList[[indTest[itest]]]$nPD2[1]
      nHC <- resultList[[indTest[itest]]]$nHC2[1]
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
      
      rp1 <- rp + 
        geom_path(data=dfPPV, mapping = aes(spec, PPV), size =1, color="red")+
        geom_hline(yintercept = PPVc, color="black")+
        annotate("segment", x = 1, xend = 1-1/slopeP, y = 0, yend = 1,
                 size =1, colour = "blue")

      print(rp1)
    }
    
    rp <- rp + 
      labs(title = col,
           color = "Integrated ROC",
           xlab = "1- Specificity",
           ylab = "Sensitivity")+
      theme(legend.text = element_text(size = font_size))+My_Theme
    
    print(rp)
    
    if (save_plot)
    {
      dev.off()
    }
  }
} 
