# scatter plot of the 2-step model
model2point <- function(Model2, df_LR, step2, indTest, path, tableStr, My_Theme,
                      save_plot = TRUE) {
  ntest      <- length(indTest)
  step2vec   <- step2[indTest]
  resultList <- Model2$resultList
  step1List  <- Model2$step1List
  
  df0 <- df_LR %>% 
    filter(EVENT_ID == "baseline") %>% 
    select(ID, group, score,names(df_LR[which(names(df_LR) %in% step2vec)])) %>% 
    filter(!is.na(score)) %>% 
    droplevels()
  
  if (save_plot)
  {
    setwd(paste(path,"/2-step correlation",sep=""))
  }
  
  for (itest in 1:length(indTest))
  {
    col <- step2vec[itest]
    df1 <- df0 %>% select(.dots = c("ID", "group", "score", all_of(col))) %>% 
      rename(subID=.dots1, group=.dots2,  score = .dots3, col =.dots4) %>% 
      filter(!is.na(col))
    
    auc_veci <- resultList[[indTest[itest]]]$auc_combine
    cut1     <- step1List[[indTest[itest]]]$thresholds
    cut2     <- resultList[[indTest[itest]]]$cut_opt_combine
    indmax   <- which(auc_veci == max(auc_veci, na.rm = TRUE))
    
    cut_s1 <- cut1[indmax]
    cut_s2 <- cut2[indmax]
    
    if (save_plot)
    {
      tiff(paste(tableStr, col, "correlation.tiff"), width = 600, height = 600)
    }
    
    # classic plot :
    p <- ggplot(df1, aes(score, col, color=group)) +
      geom_point() +
      theme(legend.position="bottom",
            legend.text = element_text(size = font_size))+
      labs(y=col,
           x="Score of step 1") + 
      scale_color_manual(name=NULL,labels = c("HC","PD"), 
                         values=c("blue","red")) +
      geom_vline(xintercept = cut_s1, color= "springgreen4", linetype = 2, size=1)+
      geom_hline(yintercept = cut_s2, color= "darkorange", linetype = 2, size=1)+
      My_Theme
    
    # Set relative size of marginal plots (main plot 10x bigger than marginals)
    p <- ggExtra::ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
    
    print(p)
    
    if (save_plot)
    {
      dev.off()
    }
  }
} 
