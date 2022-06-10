# boxplot of the 2-step model
model2box <- function(Model2, step2, indTest, path, tableStr, My_Theme,
                      save_plot = TRUE) {
  ntest      <- length(indTest)
  resultList <- Model2$resultList
  step1List  <- Model2$step1List
  
  auc_vec  <- resultList[[indTest[1]]]$auc_combine
  test_vec <- rep(step2[indTest[1]], length(auc_vec))
  auc_stepmax <- max(resultList[[1]]$auc_combine,na.rm=TRUE)
  auc_step2   <- resultList[[1]]$auc_combine[1]
  auc_step1   <- step1List[[indTest[1]]]$auc
  
  for (itest in 2:length(indTest))
  {
    auc_veci <- resultList[[indTest[itest]]]$auc_combine
    auc_vec  <- c(auc_vec, auc_veci)
    test_vec <- c(test_vec, rep(step2[indTest[itest]], length(auc_veci)))
    auc_stepmax <- c(auc_stepmax, max(resultList[[indTest[itest]]]$auc_combine,na.rm=TRUE))
    auc_step2   <- c(auc_step2, resultList[[indTest[itest]]]$auc_combine[1])
    auc_step1   <- c(auc_step1, step1List[[indTest[itest]]]$auc)
  }
  
  df_plt <- data.frame(value=auc_vec, test=test_vec)
  df_plt <- df_plt %>% 
    mutate(test = factor(test, levels = step2))
  df_plt2 <- data.frame(auc=c(auc_stepmax, auc_step2, auc_step1), 
                        type=c(rep(1,ntest),rep(2,ntest),rep(3,ntest)),
                        test=rep(step2,3))
  df_plt2 <- df_plt2 %>% 
    mutate(test = factor(test, levels = step2))
  
  if (save_plot)
  {
    setwd(paste(path,"/2-step summary box",sep=""))
    tiff(paste(tableStr, "AUC summary box.tiff"), width = 573, height = 345)
  }
  
  bp <- ggplot(df_plt, aes(x=test, y = value)) +
    geom_boxplot(color="gray65") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
    geom_point(data=df_plt2, aes(factor(test) , auc, color = factor(type)),
               size=3)+
    scale_color_manual(name="",
                       labels = c("Largest combined AUC", "AUC of step2 test only", "AUC of step1"),
                       values = c("springgreen4","purple2","black"))+ # values = c(15,8,17),
    labs(y = "combined AUC")+
    ylim(0.5, 1)+
    theme(axis.title.x = element_blank(),
          legend.text = element_text(size = 16))+My_Theme
  print(bp)
  
  if (save_plot)
  {
    dev.off()
  }
} 
