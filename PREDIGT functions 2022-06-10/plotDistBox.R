plotDistBox <- function(dfPlt, cut_BL, EVENT_vec, group_name, path, tableStr, My_Theme,
                    save_plot = TRUE, cutline = TRUE, violin = FALSE,
                    xmax = -1, adj = 1.2) {
  nvisit <- length(EVENT_vec)
  
  df <- dfPlt %>% 
    select(EVENT_ID, group, score) %>% 
    filter(!is.na(score)) %>% filter(EVENT_ID==EVENT_vec[1]) %>% 
    droplevels()
  levels(df$group) <- c("HC","PD")
  
  scrmin  <- floor(min(df$score, na.rm=TRUE))
  scrmax  <- ceiling(max(df$score, na.rm=TRUE))
  
  nHC     <- nrow(df %>% filter(group == "HC"))
  nPD     <- nrow(df %>% filter(group == "PD"))

  if (xmax == -1) xmax <- ceil(max(df$score))
  if (save_plot){
    file_name <- paste(path,"/Model 1/", tableStr," density.tiff", sep="")
    tiff(file_name, width = 1400, height = 800)
  }
  
  # density Plot
  dp <- ggplot(df, aes(x=score, fill=group)) + geom_density(alpha = 0.3, adjust = adj)
  dp <- dp + scale_fill_manual(name=NULL,labels = c(paste("HC (n=",nHC,")",sep=""), 
                                                    paste("PD (n=",nPD,")",sep="")), 
                               values=c("blue","red","green","orange"))
  dp <- dp + xlim(0,xmax)
  xmind <- layer_scales(dp)$x$range$range[1]
  xmaxd <- layer_scales(dp)$x$range$range[2]
  ymind <- layer_scales(dp)$y$range$range[1]
  ymaxd <- layer_scales(dp)$y$range$range[2]
  if (cutline)
  {
    dp <- dp + 
      geom_vline(xintercept = cut_BL, color= "red", linetype = 2, size=1) 
  }
  dp <- dp + coord_fixed(ratio=(xmaxd-xmind)/(ymaxd-ymind))
  dp <- dp + labs(x = "PREDIGT Score", y= "Density")
  
  if (save_plot)
  {
    dp <- dp + theme(legend.position = c(0.8, 0.95),
                     legend.text = element_text(size = font_size))
    dp <- dp + My_Theme
  }
  
  score1 <- df$score[df$group=="HC"]
  score2 <- df$score[df$group=="PD"]
  
  t.test(score1, score2, "two.sided", mu = 0, 
         paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  
  # box Plot
  bp <- ggplot(df,aes(x=reorder(group, score, FUN = median),y=score, fill=group)) +
    geom_boxplot(alpha = 0.3) +
    scale_fill_manual(values=c("HC"="blue",
                               "PD"="red"))
  bp <- bp + labs(x=NULL,y = "PREDIGT Score")
  bp <- bp + scale_x_discrete(labels=c("HC"=       paste("HC",  "\n",    round(mean(score1,na.rm=TRUE),2)," (",round(sd(score1,na.rm=TRUE),2),")",sep=""),
                                       "PD"=       paste("PD",  "\n",  round(mean(score2,na.rm=TRUE),2)," (",round(sd(score2,na.rm=TRUE),2),")",sep="")))
  bp <- bp + ylim(0,xmax)
  if (cutline)
  {
    bp <- bp + geom_hline(yintercept = cut_BL, color= "red", linetype = 2, size=1)
  }
  yminb <- 0 #layer_scales(bp)$y$range$range[1]
  ymaxb <- xmax #layer_scales(bp)$y$range$range[2]
  
  bp <- bp + stat_compare_means(label = "p.signif", method = "t.test",
                                label.x = 1.35, label.y = ymaxb,
                                size=label_size) # p-value
  bp <- bp + coord_fixed(ratio=2/(ymaxb-yminb))
  if (save_plot)
  {
    bp <- bp + theme(legend.position = "none")
    bp <- bp + My_Theme
  }
  
  # violin Plot
  if (violin)
  {
    vp <- ggplot(df,aes(x=reorder(group, score, FUN = median),y=score, fill=group)) +
      geom_violin(alpha = 0.3, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values=c("HC"="blue",
                                 "PD"="red"))
    vp <- vp + labs(x=NULL,y = "PREDIGT Score")
    vp <- vp + scale_x_discrete(labels=c("HC"=       paste("HC",  "\n",    round(mean(score1,na.rm=TRUE),2)," (",round(sd(score1,na.rm=TRUE),2),")",sep=""),
                                         "PD"=       paste("PD",  "\n",  round(mean(score2,na.rm=TRUE),2)," (",round(sd(score2,na.rm=TRUE),2),")",sep="")))
    if (cutline)
    {
      vp <- vp + geom_hline(yintercept = cut_BL, color= "red", linetype = 2, size=1)
    }
    yminb <- layer_scales(vp)$y$range$range[1]
    ymaxb <- layer_scales(vp)$y$range$range[2]
    vp <- vp + stat_compare_means(label = "p.signif", method = "t.test",
                                  label.x = 1.35, label.y = ymaxb,
                                  size=label_size) # p-value
    vp <- vp + coord_fixed(ratio=2/(ymaxb-yminb))
    if (save_plot)
    {
      vp <- vp + theme(legend.position = "none")
      vp <- vp + My_Theme
    }
  }
  
  # combine plots
  if (violin)
  {
    p <- ggarrange(dp, vp,ncol = 2, nrow = 1)
  } else
  {
    p <- ggarrange(dp, bp, ncol = 2, nrow = 1)
  }
  print(p)
  if (save_plot){
    dev.off()
  }
  return(list(p = p, dp = dp, bp = bp))
}