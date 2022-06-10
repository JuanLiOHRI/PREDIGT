plotBar <- function(dataBar, col, varlab, col2 = "group", 
                    level_rvs = FALSE, axis_flip = FALSE,
                    sel_contains = FALSE, levelStr = NULL,
                    or_cal = FALSE, or_combine = FALSE,
                    group_name = c("Healthy Control", "Parkinson's Disease"), ylim = c(-1.0, 1.0),
                    legend.position = "right") {
  dataBar <- dataBar %>% droplevels()
  # summary statistics
  if (sel_contains)
  {
    df1 <- dataBar %>% 
      select(all_of(col2), contains(col)) 
    df2 <- select_if(df1, is.numeric)
    df1 <- cbind(df1[,1], df2)
    names(df1)[1] <- "group"
    
    sumHC  <- colSums(df1 %>% filter(group == group_name[1]) %>% select(-group), na.rm=TRUE)
    sumHC2 <- colSums(!is.na(df1 %>% filter(group == group_name[1]) %>% select(-group)), na.rm=TRUE)
    sumPD <- colSums(df1 %>% filter(group == group_name[2]) %>% select(-group), na.rm=TRUE)
    sumPD2 <- colSums(!is.na(df1 %>% filter(group == group_name[2]) %>% select(-group)), na.rm=TRUE)
    
    d    <- data.frame(group = c(rep(group_name[1],length(sumHC)), rep(group_name[2],length(sumHC))),
                       count1 = c(sumHC2, sumPD2),
                       n = c(sumHC, sumPD),
                       col = rep(levelStr,2))
    d <- d %>% mutate(freq = ifelse(group == group_name[2],n/count1,-n/count1))
  } else
  {
    df1 <- dataBar %>% select(.dots = c("ID", all_of(col2), all_of(col))) %>% 
      rename(ID=.dots1, group=.dots2, col=.dots3) %>% 
      mutate(group = factor(group)) %>% 
      mutate(col = factor(col))
    if (!is.null(levelStr))
      levels(df1$col) <- levelStr
    lev <- levels(df1$group)
    d <- df1 %>% 
      filter(!is.na(group) & !is.na(col)) %>% 
      group_by(group,col) %>% 
      summarise(n = n())%>%
      mutate(freq = n / sum(n))
      
    if (axis_flip) {
      d <- d %>% mutate(freq = ifelse(col==levels(col)[1], freq, -freq)) 
    } else {
      d <- d %>% mutate(freq = ifelse(group==lev[2], freq, -freq)) 
    }
  }
  
  if (level_rvs)
    d$col <- factor(d$col, levels=rev(levels(d$col)))
  
  # plot
  if (axis_flip) {
    p <- ggplot(d, aes(x=group,y=freq,fill=col))+
      geom_bar(stat='identity', alpha=0.5)+
      geom_text(aes(x = group, y = freq, label = n), 
                position = position_dodge(width = 1)) +
      coord_flip()+
      scale_fill_discrete(name = varlab)+
      scale_y_continuous(name="Proportion (%)", breaks = c(-1.0, -0.5, 0.0, 0.5, 1.0), 
                         labels= c(1.0, 0.5, 0.0, 0.5, 1.0), limits=ylim)+
      xlab(varlab) +
      theme(legend.position=legend.position)
  } else {
    p <- ggplot(d, aes(x=col,y=freq,fill=group))+
      geom_bar(stat='identity', alpha=0.5)+
      geom_text(aes(x = col, y = freq, label = n), 
                position = position_dodge(width = 1)) +
      coord_flip()+
      scale_fill_discrete(name = col2)+
      scale_y_continuous(name="Proportion (%)", breaks = c(-1.0, -0.5, 0.0, 0.5, 1.0), 
                         labels= c(1.0, 0.5, 0.0, 0.5, 1.0), limits=ylim)+
      xlab(varlab) +
      theme(legend.position=legend.position)
  }
  
  if (col2 == "group")
  {
    p <- p+
      scale_fill_manual(values=c("Healthy Control"="blue",
                                 "Parkinson's Disease"="red"))
    #print(p)
  }
  
  # Calculate Odds ratio (optional)
  if (or_cal)
  {
    if (sel_contains)
    {
      for (i in 2:ncol(df1))
      {
        coli <- names(df1)[i]
        value <- df1 %>% select(all_of(coli))
        dfi  <- data.frame(group = df1$group, value = value)
        names(dfi)[2] <- "value"
        M    <- table(dfi$group, dfi$value)
        if (ncol(M) == 2)
        {
          OR   <- fisher.test(M)$estimate
          ORl  <- fisher.test(M)$conf.int[1]
          ORu  <- fisher.test(M)$conf.int[2]
          orStr1 <- paste(levelStr[i-1], ": ",
                          round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
          print(orStr1)
        } else if (ncol(M) > 2)
        {
          for(j in 2:ncol(M))
          {
            Mj   <- M[,c(1,j)]
            OR   <- fisher.test(Mj)$estimate
            ORl  <- fisher.test(Mj)$conf.int[1]
            ORu  <- fisher.test(Mj)$conf.int[2]
            orStr1 <- paste(levelStr[i-1], ": ",
                            round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
            print(orStr1)
          }
        }
      }
    } else
    {
      # OR
      M    <- table(df1$group, df1$col)
      if (ncol(M) == 2)
      {
        OR   <- fisher.test(M)$estimate
        ORl  <- fisher.test(M)$conf.int[1]
        ORu  <- fisher.test(M)$conf.int[2]
        orStr1 <- paste(varlab, " (", levels(df1$col)[2],"): ",
                        round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
        print(orStr1)
      } else if (ncol(M) > 2)
      {
        for(j in 2:ncol(M))
        {
          Mj   <- M[,c(1,j)]
          OR   <- fisher.test(Mj)$estimate
          ORl  <- fisher.test(Mj)$conf.int[1]
          ORu  <- fisher.test(Mj)$conf.int[2]
          orStr1 <- paste(varlab, " (", levels(df1$col)[j],"): ",
                          round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
          print(orStr1)
        }
        if(or_combine)
        {
          Mc <- M[,1:2]
          Mc[,2] <- rowSums(M[,2:ncol(M)])
          OR   <- fisher.test(Mc)$estimate
          ORl  <- fisher.test(Mc)$conf.int[1]
          ORu  <- fisher.test(Mc)$conf.int[2]
          orStr1 <- paste(varlab, ": ",
                          round(OR,2), " (", round(ORl,2), "-", round(ORu,2), ")", sep="")
          print(orStr1)
        }
      }
    }
  }
  return(p)
}