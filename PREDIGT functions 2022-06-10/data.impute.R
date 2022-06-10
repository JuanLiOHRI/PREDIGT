# imputate missing data using previous or next non-NA value
data.impute <- function(data_join, plt = FALSE) {
  
  if (plt)
  {
    naniar::vis_miss(data_join %>% arrange(EVENT_ID))
  }
  
  data_imp <- data_join %>% select_if(is.numeric)
  data_imp <- cbind(data_imp, data_join %>% select_if(is.character))
  data_imp <- data_imp %>% 
    purrrlyr::slice_rows("ID") %>% 
    purrrlyr::by_slice(function(x) { 
      zoo::na.locf(zoo::na.locf(x, na.rm=FALSE), fromLast=TRUE, na.rm=FALSE) }, 
      .collate = "rows") 
  data_imp <- cbind(data_imp, data_join %>% select_if(is.factor))
  #data_imp <- data_imp[,names(data_join)[! names(data_join) %in% c("time0","dateR2")]] # reorder by column name
  
  if (plt)
  {
    naniar::vis_miss(data_imp %>% arrange(EVENT_ID))
  }
  
  return(data_imp)
} 
