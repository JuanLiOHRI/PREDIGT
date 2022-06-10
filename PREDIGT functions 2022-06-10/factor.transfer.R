# transfer categorical variables into numeric
as.numeric.factor <- function(x) {
  if (length(levels(x))>0)
  {
    as.numeric(levels(x))[x]
  } else
  {
    x
  }
} 

# transfer categorical variables into numeric
as.numeric.factor2 <- function(x) {
  if (is.factor(x)) as.numeric(x)-1
  else x
} 

# transfer categorical variables into 2-levels
question <- function(x) {
  for (i in 1:length(x))
  {
    if (!is.na(x[i]))
    {
      if (x[i] > 0)
      {
        x[i] <- 1
      }
    }
  }
  return(x)
}