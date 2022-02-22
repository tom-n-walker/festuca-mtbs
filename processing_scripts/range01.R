################################################################################
#### Project: Metabolism traits
#### Title:   Function | Small function | Scale between 0 and 1
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    19 August 2020
#### ---------------------------------------------------------------------------

range01 <- function(x){
  maxX <- max(x, na.rm = T)
  minX <- min(x, na.rm = T)
  scale <- (x - minX) / (maxX - minX)
  if(sum(is.na(scale)) == length(scale)){
    scale <- rep(0, length(scale))
  }
  return(scale)
}
