################################################################################

polyg.area <- function(pol) {
  
  pol_le <- nrow(pol)
  pol <- round (pol,6)
  if (!all(pol[1,]==pol[pol_le,])) return(0.0)
  idx <- 2:pol_le
  return(abs(sum (pol[idx-1,1] * pol[idx,2] - pol[idx,1] * pol[idx-1,2]) / 2))
}






