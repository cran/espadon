.polyg.is.clockwise <- function(pol) {
  le <- nrow(pol)
  if (le==0) return(FALSE)
  
  first <- order(pol[1:(le-1),2],pol[1:(le-1),1])[1]
  pol <- rbind(pol[first:(le-1),], pol[-(first:le),])
  r <- vector.product(pol[1,]-pol[le-1,],pol[2,]-pol[1,])
  while (r[3]==0) {
    pol <- pol[-1,]
    le <- length(pol)
    if (length(pol)==0) return(FALSE)
    first <- order(pol[,2],pol[,1])[1]
    pol <- rbind(pol[first:le,], pol[-(first:le),])
    r <- vector.product(pol[1,]-pol[le-1,],pol[2,]-pol[1,])
  }
  return(r[3]<0)
}