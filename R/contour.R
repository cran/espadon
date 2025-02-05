.contour.level <- function(contour.L){
  rg<- (0:length(contour.L))[-1]
  level <- rep(0,length(rg))
  if (length(rg)>1) {
    level <- sapply (rg, function(j){
      ptj<- contour.L[[j]]
      roi.index.k <- rg[-j]
      # if (length(roi.index.z)!=0) {
      r <- unique (sapply (roi.index.k, function (k) {
        ptk <- contour.L[[k]]
        keep <- .pt.in.polygon (ptj[1,1], ptj[1,2], ptk[ ,1], ptk[ ,2], 1e-20) > 0.5
        return (ifelse (all(keep), k,NA))}))
      r <- r[!is.na (r)]
      return(length(r))
    }) 
  }
  return(level)
}
