
# Distance to conformity
dist.to.conformity <- function (ref_inter,test_inter, inter.vol){
  under = NULL
  over = NULL
  
  ijk.over <- get.ijk.from.index(idx=which(test_inter$vol3D.data), vol=test_inter)
  
  ijk.under <- get.ijk.from.index(which(ref_inter$vol3D.data), ref_inter)
  
  if(!is.null(ijk.over)) {
    if (is.matrix(ijk.over)) {
      over <- .mdcC(as.numeric(inter.vol$vol3D.data),as.numeric(inter.vol$n.ijk),as.numeric(inter.vol$dxyz),
                    ijk.over[,1],ijk.over[,2],ijk.over[,3])
    } else {
      over <- .mdcC(as.numeric(inter.vol$vol3D.data),as.numeric(inter.vol$n.ijk),as.numeric(inter.vol$dxyz),
                    ijk.over[1],ijk.over[2],ijk.over[3])
    }
  }
  
  if(!is.null(ijk.under)) {
    if (is.matrix(ijk.under)) {
      under <- .mdcC(as.numeric(inter.vol$vol3D.data),as.numeric(inter.vol$n.ijk),as.numeric(inter.vol$dxyz),
                     ijk.under[,1],ijk.under[,2],ijk.under[,3])
    } else {
      under <- .mdcC(as.numeric(inter.vol$vol3D.data),as.numeric(inter.vol$n.ijk),as.numeric(inter.vol$dxyz),
                     ijk.under[1],ijk.under[2],ijk.under[3])
    }
  }
  
  # under.MDC <- ifelse(is.null(under),0, mean(under, na.rm = TRUE))
  # over.MDC <- ifelse(is.null(over),0, mean(over, na.rm = TRUE))
  
  return(list(under.contouring = under,over.contouring = over))
}

