.vol.border.tuning <- function (vol, pre.nijk = c (0, 0, 0), post.nijk = pre.nijk) {
  if (all(pre.nijk == c(0,0,0)) &  all(post.nijk == c(0,0,0))) return(vol)
  
  new.vol <- vol
  
  new.vol$n.ijk <- vol$n.ijk + pre.nijk + post.nijk
  pt000 <- c(-pre.nijk[1:2],vol$k.idx[1]-pre.nijk[3], 1) %*% t(vol$xyz.from.ijk)
  new.vol$xyz.from.ijk[ ,4] <- pt000
  
  new.vol$k.idx <- (vol$k.idx[1]-pre.nijk[3]):(vol$k.idx[length(vol$k.idx)]+post.nijk[3])
  in.k <- new.vol$k.idx >=vol$k.idx[1] &  new.vol$k.idx <= vol$k.idx[length(vol$k.idx)]
  to.suppress <- which(is.na(match (new.vol$k.idx[in.k], vol$k.idx)))
  if (length(to.suppress)>0){
    new.vol$k.idx <- sort(c(new.vol$k.idx [!in.k], new.vol$k.idx [in.k][-to.suppress]))
    new.vol$missing.k.idx <- TRUE
  }
  new.vol$k.idx <-  new.vol$k.idx + pre.nijk[3]
  
  new.vol$xyz0  <- matrix ((as.matrix (expand.grid (0, 0, new.vol$k.idx,1))%*% t(new.vol$xyz.from.ijk))[ ,1:3],ncol=3)
  new.vol$cube.idx <- matrix ( c(0, 0, 0, 1, new.vol$n.ijk[1]-1, 0, 0, 1,
                                 new.vol$n.ijk[1]-1, new.vol$n.ijk[2]-1, 0, 1, 0, new.vol$n.ijk[2]-1, 0, 1,
                                 0, 0, new.vol$n.ijk[3]-1, 1, new.vol$n.ijk[1]-1, 0, new.vol$n.ijk[3]-1, 1,
                                 new.vol$n.ijk[1]-1, new.vol$n.ijk[2]-1, new.vol$n.ijk[3]-1, 1, 0, new.vol$n.ijk[2]-1, new.vol$n.ijk[3]-1, 1),
                               nrow=4, byrow= FALSE)
  
  new.vol$vol3D.data <- array(NA, dim=new.vol$n.ijk)
  
  flag.i <- max(1, 1 - pre.nijk[1]) : min (vol$n.ijk[1], vol$n.ijk[1] + post.nijk[1])
  flag.j <- max(1, 1 - pre.nijk[2]) : min (vol$n.ijk[2], vol$n.ijk[2] + post.nijk[2])
  flag.k <- max(1, 1 - pre.nijk[3]) : min (vol$n.ijk[3], vol$n.ijk[3] + post.nijk[3])
  
  new.vol$vol3D.data[flag.i + pre.nijk[1], 
                     flag.j + pre.nijk[2], 
                     flag.k + pre.nijk[3]] <- vol$vol3D.data[flag.i, flag.j, flag.k] 
  if (any(!is.na(new.vol$vol3D.data))){
    new.vol$min.pixel <- min(new.vol$vol3D.data, na.rm = TRUE)
    new.vol$max.pixel <- max(new.vol$vol3D.data, na.rm = TRUE)
  } else {
    new.vol$min.pixel <- NA
    new.vol$max.pixel <- NA
  }	
  return(new.vol)
}