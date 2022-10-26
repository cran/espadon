#' @importFrom stats pnorm 
#' @export
#' 
#' 
.toy.rtdose.from.roi <- function (D.max, DSA, beam.nb, vol, struct, roi.name = NULL, roi.sname = NULL, 
                         roi.idx = NULL, alias = "", description = "") {
  
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  ct.ptv <- nesting.roi (vol, struct, roi.idx = roi.idx, 
                         xyz.margin = rep(4 * max(vol$dxyz),3))
  bin.PTV <- bin.from.roi(ct.ptv,struct, roi.sname = "ptv")
  bin.PTV.d <- bin.dilation(bin.PTV, radius = 2 * max(abs(vol$dxyz)))
  # PTV.xyz <- get.xyz.from.index(which(bin.PTV$vol3D.data),bin.PTV)
  PTV.G <- as.numeric(struct$roi.info[roi.idx, c ("Gx","Gy","Gz")])
  alpha <- atan(max(c(abs (as.numeric(struct$roi.info[roi.idx, c ("max.x","max.y","max.z")] - 
                                   PTV.G)) + abs(vol$dxyz),
                      abs (as.numeric(struct$roi.info[roi.idx, c ("min.x","min.y","min.z")] - 
                                   PTV.G)) + abs(vol$dxyz)))/DSA)*180/pi
  # PTV.size <- apply(PTV.xyz,2,max)  - apply(PTV.xyz,2,min) + bin.PTV$dxyz
  # alpha <- max(atan(PTV.size/DSA/2))*180/pi
  theta <- rev (seq (360, 0, length.out = beam.nb + 1)[-1] * pi / 180)
  
  direction <- cbind(-sin (theta), cos (theta), 0)
  orientation <- as.matrix(do.call (rbind,
                                    lapply(1:beam.nb, function(i) 
                                      c(c(0,0,1), vector.product (direction[i,],c(0,0,1))))))
  src <- sweep(-DSA *direction,2,PTV.G,FUN="+")
  
  src.extrem.pt.dist <- max(apply(src,1, function(sc)  
    max(.fnorm(sweep(t((vol$xyz.from.ijk %*% vol$cube.idx)[1:3,]),2, sc)))))
  dalpha <- 180 * min (vol$dxyz)/ src.extrem.pt.dist / pi
  alpha <- alpha + dalpha
  
  density <- vol.copy(vol)
  density$vol3D.data[density$vol3D.dat<=-900]  <- -1000
  density$vol3D.data[]  <- (density$vol3D.data[] / 1000 + 1)
  density$max.pixel <- (density$max.pixel / 1000 + 1)
  density$min.pixel <-(density$min.pixel / 1000 + 1)
  
  D.list <- lapply( 1:beam.nb, function(theta.idx){
    
    fan <- fan.beam (alpha =alpha,dalpha, origin=as.numeric(src[theta.idx,]),
                     orientation =  as.numeric(orientation[theta.idx,]), 
                     ref.pseudo = bin.PTV.d$ref.pseudo)
    vox <- fan.to.voxel(bin.PTV.d,fan,restrict = T )
    ray.index <- sort( unique(vox$ray.index))
    fan$xyz <- fan$xyz[ray.index, ]
    fan$local.coord <- fan$local.coord[ray.index, ]
    r <- as.numeric(DSA*sqrt((fan$xyz %*% fan$orientation[1:3])^2 + (fan$xyz %*% fan$orientation[4:6])^2))
    flou <- 1-pnorm(r,tan(alpha*pi/180) * DSA * 0.9 ,1)
    voxel <- .fan.to.voxel.with.att(density,fan, att=TRUE)
    voxel$att  <-  pnorm(voxel$att, 2, 5) * exp (-voxel$att / 150)
    vol.idx <- unique(voxel$vol.index)
    D_<- vol.copy(vol,modality ="rtdose")
    D_$vol3D.data[] <-0
    D_$vol3D.data <-  .mean_voxC(as.numeric(D_$vol3D.data), vol.idx-1, 
                                 voxel$vol.index-1,voxel$att, flou[voxel$ray.index])
    D_$vol3D.data  <- array(D_$vol3D.data ,dim=vol$n.ijk)
    D_$vol3D.data  <- D_$vol3D.data / max(D_$vol3D.data,na.rm=TRUE)
    D_$max.pixel <- 1
    D_$min.pixel <- 0
    # display.plane(top = D_, view.coord = PTV.G[3], struct=struct, interpolate = F)
    D_
  })
  
  D <- vol.copy(vol,modality ="rtdose" )
  D$vol3D.data[] <-0
  D$min.pixel <- 0
  D$max.pixel <- 0
  for(i in 1:beam.nb) D <- vol.sum(D,D.list[[i]],alias = alias, description = description)
  
  D$vol3D.data <- D.max*D$vol3D.data / nesting.roi (D, struct, 
                                                    roi.idx = roi.idx)$max.pixel
  D$max.pixel  <- max(D$vol3D.data,na.rm=TRUE)
  D$min.pixel  <- min(D$vol3D.data,na.rm=TRUE)
  return(D)
}

.toy.rtdose.from.bin <- function (D.max, DSA, beam.nb, vol, bin , alias = "", description = "") {
  
  PTV.xyz <- get.xyz.from.index(which(bin$vol3D.data),bin)
  margin <- max(abs(vol$dxyz))
  pt.min <- apply(PTV.xyz,2,min)
  pt.max <- apply(PTV.xyz,2,max)
  
  bin.PTV.d <- nesting.cube (bin, pt.min = pt.min - 4*margin,
                             pt.max =  pt.max+ 4*margin)
  
  PTV.G <-  apply(PTV.xyz,2,mean)
  alpha <- atan(max(c(abs (as.numeric(pt.max - PTV.G)) + abs(vol$dxyz),
                      abs (as.numeric(pt.min - PTV.G)) + abs(vol$dxyz)))/DSA)*180/pi
  theta <- rev (seq (360, 0, length.out = beam.nb + 1)[-1] * pi / 180)
  
  direction <- cbind(-sin (theta), cos (theta), 0)
  orientation <- as.matrix(do.call (rbind,
                                    lapply(1:beam.nb, function(i) 
                                      c(c(0,0,1), vector.product (direction[i,],c(0,0,1))))))
  src <- sweep(-DSA *direction,2,PTV.G,FUN="+")
  
  src.extrem.pt.dist <- max(apply(src,1, function(sc)  
    max (.fnorm(sweep(t((vol$xyz.from.ijk %*% vol$cube.idx)[1:3,]),2, sc)))))
  dalpha <- 180 * min (vol$dxyz)/ src.extrem.pt.dist / pi
  alpha <- alpha + dalpha
  
  density <- vol.copy(vol)
  density$vol3D.data[density$vol3D.dat<=-900]  <- -1000
  density$vol3D.data[]  <- (density$vol3D.data[] / 1000 + 1)
  density$max.pixel <- (density$max.pixel / 1000 + 1)
  density$min.pixel <-(density$min.pixel / 1000 + 1)
  
  D.list <- lapply( 1:beam.nb, function(theta.idx){
    
    fan <- fan.beam (alpha =alpha,dalpha, origin=as.numeric(src[theta.idx,]),
                     orientation =  as.numeric(orientation[theta.idx,]),
                     ref.pseudo = bin.PTV.d$ref.pseudo)
    vox <- fan.to.voxel(bin.PTV.d,fan,restrict = T )
    ray.index <- sort( unique(vox$ray.index))
    fan$xyz <- fan$xyz[ray.index, ]
    fan$local.coord <- fan$local.coord[ray.index, ]
    r <- as.numeric(DSA*sqrt((fan$xyz %*% fan$orientation[1:3])^2 + (fan$xyz %*% fan$orientation[4:6])^2))
    flou <- 1-pnorm(r,tan(alpha*pi/180) * DSA * 0.9 ,1)
    voxel <- .fan.to.voxel.with.att(density,fan, att=TRUE)
    voxel$att  <-  pnorm(voxel$att, 2, 5) * exp (-voxel$att / 150)
    vol.idx <- unique(voxel$vol.index)
    D_<- vol.copy(vol,modality ="rtdose")
    D_$vol3D.data[] <-0
    D_$vol3D.data <- .mean_voxC(as.numeric(D_$vol3D.data), vol.idx-1, 
                                voxel$vol.index-1,voxel$att, flou[voxel$ray.index])
    D_$vol3D.data  <- array(D_$vol3D.data ,dim=vol$n.ijk)
    D_$vol3D.data  <- D_$vol3D.data / max(D_$vol3D.data,na.rm=TRUE)
    D_$max.pixel <- 1
    D_$min.pixel <- 0
    # display.plane(top = D_, view.coord = PTV.G[3], interpolate = F)
    D_
  })
  
  D <- vol.copy(vol,modality ="rtdose" )
  D$vol3D.data[] <-0
  D$min.pixel <- 0
  D$max.pixel <- 0
  for(i in 1:beam.nb) D <- vol.sum(D,D.list[[i]],alias = alias, description = description)
  
  D$vol3D.data <- D.max*D$vol3D.data / nesting.cube (D, pt.min = pt.min,
                                                     pt.max =  pt.max)$max.pixel
  D$max.pixel  <- max(D$vol3D.data,na.rm=TRUE)
  D$min.pixel  <- min(D$vol3D.data,na.rm=TRUE)
  return(D)
}


# 
# .toy.rtdose <- function (CT.vol, PTV.radius, PTV.G, beam.nb) {
#   
#   sim.r <- PTV.radius + 5
#   D.tot <- CT.vol
#   D.tot$vol3D.data[] <- 0
#   D <- CT.vol
#   D$vol3D.data[] <- 0
#   N <- CT.vol
#   N$vol3D.data[] <- 0
#   DSA <- 600
#   theta <- rev (seq (360, 0, length.out = beam.nb + 1)[-1] * pi / 180)
#   
#   cone.b.a <- sim.r / DSA #angle du cone (= 1cm @ 1m)
#   
#   # distances min et max depuis le point de foc
#   DMAX <- DSA + 256 * sqrt (2)
#   d.theta <- min (CT.vol$dxyz) / DMAX / 1.5 # angle de simulation pour 1mm en distal
#   uv <- as.matrix (expand.grid (seq (-cone.b.a, cone.b.a, by=d.theta),
#                                 seq (-cone.b.a, cone.b.a, by=d.theta)))
#   r <- sqrt(uv[, 1]^2 + uv[, 2]^2) * DSA
#   transverse.int <- 1 - pnorm (r, (PTV.radius + sim.r)/2, 1)
#   
#   uvt <- cbind (uv[r <= cone.b.a * DSA, ], transverse.int[r <= cone.b.a * DSA])
#   colnames (uvt) <- c("u", "v", "t")
#   
#   
#   for (angle in theta) {
#     
#     # point de foc
#     foc <- PTV.G + DSA * c(sin (angle), -cos (angle), 0)
#     
#     PTV2foc <- foc - PTV.G
#     PTV2foc <- PTV2foc / sqrt (sum (PTV2foc^2))
#     
#     n.b <- -PTV2foc # base du faisceau
#     v.b <- c(n.b[2], n.b[3], n.b[1])
#     u.b <- vector.product (n.b, v.b)
#     u.b <- u.b / sqrt (sum (u.b^2))
#     v.b <- vector.product (n.b, u.b)
#     
#     d.mM <- Re (polyroot (c (sum (foc^2) - 2 * 128^2, 2 * sum (foc * n.b), sum (n.b^2))))
#     d.min <- d.mM[1]
#     d.max <- d.mM[2]
#     
#     b.dirs <- t (apply (uvt, 1, function (uvt) {
#       n.beamlet <- n.b + uvt[1] * u.b + uvt[2] * v.b
#       c (n.beamlet / sqrt (sum (n.beamlet^2)), uvt)
#     }))
#     
#     beamlets <- apply (b.dirs, 1, function (bl) {
#       ds <- .95 * min (CT.vol$dxyz)
#       L <- get.line (CT.vol, foc, bl[1:3], seq (d.min, d.max, by=ds))
#       idx <- range (which (L$value > -900))
#       L <- L[idx[1]:idx[2], ]
#       L$mu <- L$value / 1000 + 1
#       L$s_ <- cumsum (L$mu * ds)
#       D <- pnorm(L$s_, 2, 5) * exp (-L$s_ / 150)
#       L$D <- D / max (D) * bl[6]
#       L$PTV <- (L$x - PTV.G[1])^2 + (L$y - PTV.G[2])^2 +
#         (L$z - PTV.G[3])^2 <= PTV.radius^2
#       return (L)
#     })
#     
#     D$vol3D.data[] <- 0
#     N$vol3D.data[] <- 0
#     for (bl in beamlets) {
#       ijk <- round (get.ijk.from.xyz (cbind (bl$x, bl$y, bl$z), CT.vol)) + 1
#       D$vol3D.data[ijk] <- D$vol3D.data[ijk] + bl$D
#       N$vol3D.data[ijk] <- N$vol3D.data[ijk] + 1
#     }
#     avg <- which (N$vol3D.data != 0)
#     D$vol3D.data[avg] <- D$vol3D.data[avg] / N$vol3D.data[avg]
#     D.tot$vol3D.data[avg] <- D.tot$vol3D.data[avg] + D$vol3D.data[avg]
#     
#   }
#   
#   D.tot$vol3D.data <- D.tot$vol3D.data / max (D.tot$vol3D.data) * 52
#   D.tot$min.pixel <- 0
#   D.tot$max.pixel <- 52
#   
#   D.tot$modality <- "rtdose"
#   D.tot$description <- "IMRT|PTV52"
#   
#   return (D.tot)
# }