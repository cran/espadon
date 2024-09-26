#' @name bin.cuboid
#' @aliases bin.cylinder
#' @aliases bin.ellipsoid

#' @title Create a volume containing predefined shapes
#' @description These functions create espadon objects of class “volume”, and 
#' of modality “binary” or “weight”, by selecting the voxels defining a rectangular 
#' cuboid, an elliptical cylinder or an ellipsoid.
#' @param back.vol "volume" class object.
#' @param pt.min,pt.max Minimum and maximum x, y, z coordinates of the vertices 
#' of the rectangular cuboid.
#' @param center Numeric vector of length 3, representing the xyz-center of the 
#' cylinder or the ellipsoid.
#' @param radius Positive number, or xy-vector or xyz-vector of 2 or 3 positive 
#' numbers, representing the radius of the cylinder or the ellipsoid.
#' @param height Positive number representing the height of the cylinder.
#' @param ref.shape Character string. Pseudonym of the frame of reference in which 
#' the requested shape is defined.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If\code{T.MAT = NULL}, \code{ref.shape} must be 
#' equal to \code{back.vol$ref.pseudo}.
#' @param modality modality ("binary" or "weight") of the generated object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. 
#' @param ... Additional arguments, such as \code{error}, representing the relative 
#' error on the volume of the ellipsoid, if \code{back.vol$orientation = c(1, 0, 0, 0, 1, 0)}.
#' @return Returns a "volume" class object of "binary" or "weight" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{back.vol}. 
#' \itemize{
#' \item In the “binary” modality, voxels with 50 percent of their volume within  
#' the requested shape are set to \code{TRUE}.
#' \item In the “weight” modality, the value of each voxel is its volume fraction 
#' included in the requested shape. 
#' }
#' @seealso \link[espadon]{add.shape}
#' @examples
#' # Creation of back.vol
#' CT <- vol.create (c(80, 80,40), c(1.2, 1.2, 2), 
#'                   pt000 = c(-50.4,-50.4,-75), modality = "ct", 
#'                   default.value = as.integer(-997), value.sd = 1)
#'                    
#' # Creation of a cuboid
#' cuboid <- bin.cuboid(CT, pt.min = c(-45.3, -40.2, -60.4), 
#'                      pt.max =  c(-15.6, -20.2, -40.2))
#'                        
#' # Definition of the shape frame of reference
#' angle.yz =30
#' Myz <- cbind(c(1, 0, 0, 0),
#'              c(0, cos(pi * angle.yz / 180), sin(pi * angle.yz / 180), 0) ,
#'              c(0, -sin(pi * angle.yz / 180), cos(pi * angle.yz / 180), 0),
#'              c(0, 0, 0, 1 )) 
#' angle.xy <- 59
#' Mxy <- cbind(c(cos(pi * angle.xy / 180), sin(pi * angle.xy / 180), 0, 0) ,
#'              c(-sin(pi * angle.xy / 180), cos(pi * angle.xy / 180), 0, 0), 
#'              c(0 ,0, 1, 0), c(0, 0, 0, 1)) 
#'
#' t.mat <- ref.srctodest.add("ref1", "shaperef", TM = Myz %*% Mxy) 
#' 
#' # Creation of a cylinder 
#' cylinder <- bin.cylinder(CT, center = c(10.3, 30.6, -30.7), 
#'                          radius =  c(10, 20), height = 50, 
#'                          ref.shape ="shaperef", T.MAT = t.mat) 
#'                       
#' # Creation of an ellipsoid 
#' ellipsoid <- bin.ellipsoid(CT, center = c(-20.1, 0.1, -15), 
#'                            radius =  c(19.3, 20.2, 15.3))  
#'                        
#'  # Display                         
#'  k.idx <- unique(which(cuboid$vol3D.data>0, arr.ind = TRUE)[,3]) - 1
#'  display.3D.stack(cuboid, k.idx, border = FALSE,
#'                   col = c("#FFFFFF00", "#EBDFDFFF", "#D8BFBFFF", "#C59F9FFF", 
#'                           "#B27F7FFF", "#9F5F5FFF", "#8C3F3FFF", "#791F1FFF", 
#'                           "#660000FF"))   
#'                            
#'  k.idx <- unique(which(cylinder$vol3D.data>0, arr.ind = TRUE)[,3]) - 1
#'  display.3D.stack(cylinder, k.idx, border = FALSE,
#'                   col = c("#FFFFFF00", "#DFEBDFFF", "#BFD8BFFF", "#9FC59FFF", 
#'                           "#7FB27FFF", "#5F9F5FFF", "#3F8C3FFF", "#1F791FFF", 
#'                           "#006600FF"))    
#'
#'  k.idx <- unique(which(ellipsoid$vol3D.data>0, arr.ind = TRUE)[,3]) - 1
#'  display.3D.stack(ellipsoid, k.idx, border = FALSE,
#'                   col = c("#FFFFFF00", "#DFDFEBFF", "#BFBFD8FF", "#9F9FC5FF", 
#'                           "#7F7FB2FF", "#5F5F9FFF", "#3F3F8CFF", "#1F1F79FF", 
#'                           "#000066FF"))                           
#' @export
bin.cuboid <- function(back.vol, pt.min, pt.max, ref.shape = back.vol$ref.pseudo, T.MAT = NULL,
                       modality ="weight", alias = "", description = NULL){
  if (length(pt.min)!=3 | length(pt.max)!=3) stop ("pt.min and pt.max must be numerical vectors of length 3.")
  args <- list()
  args[["back.vol"]] <- back.vol
  
  rg <- apply(cbind(pt.min,pt.max), 1,range)
  side <- args[["pt.max"]] - args[["pt.min"]]
  if (any(side<1e-6)) stop("the sides of the cuboid must not be equal to zero")
  args[["pt.min"]] <- rg[1,]
  args[["pt.max"]] <- rg[2,]
  args[["ref.shape"]] <- ref.shape
  args[["T.MAT"]] <- T.MAT
  args[["modality"]] <- modality[1]
  args[["alias"]] <- alias
  args[["description"]] <- description
  args[["shape"]] <- "cuboid"
  if (is.null(description)) {
    
    if (length(unique(side))==1) {description <- paste(side[1], "mm square cube")
    } else {description <- paste0("rectangular cuboid with (", paste(side, collapse=","), ") mm sides")}
  }
  do.call(bin.shape, args)
}

#' @rdname bin.cuboid   
#' @export
bin.cylinder <- function(back.vol, center, radius, height, ref.shape = back.vol$ref.pseudo, T.MAT = NULL,
                         modality ="weight", alias = "", description = NULL){
  if (length(radius)>2) stop ("radius must be a numerical vector of length 1 or 2.")
  args <- list()
  args[["back.vol"]] <- back.vol
  args[["center"]] <- center
  args[["radius"]] <- abs(rep(radius,2)[1:2])
  args[["height"]] <- abs(height[1])
  args[["ref.shape"]] <- ref.shape
  args[["T.MAT"]] <- T.MAT
  args[["modality"]] <- modality[1]
  args[["alias"]] <- alias
  args[["description"]] <- description
  args[["shape"]] <- "cylinder"
  if (is.null(description)) {
    if (length(unique(args[["radius"]] ))==1) {description <- paste("cylinder with ",  args[["radius"]][1], 
                                                                    "mm radius and ",  args[["height"]], "mm height")
    } else {description <- paste0("cylinder with radius (",  paste(args[["radius"]], collapse=","), ") mm and height ",
                                  args[["height"]], " mm" )}
  }
  do.call(bin.shape, args)
}

#' @rdname bin.cuboid   
#' @export
bin.ellipsoid <- function(back.vol, center, radius, ref.shape = back.vol$ref.pseudo, T.MAT = NULL,
                          modality ="weight", alias = "", description = NULL,...){
  if (length(radius)>3) stop ("radius must be a numerical vector of length less than 3.")
  args <- tryCatch(list(...), error = function(e)list())
  args[["back.vol"]] <- back.vol
  args[["center"]] <- center
  args[["radius"]] <- abs(c(radius,rep(radius[length(radius)],2)))[1:3]
  args[["ref.shape"]] <- ref.shape
  args[["T.MAT"]] <- T.MAT
  args[["modality"]] <- modality[1]
  args[["alias"]] <- alias
  if (is.null(description)) {
    if (length(unique(args[["radius"]] ))==1) {description <- paste(args[["radius"]][1], "mm radius sphere")
    } else {description <- paste0("ellipsoid of radius (",  paste(args[["radius"]], collapse=","), ") mm")}
  }
  args[["description"]] <- description
  args[["shape"]] <- "ellipsoid"
  do.call(bin.shape, args)
}

bin.shape <- function(back.vol, shape, ref.shape = back.vol$ref.pseudo, T.MAT = NULL,
                      modality ="weight", alias = "", description = NULL, ...) {
  args <- tryCatch(list(...), error = function(e)list())
  
  if (!is (back.vol, "volume")) 
    stop ("back.vol should be a volume class object.", call.=FALSE)
  
  rigid.M <- get.rigid.M (T.MAT, src.ref = ref.shape, dest.ref = back.vol$ref.pseudo)
  if (is.null(rigid.M)) 
    stop("back.vol and the box do not have the same frame of reference: define a valid T.MAT", call.=FALSE)
  
  
  if (is.null(description)) description <- args[["shape"]]
  vol.out <- vol.copy (back.vol, alias = alias, modality=modality, description=description)
  vol.out$vol3D.data[] <- 0
  vol.out$min.pixel <- 0
  vol.out$max.pixel <- 0
  
  #2D
  idx.c <- which(apply(abs(back.vol$xyz.from.ijk[1:3,1:3]),2,sum)==0)
  idx.r <-  which(apply(abs(back.vol$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    if (abs(back.vol$xyz0[1,idx.r])>1e-6) return(NULL)
    u <- back.vol$xyz.from.ijk
    u[idx.r,idx.c]<- 1
    ijk.from.shape <- solve(u)
    ijk.from.shape[idx.r,idx.c] <- 0
  } else {#3D
    ijk.from.shape <- solve(back.vol$xyz.from.ijk)}
  
  ijk.from.shape <- ijk.from.shape %*% rigid.M
  shape.from.ijk <- get.rigid.M (T.MAT, src.ref=back.vol$ref.pseudo, dest.ref = ref.shape) %*% back.vol$xyz.from.ijk
  
  # on crée un volume avec un pas au plus proche de back.vol
  p <- abs(shape.from.ijk %*% cbind(c(1,0,0,0),c(0,1,0,0), c(0,0,1,0)))
  Vb.dxyz <- apply(p[1:3,],1,max)
  Vb.pt000 <- (apply(shape.from.ijk %*% back.vol$cube.idx,1,range))[,1:3]
  Vb.n.ijk <- ceiling((Vb.pt000[2,]- Vb.pt000[1,])/Vb.dxyz + 1)
  Vb <- vol.create(n.ijk = Vb.n.ijk, dxyz = Vb.dxyz, modality ="weight",pt000 = Vb.pt000[1,],
                   default.value = 0,ref.pseudo = ref.shape)
  
  if (shape == "cuboid"){
    pt.min <-  args[["pt.min"]]
    pt.max <-  args[["pt.max"]]
    base_m <- matrix(c(pt.min, pt.min[1],pt.max[2],pt.min[3], pt.max[1:2],
                       pt.min[3],pt.max[1],pt.min[2:3],pt.min), ncol = 3, byrow = TRUE)
    base_M <- base_m
    base_M[,3] <- pt.max[3]
    # on regarde à quelle indice de Vb correpond les points du struct
    Mat <- solve (Vb$xyz.from.ijk)
    
    base_m.ijk <- cbind(base_m,1) %*% t(Mat)
    base_M.ijk <- cbind(base_M,1) %*% t(Mat)
    
    range.s <- t(apply(rbind(base_m.ijk[,1:3],base_M.ijk[,1:3]),2,range))
    range.s <- cbind(floor(range.s[,1]),ceiling(range.s[,2]))
    vol3D.range  <- cbind( c(range.s[1,1]-1,max(range.s[2,1]-1,-1),max(range.s[3,1],0)),  
                           apply(cbind(range.s[,2]+1, Vb$n.ijk - 1),1,min))
    common.range <- vol3D.range
    common.range[,1] <-apply(cbind(range.s[,1], c(0,0,0)),1,max)
    common.range[,2] <- apply(cbind(range.s[,2], Vb$n.ijk - 1),1,min)
    if (any(common.range[,1]>common.range[,2])) return (vol.out)
    
    vol3D.nijk <-vol3D.range[,2]-vol3D.range[,1] +1
    vol3D <- array(0,dim = vol3D.nijk)
    base_ijk <- base_m.ijk[, 1:3]
    base_ijk[,3] <- vol3D.range[3,1]
    u <- base_ijk[2:5,1:3]-base_ijk[1:4,1:3]
    dist.u<- sqrt(u[,1]^2 +u[,2]^2 + u[,3]^2)
    u[ ,1] <- u[ ,1]/dist.u; 
    u[ ,2] <- u[ ,2]/dist.u;
    pO <- sweep(base_ijk[,1:3] ,2,vol3D.range[,1])
    vol3D_ <- array(.ouline2voxC(as.numeric(t(u)),c(vol3D.nijk[1:2],1) , 0,  0, as.numeric(t(pO[1:4, ])),
                                 lambda_max = dist.u, ntab = sum(ceiling(dist.u)+5)),dim=c(vol3D.nijk[1:2],1))
    vol3D_[vol3D_>1] <- 1
    vol3D_[vol3D_<0] <- 0
    rg.k <- vol3D.range[3,1]:vol3D.range[3,2]
    k.min <- floor(base_m.ijk[1,3] + 0.5)
    k.max <- floor(base_M.ijk[1,3] + 0.5)
    vol3D[,,(k.min:k.max)-vol3D.range[3,1] + 1] <- vol3D_[,,rep(1,k.max-k.min + 1)]
    
    if (k.min == k.max){
      vol3D[,,k.min-vol3D.range[3,1] + 1] <- vol3D[,,k.min-vol3D.range[3,1] + 1] * (base_M.ijk[1,3]- base_m.ijk[1,3])
    }else{
      vol3D[,,k.min-vol3D.range[3,1] + 1] <- vol3D[,,k.min-vol3D.range[3,1] + 1] * (k.min + 0.5 - base_m.ijk[1,3])
      vol3D[,,k.max-vol3D.range[3,1] + 1] <- vol3D[,,k.max-vol3D.range[3,1] + 1] * (base_M.ijk[1,3]  - k.max + 0.5)
    }
  } else if (shape == "cylinder") {
    
    Rijk <- args[["radius"]]/Vb$dxyz[1:2]
    H <- args[["height"]]/Vb$dxyz[3]
    Cijk <-get.ijk.from.xyz(args[["center"]],Vb)
    x.breaks <- floor(Cijk[1]-Rijk[1]):(ceiling(Cijk[1]+Rijk[1]) +1)  - 0.5 #0:Vb$n.ijk[1] - 0.5
    y.breaks <- floor(Cijk[2]-Rijk[2]):(ceiling(Cijk[2]+Rijk[2]) +1)  - 0.5
    
    costheta <- (x.breaks-Cijk[1])/Rijk[1]
    sintheta <- (y.breaks-Cijk[2])/Rijk[2]
    fx <- abs(costheta) <=1
    fy <- abs(sintheta) <=1
    costheta <- costheta[fx]
    x.breaks <- x.breaks[fx]
    x.theta <- acos(costheta)
    x.theta <- (c(2*pi-x.theta,x.theta) + 2*pi) %% (2*pi)
    sintheta <- sintheta[fy]
    y.breaks <- y.breaks[fy]
    y.theta <- asin(sintheta)
    y.theta <- (c(pi-y.theta,y.theta) + 2*pi) %% (2*pi)
    theta <- sort(unique (c(x.theta, y.theta)),decreasing = TRUE)
    theta <- c(theta,theta[1]-2*pi)
    
    
    vol3D.range  <- cbind( c(x.breaks[1]-1.5, max(y.breaks[1]-1.5,-1),floor(Cijk[3]-H/2) -1),  
                           apply(cbind(c(x.breaks[length(x.breaks)]+0.5,
                                         y.breaks[length(y.breaks)]+0.5,
                                         ceiling(Cijk[3]+H/2)+1),  Vb$n.ijk - 1),1,min))

    
    common.range <- vol3D.range
    common.range[,1] <-apply(cbind(common.range[,1], c(0,0,0)),1,max)
    common.range[,2] <- apply(cbind(common.range[,2], Vb$n.ijk - 1),1,min)
    le <- length(theta)
    if (any(common.range[,1]>common.range[,2]) | le==1) return (vol.out)
    vol3D.nijk <-vol3D.range[,2]-vol3D.range[,1] +1
    vol3D <- array(0,dim = vol3D.nijk)
    
    
    O_ijk <-cbind(Rijk[1]*cos(theta) + Cijk[1] - vol3D.range[1,1], 
                  Rijk[2]*sin(theta) + Cijk[2] - vol3D.range[2,1], 
                  0)
    mid.theta <- 0.5*(theta[1:(le-1)] + theta[2:le])
    pt <- floor(cbind(Rijk[1]*cos(mid.theta) + Cijk[1] - vol3D.range[1,1], 
                      Rijk[2]*sin(mid.theta) + Cijk[2] - vol3D.range[2,1], 
                      0)+0.5) + 1
    u <- O_ijk[2:le, ] - O_ijk[1:(le-1), ]
    dist.u<- sqrt(u[,1]^2 +u[,2]^2 + u[,3]^2)
    u[ ,1] <- u[ ,1]/dist.u; 
    u[ ,2] <- u[ ,2]/dist.u;
    
    vol3D_ <- array(.ouline2voxC(as.numeric(t( u)),c(vol3D.nijk[1:2],1) , 0,  0, as.numeric(t(O_ijk[1:(le-1), ])),
                                 lambda_max = dist.u, ntab = sum(ceiling(dist.u)+5)),dim=c(vol3D.nijk[1:2],1))
    
    theta_ <- theta
    f <- cos(theta)!=0
    theta_[f] <- atan(Rijk[2] * tan(theta[f])/ Rijk[1])
    
    S.triangle <- 0.5*abs(Rijk[1] * Rijk[2] * sin(diff(theta))) 
    zone <- theta /(pi/2)
    f <- round(zone) == zone
    S.ellipse <- (as.integer(zone) + (as.integer(zone) %% 2)) * Rijk[1]*Rijk[2]*pi * 0.25
    S.ellipse[!f] <- S.ellipse[!f]  + 0.5*Rijk[1]*Rijk[2]* atan(tan(theta_[!f])*Rijk[1]/Rijk[2])
    
    pt <- cbind(pt,abs(diff(S.ellipse)) - S.triangle)
    tab<- as.matrix(do.call (rbind,by(pt, pt[,1]+ 1i*pt[,2],function(tab) cbind(tab[1,1:3],sum(tab[,4])))))
    f <- (tab[,3] == 1) & (tab[,2] >0) & (tab[,1] >0) &  (tab[,2] <= vol3D.nijk[2]) & (tab[,1] <= vol3D.nijk[1])
    vol3D_[tab[f,1:3 ]] <- vol3D_[tab[f,1:3]] + tab[f,4]
    vol3D_[vol3D_>1] <- 1
    vol3D_[vol3D_<0] <- 0  
    
    rg.k <- vol3D.range[3,1]:vol3D.range[3,2]
    k.min <- floor(Cijk[3]-H/2 + 0.5)
    k.max <- floor(Cijk[3]+H/2 + 0.5)
    vol3D[,,(k.min:k.max)-vol3D.range[3,1] + 1] <- vol3D_[,,rep(1,k.max-k.min + 1)]
    
    if (k.min == k.max){
      vol3D[,,k.min-vol3D.range[3,1] + 1] <- vol3D[,,k.min-vol3D.range[3,1] + 1] * H
    }else{
      vol3D[,,k.min-vol3D.range[3,1] + 1] <- vol3D[,,k.min-vol3D.range[3,1] + 1] * (k.min + 0.5 - (Cijk[3]-H/2))
      vol3D[,,k.max-vol3D.range[3,1] + 1] <- vol3D[,,k.max-vol3D.range[3,1] + 1] * (Cijk[3]+H/2  - k.max + 0.5)
    }
  } else if (shape == "ellipsoid"){
    #------------------------------
    
    Rijk <-  args[["radius"]]/Vb$dxyz
    Cijk <-get.ijk.from.xyz(args[["center"]],Vb)
    
    
    k.min <- floor(Cijk[3] - Rijk[3])
    k.max <- ceiling(Cijk[3] + Rijk[3])
    samp.nb <- -1
    error_abs <- +Inf
    error_rel <- +Inf
    
    error_rel0 <- args[["error"]]
    
    if (!is.numeric(error_rel0)) error_rel0 <- 0.01
    error_rel0 <- abs(error_rel0[1])
    error_abs0 <- error_rel0/10/( pi*Rijk[1]*Rijk[2])
    while((abs(error_rel) > error_rel0) & (abs(error_abs*pi*Rijk[1]*Rijk[2])> error_abs0) & ((k.max - k.min+1) * samp.nb < 500)){
      samp.nb <- samp.nb + 2
      mid.samp <- floor(0.5*samp.nb)
      rg.k <- seq(k.min - mid.samp/samp.nb,k.max + mid.samp/samp.nb,1/samp.nb)
      sinphi <- (rg.k - Cijk[3])/Rijk[3]
      fz <- abs(sinphi) < 1 
      error_abs <- sum(1-sinphi[fz]^2)/samp.nb - (4*Rijk[3]/3)
      error_rel <- 100 * error_abs /(4*Rijk[3]/3) 
    }
    
    cosphi <- suppressWarnings(sqrt(1-sinphi^2))
    
    x.breaks <- floor(Cijk[1] - Rijk[1]):(ceiling(Cijk[1]+Rijk[1]) +1)  - 0.5 #0:Vb$n.ijk[1] - 0.5
    y.breaks <- floor(Cijk[2] - Rijk[2]):(ceiling(Cijk[2]+Rijk[2]) +1)  - 0.5
    
    vol3D.range <- rbind(c(x.breaks[1]+0.5, x.breaks[length(x.breaks)]-0.5), 
                         c(y.breaks[1]+0.5, y.breaks[length(y.breaks)]-0.5),c(k.min,k.max))
    
    common.range <- vol3D.range
    common.range[,1] <-apply(cbind(common.range[,1], c(0,0,0)),1,max)
    common.range[,2] <- apply(cbind(common.range[,2], Vb$n.ijk - 1),1,min)
    if (any(common.range[,1]>common.range[,2])) return (vol.out)
    
    
    vol3D.nijk <-vol3D.range[,2]-vol3D.range[,1] +1
    vol3D.nijk[3] <- vol3D.nijk[3] * samp.nb
    vol3D <- array(0,dim = vol3D.nijk)
    
    O_ijk <- dist.u <- u <- list()
    for (cp.idx in which(fz)){
      
      costheta <- (x.breaks-Cijk[1])/Rijk[1]/cosphi[cp.idx]
      sintheta <- (y.breaks-Cijk[2])/Rijk[2]/cosphi[cp.idx]
      x.theta <- acos(costheta[abs(costheta)<=1])
      x.theta <- (c(2*pi-x.theta,x.theta) + 2*pi) %% (2*pi)
      y.theta <- asin(sintheta[abs(sintheta)<=1])
      y.theta <- (c(pi-y.theta,y.theta) + 2*pi) %% (2*pi)
      theta <- sort(unique (c(x.theta, y.theta)),decreasing = TRUE)
      theta <- c(theta,theta[1]-2*pi)
      le <- length(theta)
      
      if (le>1){
        O_ijk[[cp.idx]] <-cbind(Rijk[1]*cosphi[cp.idx]*cos(theta) + Cijk[1] - vol3D.range[1,1], 
                                Rijk[2]*cosphi[cp.idx]*sin(theta) + Cijk[2] - vol3D.range[2,1], 
                                cp.idx -1)
        
        mid.theta <- 0.5*(theta[1:(le-1)] + theta[2:le])
        pt <- floor(cbind(Rijk[1]*cosphi[cp.idx]*cos(mid.theta) + Cijk[1] - vol3D.range[1,1], 
                          Rijk[2]*cosphi[cp.idx]*sin(mid.theta) + Cijk[2]- vol3D.range[2,1], 
                          cp.idx-1)+0.5) + 1
        
        
        u[[cp.idx]] <- O_ijk[[cp.idx]][2:le, ] - O_ijk[[cp.idx]][1:(le-1), ]
        dist.u[[cp.idx]] <- sqrt(u[[cp.idx]][,1]^2 +u[[cp.idx]][,2]^2 + u[[cp.idx]][,3]^2)
        u[[cp.idx]][ ,1] <- u[[cp.idx]][ ,1]/dist.u[[cp.idx]]; 
        u[[cp.idx]][ ,2] <- u[[cp.idx]][ ,2]/dist.u[[cp.idx]];
        
        theta_ <- theta
        f <- cos(theta)!=0
        theta_[f] <- atan(Rijk[2] * tan(theta[f])/ Rijk[1] ) 
        
        S.triangle <- 0.5*abs(Rijk[1] * Rijk[2] * sin(diff(theta))) 
        zone <- theta /(pi/2)
        f <- round(zone) == zone
        S.ellipse <- (as.integer(zone) + (as.integer(zone) %% 2)) * Rijk[1]*Rijk[2]*pi * 0.25
        S.ellipse[!f] <- S.ellipse[!f]  + 0.5*Rijk[1]*Rijk[2]* atan(tan(theta_[!f])*Rijk[1]/Rijk[2]) 
        pt <- cbind(pt,(abs(diff(S.ellipse)) - S.triangle) * cosphi[cp.idx]^2)
        tab<- as.matrix(do.call (rbind,by(pt, pt[,1]+ 1i*pt[,2],function(tab) cbind(tab[1,1:3],sum(tab[,4])))))
        vol3D [tab[,1:3]] <-tab[,4]
        
        O_ijk[[cp.idx]] <- O_ijk[[cp.idx]][1:(le-1),]
      # vol3D <- vol3D +
      #   array(.ouline2voxC(as.numeric(t( u[[cp.idx]])),vol3D.nijk , 0:(vol3D.nijk[3]-1),
      #                                0:(vol3D.nijk[3]-1),as.numeric(t( O_ijk[[cp.idx]] )),
      #                                lambda_max = dist.u[[cp.idx]], ntab = sum(ceiling(dist.u[[cp.idx]])+5)),dim=vol3D.nijk)
      #  print(sum(vol3D[,,O_ijk[[cp.idx]][1,3] + 1]) - pi*Rijk[1]*Rijk[2]* cosphi[cp.idx]^2)
      }
    }
    dist.u <- unlist(dist.u)
    O_ijk <- do.call(rbind,O_ijk)
    u <- do.call(rbind,u)
    
    vol3D <- vol3D + array(.ouline2voxC(as.numeric(t( u)),vol3D.nijk , 0:(vol3D.nijk[3]-1), 
                                        0:(vol3D.nijk[3]-1),as.numeric(t(O_ijk)),
                                        lambda_max = dist.u, ntab = sum(ceiling(dist.u)+5)),dim=vol3D.nijk)
    
    
    sm <- 1:(k.max-k.min + 1)*samp.nb- samp.nb + 1
    for(idx in (0:(samp.nb-1))[-1]) vol3D[,,sm] <- vol3D[,,sm] + vol3D[,,sm+idx]
    vol3D <- vol3D[,,sm]/samp.nb
    
  }
  
  
  rgi <- (common.range[1,1]:common.range[1,2])
  rgj <- (common.range[2,1]:common.range[2,2])
  rgk <- (common.range[3,1]:common.range[3,2]) 
  Vb$vol3D.data[ rgi + 1, rgj + 1, rgk + 1] <- 
    vol3D[ rgi-vol3D.range[1,1]+1, rgj-vol3D.range[2,1]+1, rgk-vol3D.range[3,1]+1]
  
  # Vb$max.pixel <- max(Vb$vol3D.data, na.rm=TRUE)
  # Vb$min.pixel <- min(Vb$vol3D.data, na.rm=TRUE)
  #   display.3D.stack(Vb,Vb$k.idx, border =F)
  # rgl::bg3d("black")
  
  if (grid.equal (vol.out, Vb)) {
    vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1] <- Vb$vol3D.data[rgi + 1, rgj + 1, rgk + 1]
    vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1][vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1] >1] <- 1
    vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1][vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1] <0] <- 0
    vol.out$max.pixel <- max(vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1], na.rm=TRUE)
    vol.out$min.pixel <- min(vol.out$vol3D.data[rgi + 1, rgj + 1, rgk + 1], na.rm=TRUE)
    
  } else {
    f <-common.range[,1]>1 
    common.range[f,1] <- common.range[f,1]-1
    f <-common.range[,2]< Vb$n.ijk -1 
    common.range[f,2] <- common.range[f,2]+1
    cube.idx <- Vb$cube.idx
    cube.idx [1,] <- common.range[1,c(1,2,2,1,1,2,2,1)]
    cube.idx [2,] <- common.range[2,c(1,1,2,2,1,1,2,2)]
    cube.idx [3,] <- common.range[3,c(1,1,1,1,2,2,2,2)]
    rg.cube <- ijk.from.shape %*% Vb$xyz.from.ijk %*% cube.idx
    rgi <- min(floor(rg.cube[1,])):max(ceiling(rg.cube[1,]))
    rgj <- min(floor(rg.cube[2,])):max(ceiling(rg.cube[2,]))
    rgk <- min(floor(rg.cube[3,])):max(ceiling(rg.cube[3,]))
    rgi <- rgi[rgi>=0 & rgi<vol.out$n.ijk[1]]
    rgj <- rgj[rgj>=0 & rgj<vol.out$n.ijk[2]]
    back.rgk.loc <- match (rgk, vol.out$k.idx)
    fna <- is.na(back.rgk.loc)
    rgk <- rgk[!fna]
    back.rgk.loc <-  back.rgk.loc[!fna]
    
    ijk.out <- as.matrix(expand.grid(rgi,rgj,rgk,1))
    R.ijk.out <- as.matrix(expand.grid(rgi+1,rgj+1,back.rgk.loc))
    Mref <- t(solve(Vb$xyz.from.ijk) %*% shape.from.ijk)
    ijk <- (ijk.out %*% Mref)[,1:3,drop =FALSE]
    
    s_ijk <-apply(abs(rbind(c(1,0,0,0),c(0,1,0,0), c(0,0,1,0)) %*% Mref)[,1:3, drop =FALSE],1,max)
    vol.out$vol3D.data[R.ijk.out] <- .getvaluefromijkC (vol3D = as.numeric(Vb$vol3D.data),
                                                        interpolate = TRUE,
                                                        i = as.numeric(ijk[ ,1]),
                                                        j = as.numeric(ijk[ ,2]),
                                                        k = as.numeric(ijk[ ,3]),
                                                        k_idx = Vb$k.idx,
                                                        k_loc = Vb$k.idx, n_ijk=Vb$n.ijk,
                                                        s_ijk = s_ijk)
    vol.out$vol3D.data[R.ijk.out][is.na(vol.out$vol3D.data[R.ijk.out])] <- 0
    vol.out$vol3D.data[R.ijk.out][vol.out$vol3D.data[R.ijk.out] >1] <- 1
    vol.out$vol3D.data[R.ijk.out][vol.out$vol3D.data[R.ijk.out] <0] <- 0
    vol.out$max.pixel <- max(vol.out$vol3D.data[R.ijk.out], na.rm=TRUE)
    vol.out$min.pixel <- min(vol.out$vol3D.data[R.ijk.out], na.rm=TRUE)
    
  }
  
  if (modality == "binary"){
    vol.out$max.pixel <- !as.logical(round(1-vol.out$max.pixel))
    vol.out$min.pixel <- !as.logical(round(1-vol.out$min.pixel))
    vol.out$vol3D.data <- array(!as.logical (round(1-vol.out$vol3D.data)),dim =vol.out$n.ijk)
  }
  return(vol.out)
}