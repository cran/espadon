####################################################################
#' Creation of a binary volume according to RoI
#' @description The \code{bin.from.roi} function creates a "volume" class object,
#' of modality "binary" or "weight", by selecting the voxels defined by the RoI.
#' @param vol "volume" class object.
#' @param struct "struct" class object.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object.
#' By default \code{roi.idx = NULL}. See Details.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{struct$ref.pseudo} must be equal to 
#' \code{vol$ref.pseudo} or set to \code{NULL}.
#' @param within Boolean, defaults to \code{TRUE}. If \code{within = TRUE}, the
#' contours included in a RoI are managed,
#' depending on their \code{$level} field. If \code{within = FALSE}, only the 
#' \code{$level = 0} fields of the RoI are used (i.e. only the external outlines).
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL}
#' (default value), it will be set to \code{struct$roi.info$roi.pseudo[roi.idx]}.
#' @param modality modality ("binary" or "weight") of the generated object.
#' @param ... additional argument. 
#' @return Returns a "volume" class object of "binary" or "weight" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol}. 
#' \itemize{
#' \item In the “binary” modality, voxels with 50 percent of their volume within  
#' the ROI are set to \code{TRUE}.
#' \item In the “weight” modality, the value of each voxel is its volume fraction 
#' included in the ROI. 
#' }
#' 
#' @details \code{roi.name}, \code{roi.sname}, and \code{roi.idx} must select
#' only one RoI.
#' @seealso \link[espadon]{bin.from.vol}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 3
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name =  c("eye", "optical nerve", "brain"), 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' # "optical nerve" binary, with  modality "binary" and "weight"
#' binl <- bin.from.roi (CT, struct = S, roi.sname = "left optical",
#'                       alias = "left_optical_nerve", verbose = FALSE,
#'                       modality = "binary")
#' binr <- bin.from.roi (CT, struct = S, roi.sname = "right optical",
#'                       alias = "right_optical_nerve", verbose = FALSE,
#'                       modality = "weight")
#' 
#' view.coord <- S$roi.info[S$roi.info$roi.pseudo == "leftopticalnerve",]$Gz
#' palette <- grey.colors(100, start = 0, end = 1, 
#'                        alpha = c(0, rep(1,99)))
#' plot(S, view.coord = view.coord, main ="optical nerves")
#' plot(binl, view.coord = view.coord, col = palette, 
#'     cut.interpolate = FALSE, add = TRUE)
#' plot(binr, view.coord = view.coord, col =palette, 
#'      cut.interpolate = FALSE, add = TRUE)
#' plot(S, view.coord = view.coord, lwd = 2, add= TRUE)
#' 
#' \dontrun{
#' # with a smaller step
#' step <- 1
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name =  c("eye", "optical nerve", "brain"), 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'#' view.coord <- S$roi.info[S$roi.info$roi.pseudo == "leftopticalnerve",]$Gz
#'
#' # "optical nerve" binary without inclusions management
#' bin <- bin.from.roi (CT, struct = S, roi.sname = "left optical", 
#'                      alias = "left_optical_nerve", within = FALSE,
#'                      verbose = FALSE)
#' display.plane (CT, top = bin, struct = S, view.coord = view.coord, 
#'                legend.shift = -80, interpolate = FALSE, main = "Left nerve")
#'
#' # "optical nerve" binary with inclusions management
#' bin <- bin.from.roi (CT, struct = S, roi.sname = "left optical", 
#'                      alias = "left_optical_nerve", within = TRUE,
#'                      verbose = FALSE)
#' display.plane (CT, top = bin, struct = S, view.coord = view.coord, 
#'                legend.shift = -80, interpolate = FALSE, main = "Left nerve") 
#' }

#' @importFrom methods is
#' @export
bin.from.roi <- function (vol, struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                           T.MAT = NULL,  within = TRUE, alias = "", description = NULL,
                           modality=c("binary", "weight"), ...){
  
  if (is.null(vol)) return (NULL)
  args <- tryCatch(list(...), error = function(e)list())
  verbose <- args[["verbose"]]
  eps <- args[["eps"]]
  
  if (is.null(verbose)) verbose <- TRUE
  if (is.null(eps)) eps <- 1e-9
  
  
  modality <- modality[1]
  if (!is (struct, "struct")) stop ("struct should be a struct class object.")
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (length (roi.idx) > 1) stop ("multiple rois forbidden.")
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if(is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  if (length (vol$k.idx)>1) {if (!all (diff(vol$k.idx)==1)) stop ("planes must be contiguous.")}
  if(is.null(struct$roi.data)) stop ("empty roi.data.")
  rigid.M <- get.rigid.M (T.MAT, src.ref=struct$ref.pseudo, dest.ref = vol$ref.pseudo)
  if (is.null(rigid.M)) stop("vol and struct do not have the same frame of reference: define a valid T.MAT")
  
  warning.f <- FALSE
  if (length (roi.idx) == 0) {
    warning ("no roi selected.")
    warning.f <- TRUE
  }
  l.data <- struct$roi.data[[roi.idx]]
  if (length(l.data)==0) {
    warning ("no cutting plan in roi.")
    warning.f <- TRUE
  }
  
  type <- castlow.str(sapply(l.data,function(l) l$type))
  if (!all(type =="closedplanar" | type =="point")){
    warning ("roi is not closed_planar or point")
    warning.f <- TRUE
  }
  l.data <- l.data[type!="point"]
  if (length(l.data)==0) warning.f <- TRUE
  # if (verbose) cat(struct$roi.info$name[roi.idx]," ")
  if (is.null(description)) description <- struct$roi.info$roi.pseudo[roi.idx]
  vol.out <- vol.copy (vol, alias = alias, modality=modality, description=description, number=roi.idx)
  vol.out$vol3D.data[] <- 0
  vol.out$min.pixel <- 0
  vol.out$max.pixel <- 0
  
  if (warning.f) {
    if (alias=="") return (vol.out)
    return(.set.ref.obj(vol.out,list(struct)))
  }
  
  #2D
  idx.c <- which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),2,sum)==0)
  idx.r <-  which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    if (abs(vol$xyz0[1,idx.r])>1e-6) return(NULL)
    u <- vol$xyz.from.ijk
    u[idx.r,idx.c]<- 1
    ijk.from.contour <- solve(u)
    ijk.from.contour[idx.r,idx.c] <- 0
  } else {#3D
    ijk.from.contour <- solve(vol$xyz.from.ijk)}
  
  
  ijk.from.contour <- ijk.from.contour %*% rigid.M %*% struct$ref.from.contour
  contour.from.ijk <- solve(struct$ref.from.contour) %*%
    get.rigid.M (T.MAT, src.ref=vol$ref.pseudo, dest.ref = struct$ref.pseudo) %*% vol$xyz.from.ijk
  
  # points extremes struct
  levels <- as.numeric (sapply (l.data,function(d)d$level))
  level.idx <- 0
  level.range <- sort(unique(levels))
  if (!within) level.range <- 0
  pt.l <-lapply(level.range, function(lev){
    L <- lapply(l.data[levels==lev], function(l) {
      # pt <- polyg.sort(l$pt)
      pt <- matrix(.polygcleanC(as.vector(t(as.matrix(l$pt))), sort = TRUE), byrow = TRUE, ncol=3)
      list(as.matrix(pt[-nrow(pt),]),cbind(diff(pt[,1]),diff(pt[,2]),diff(pt[,3])))})
    list(O = do.call(rbind,lapply(L, function(l) l[[1]])),
         u = do.call(rbind,lapply(L, function(l) l[[2]])))
  })
  
  
  # on crée un volume avec un pas au plus proche de vol, dans le référentiel des plans de coupes du struct
  p <- abs(contour.from.ijk %*% cbind(c(1,0,0,0),c(0,1,0,0), c(0,0,1,0)))
  Vb.dxyz <- c(max(p[1,]), max(p[2,]),struct$thickness)
  Vb.pt000 <- (apply(contour.from.ijk %*% vol$cube.idx,1,range))[,1:3]
  Vb.pt000[2,3] <- ceiling((Vb.pt000[2,3]-pt.l[[1]]$O[1,3])/struct$thickness)*struct$thickness + pt.l[[1]]$O[1,3]
  Vb.pt000[1,3] <- floor((Vb.pt000[1,3]-pt.l[[1]]$O[1,3])/struct$thickness)*struct$thickness + pt.l[[1]]$O[1,3]
  Vb.n.ijk <- ceiling((Vb.pt000[2,]- Vb.pt000[1,])/Vb.dxyz + 1)
  Vb <- vol.create(n.ijk = Vb.n.ijk, dxyz = Vb.dxyz, modality ="weight",pt000 = Vb.pt000[1,],
                   default.value = 0,ref.pseudo = struct$ref.pseudo)
  
  
  
  # on regarde à quelle indice de Vb correpond les points du struct
  Mat <- solve (Vb$xyz.from.ijk)
  O_ijk <- lapply( 1:length(level.range),function(level.idx) cbind(pt.l[[level.idx]]$O,1) %*% t(Mat))
  range.s <- t(apply(O_ijk[[1]][,1:3],2,range))
  range.s <- cbind(floor(range.s[,1]),ceiling(range.s[,2]))
  vol3D.range  <- cbind( c(range.s[1,1]-1,max(range.s[2,1]-1,-1),max(range.s[3,1],0)),  
                         apply(cbind(range.s[,2], Vb$n.ijk - 1),1,min))
  common.range <- vol3D.range
  common.range[,1] <-apply(cbind(range.s[,1], c(0,0,0)),1,max)
  common.range[,2] <- apply(cbind(range.s[,2], Vb$n.ijk - 1),1,min)
  # vol3D.range[,2] <- common.range[,2]
  if (any(common.range[,1]>common.range[,2])) {
    if (alias=="") return (vol.out)
    return(.set.ref.obj(vol.out,list(struct)))
  }
  vol3D.nijk <-vol3D.range[,2]-vol3D.range[,1] +1
  # Vb représente ici le vol3D du référentiel des plans de coupes du struct, incluant tout le vol3D vol.
  
  
  vol3D <- array(0,dim = vol3D.nijk)
  level.idx <- 1
  
  # A <- array(1:prod(vol3D.nijk),dim = vol3D.nijk)
  # do.call(rbind, lapply(1:vol3D.nijk[3], function(k) A[vol3D.nijk[1],,k]-1))
  for (level.idx in 1:length(level.range)){
    sens <- 1-2 * (level.range[level.idx]%%2)
    u <- (cbind(pt.l[[level.idx]]$u,0)  %*% t(Mat))[,1:3]
    
    dist.u<- sqrt(u[,1]^2 +u[,2]^2 + u[,3]^2)
    u[ ,1] <- u[ ,1]/dist.u; 
    u[ ,2] <- u[ ,2]/dist.u;
    
    pO <- as.numeric(t(sweep(O_ijk[[level.idx]][,1:3] ,2,vol3D.range[,1])))
    
    vol3D  <- vol3D + sens * array(.ouline2voxC(as.numeric(t( u)),vol3D.nijk , 0:(vol3D.nijk[3]-1),  0:(vol3D.nijk[3]-1),pO,
                                                          lambda_max = dist.u, ntab = sum(ceiling(dist.u/min(Vb.dxyz))+5)),dim=vol3D.nijk)
    
  }
  
  rgi <- (common.range[1,1]:common.range[1,2])
  rgj <- (common.range[2,1]:common.range[2,2])
  rgk <- (common.range[3,1]:common.range[3,2]) 
  Vb$vol3D.data[ rgi + 1, rgj + 1, rgk + 1] <-
    vol3D[ rgi-vol3D.range[1,1]+1, rgj-vol3D.range[2,1]+1, rgk-vol3D.range[3,1]+1]
  
  # Vb$max.pixel <- max(Vb$vol3D.data, na.rm=TRUE)
  # Vb$min.pixel <- min(Vb$vol3D.data, na.rm=TRUE)
  # display.3D.stack(Vb,Vb$k.idx, border =F)
  # display.3D.contour(struct, roi.idx =roi.idx)
  # bg3d("black")
  
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
    rg.cube <- ijk.from.contour %*% Vb$xyz.from.ijk %*% cube.idx
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
    Mref <- t(solve(Vb$xyz.from.ijk) %*% contour.from.ijk)
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
  
  # display.3D.stack(vol.out,Vb$k.idx, border =F)
  # display.3D.contour(struct, roi.idx =roi.idx)
  # bg3d("black")
  # get.volume.from.bin(vol.out)
  # get.volume.from.roi(struct, roi.idx =roi.idx)
  
  if (alias=="") return (vol.out)
  return(.set.ref.obj(vol.out,list(struct)))
}

######################################################################################
#' @import progress
.bin.from.roi <- function (vol, struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                          T.MAT = NULL,  within = TRUE, alias = "", description = NULL,...){
  if (is.null(vol)) return (NULL)
  args <- tryCatch(list(...), error = function(e)list())
  verbose <- args[["verbose"]]
  eps <- args[["eps"]]
  if (is.null(verbose)) verbose <- TRUE
  if (is.null(eps)) eps <- 1e-9
  
  if (!is (struct, "struct")) stop ("struct should be a struct class object.")
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (length (roi.idx) > 1) stop ("multiple rois forbidden.")
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if(is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  if (length (vol$k.idx)>1) {if (!all (diff(vol$k.idx)==1)) stop ("planes must be contiguous.")}
  if(is.null(struct$roi.data)) stop ("empty roi.data.")
  transfert.M <- get.rigid.M (T.MAT, src.ref=vol$ref.pseudo, dest.ref = struct$ref.pseudo)
  if (is.null (transfert.M)){
    if (vol$ref.pseudo!=struct$ref.pseudo) {stop ("different ref.pseudo. Enter T.MAT")
    } else {transfert.M  <- diag (4)}
  }
  
  warning.f <- FALSE
  if (length (roi.idx) == 0) {
    warning ("no roi selected.")
    warning.f <- TRUE
  }
  if (length(struct$roi.data[[roi.idx]])==0) {
    warning ("no cutting plan in roi.")
    warning.f <- TRUE
  }
  
  type <- castlow.str(sapply(struct$roi.data[[roi.idx]],function(l) l$type))
  if (!all(type =="closedplanar" | type =="point")){
    warning ("roi is not closed_planar or point")
    warning.f <- TRUE
  }
  
  if (is.null(description)) description <- struct$roi.info$roi.pseudo[roi.idx]
  Vb <- vol.copy (vol, alias = alias, modality="binary", description=description, number=roi.idx)
  
  Vb$min.pixel <- FALSE
  Vb$max.pixel <- FALSE
  Vb$vol3D.data <- array (FALSE, Vb$n.ijk)
  
  if (warning.f) {
    if (alias=="") return (Vb)
    return(.set.ref.obj(Vb,list(struct)))
  }
  
  dz <- round(struct$thickness,6)
  cont.z <- round(sapply (struct$roi.data[[roi.idx]], function(co) co$pt[1,3]),3)
  # z.corner <- cont.z_[1]-dz/2
  # cont.z <- cont.z_
  # if (dz!=0) cont.z <- floor ((cont.z_-z.corner)/dz  )*dz + cont.z_[1]
  
  l.data <- struct$roi.data[[roi.idx]]
  levels <- as.numeric (sapply (l.data,function(d)d$level))
  M.ijk.from.contour <- solve (transfert.M %*% Vb$xyz.from.ijk ) %*% struct$ref.from.contour
  M.contour.from.ijk <- solve(M.ijk.from.contour)
  
  if (verbose) 
    pb <- progress_bar$new(format = paste(struct$roi.info$roi.pseudo[roi.idx],"[:bar] :percent"),
                           total = length(l.data), width= 60)
  
  for (level.idx in sort(unique(levels))){
    
    for (j in which (levels==level.idx)) {
      
      z <- cont.z[j] #struct$roi.data[[roi.idx]][[j]]$pt$z[1]
      
      #par plan contour, on cherche les indices du volume qui pourraient être concernés
      roi.x <- c(min(l.data[[j]]$pt[,1]), max(l.data[[j]]$pt[,1]))
      roi.y <- c(min(l.data[[j]]$pt[,2]), max(l.data[[j]]$pt[,2]))
      roi.cube <- matrix(c(roi.x[1],roi.y[1], z-dz/2, 1, roi.x[2],roi.y[1],z-dz/2, 1,
                           roi.x[2],roi.y[2], z-dz/2, 1, roi.x[1],roi.y[2],z-dz/2, 1,
                           roi.x[1],roi.y[1], z+dz/2, 1, roi.x[2],roi.y[1],z+dz/2, 1,
                           roi.x[2],roi.y[2], z+dz/2, 1, roi.x[1],roi.y[2],z+dz/2, 1),
                         ncol=8, byrow = FALSE)
      roi.idx.in.map <- M.ijk.from.contour %*% roi.cube
      rg.i <- floor(min(roi.idx.in.map[1,])):ceiling(max(roi.idx.in.map[1,]))
      rg.j <- floor(min(roi.idx.in.map[2,])):ceiling(max(roi.idx.in.map[2,]))
      rg.k <- floor(min(roi.idx.in.map[3,])):floor(max(roi.idx.in.map[3,]))
      
      idx.matrix <- as.matrix(expand.grid(i= rg.i[!(is.na (match(rg.i, 0:(Vb$n.ijk[1]-1))))],
                                          j= rg.j[!(is.na (match(rg.j, 0:(Vb$n.ijk[2]-1))))],
                                          k= rg.k[!(is.na (match(rg.k, Vb$k.idx)))], t=1) )
      Ridx<- idx.matrix[,1:3] + 1
      
      
      coords <- round(idx.matrix %*%  t(M.contour.from.ijk),6)
      if (dz!=0) {
        # plan.z.flag <- ceiling(round(((coords[,3]-z-dz/2)/dz ) +0.5,3)) == 0
        plan.z.flag <- floor((coords[,3]-z)/dz  +  0.5) == 0
        # plan.z.flag <- (coords[,3]>=z-dz) & (coords[,3]<z+dz)
      } else {
        plan.z.flag <- coords[,3]==z
      }
      keep <- .pt.in.polygon (coords[plan.z.flag, 1], coords[plan.z.flag, 2],
                              l.data[[j]]$pt[,1],
                              l.data[[j]]$pt[,2]) > 0.5
      # plot(l.data[[j]]$pt$x, l.data[[j]]$pt$y, type="l", lwd=2, xlim=range(coords[plan.z.flag,1]), ylim=range(coords[plan.z.flag,2]))
      # points (coords[plan.z.flag,1],coords[plan.z.flag,2],col=keep+1, pch=16)
      
      if (within){ value <- (level.idx %% 2) == 0
      }else {value <- TRUE}
      Vb$vol3D.data [Ridx[plan.z.flag, ]][keep] <- value #Vb$vol3D.data [Ridx[plan.z.flag, ]] | keep
      if (verbose) pb$tick()
    }
  }
  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  
  if (alias=="") return (Vb)
  return(.set.ref.obj(Vb,list(struct)))
}


##############################################################################################
