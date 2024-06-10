####################################################################
#' Creation of a binary volume according to RoI
#' @description The \code{bin.from.roi} function creates a "volume" class 
#' object, of "binary" modality, in which all the voxels of a RoI are set to \code{TRUE}.
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
#' @param ... additional argument such as verbose
#' @return Returns a "volume" class object of "binary" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol}, in which the voxels in the RoI are set to TRUE.
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
#' # "optical nerve" binary without inclusions management
#' bin <- bin.from.roi (CT, struct = S, roi.sname = "left optical", 
#'                      alias = "left_optical_nerve", verbose = FALSE)
#' display.plane (CT, top = bin, struct = S,
#'                view.coord = S$roi.info[S$roi.info$roi.pseudo == "leftopticalnerve",]$Gz, 
#'                legend.shift = -80, interpolate = FALSE, main = "Left nerve selection")
#' 
#' \dontrun{
#' # with a smaller step
#' step <- 1
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name =  c("eye", "optical nerve", "brain"), 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' # "optical nerve" binary without inclusions management
#' bin <- bin.from.roi (CT, struct = S, roi.sname = "left optical", 
#'                      alias = "left_optical_nerve", within = FALSE,
#'                      verbose = FALSE)
#' display.plane (CT, top = bin, struct = S,
#'                view.coord = S$roi.info[S$roi.info$roi.pseudo == "leftopticalnerve",]$Gz, 
#'                legend.shift = -80, interpolate = FALSE, main = "Left nerve selection")
#'
#' # "optical nerve" binary with inclusions management
#' bin <- bin.from.roi (CT, struct = S, roi.sname = "left optical", 
#'                      alias = "left_optical_nerve", within = TRUE,
#'                      verbose = FALSE)
#' display.plane (CT, top = bin, struct = S,
#'                view.coord = S$roi.info[S$roi.info$roi.pseudo == "leftopticalnerve",]$Gz, 
#'                legend.shift = -80, interpolate = FALSE, main = "Left nerve selection") 
#' }

#' @importFrom methods is
#' @import progress


  
#' @export
bin.from.roi <- function (vol, struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
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
.bin.from.roi <- function (vol, struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                           T.MAT = NULL,  within = TRUE, alias = "", description = NULL,...){
  args <- tryCatch(list(...), error = function(e)list())
  verbose <- args[["verbose"]]
  eps <- args[["eps"]]
  if (is.null(verbose)) verbose <- TRUE
  if (is.null(eps)) eps <- 1e-9
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  
  if (length (roi.idx) == 0) {
    warning ("no roi selected.")
    return (NULL)
  }
  if (length (roi.idx) > 1) {
    warning ("multiple rois forbidden.")
    return (NULL)
  }
  
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if (!is (struct, "struct")) {
    warning ("struct should be a struct class object.")
    return (NULL)
  }
  
  if (length(struct$roi.data[[roi.idx]])==0) {
    warning ("no cut in roi.")
    return (NULL)
  }
  
  type <- castlow.str(sapply(struct$roi.data[[roi.idx]],function(l) l$type))
  if (!all(type =="closedplanar" | type =="point")){
    warning ("roi is not closed_planar or point")
    return (NULL)
  }
  
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  if (length (vol$k.idx)>1){
    if (!all (diff(vol$k.idx)==1)) {
      warning ("planes must be contiguous.")
      return (NULL)
    }
  }
  
  if(is.null(struct$roi.data)){
    warning ("empty roi.data.")
    return (NULL)
  }
  
  
  transfert.M <- get.rigid.M (T.MAT, src.ref=vol$ref.pseudo, dest.ref = struct$ref.pseudo)
  
  if (is.null (transfert.M)){
    if (vol$ref.pseudo!=struct$ref.pseudo) {
      warning ("different ref.pseudo. Enter T.MAT")
      return (NULL)
    } else {
      transfert.M  <- diag (4)
    }
  }
  if (is.null(description)) description <- struct$roi.info$roi.pseudo[roi.idx]
  Vb <- vol.copy (vol, alias = alias, modality="binary", description=description, number=roi.idx)
  
  Vb$min.pixel <- 0
  Vb$max.pixel <- 1
  Vb$vol3D.data <- array (FALSE, Vb$n.ijk)
  
  dz <- round(struct$thickness,6)
  thickness <- struct$thickness
  if (thickness ==0)  thickness <- Vb$slice.thickness
  l.data <- roidata.keyhole.to.xor(struct$roi.data[[roi.idx]], thickness = thickness)
  
  cont.z <- round(sapply (l.data, function(co) co$pt[1,3]),3)
  
  # z.corner <- cont.z_[1]-dz/2
  # cont.z <- cont.z_
  # if (dz!=0) cont.z <- floor ((cont.z_-z.corner)/dz  )*dz + cont.z_[1]
  
  levels <- as.numeric (sapply (l.data,function(d)d$level))
  M.ijk.from.contour <- solve (transfert.M %*% Vb$xyz.from.ijk ) %*% struct$ref.from.contour
  M.contour.from.ijk <- solve(M.ijk.from.contour)
  if (verbose) 
    pb <- progress_bar$new(format = paste(struct$roi.info$roi.pseudo[roi.idx],"[:bar] :percent"),
                           total = length(l.data), width= 60)
  for (level.idx in sort(unique(levels))){
    for (j in which (levels==level.idx)) {
      z <- round(cont.z[j],6)#struct$roi.data[[roi.idx]][[j]]$pt$z[1]
      #par plan contour, on cherche les indices du volume qui pourraient être concernés
      roi.L <- list(x = c(min(l.data[[j]]$pt[,1]), max(l.data[[j]]$pt[,1])),
                    y = c(min(l.data[[j]]$pt[,2]), max(l.data[[j]]$pt[,2])))
      
      roi.cube <- matrix(c(roi.L$x[1],roi.L$y[1], z-dz/2, 1, roi.L$x[2],roi.L$y[1],z-dz/2, 1,
                           roi.L$x[2],roi.L$y[2], z-dz/2, 1, roi.L$x[1],roi.L$y[2],z-dz/2, 1,
                           roi.L$x[1],roi.L$y[1], z+dz/2, 1, roi.L$x[2],roi.L$y[1],z+dz/2, 1,
                           roi.L$x[2],roi.L$y[2], z+dz/2, 1, roi.L$x[1],roi.L$y[2],z+dz/2, 1),
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
      
      if (any(plan.z.flag)){
        keep <- rep(FALSE, nrow(coords[plan.z.flag, ,drop=FALSE]))
        coord.L <- list(sort(unique(coords[plan.z.flag, 1])), sort(unique(coords[plan.z.flag, 2])))
        rg.L <- lapply(1:2, function(i) floor(range(c(coord.L[[i]],roi.L[[i]]))+c(0,0.5)) + c(-1,+1))
        # ,roi.L[[i]]
        le <- sapply(coord.L, length)
        x.sup.y <- as.numeric(length(coord.L[[1]])>length(coord.L[[2]]))
        
        M <- matrix(.roiinterC(pt1_x = ((rep(rg.L[[2 - x.sup.y]][1],le[x.sup.y+1]) * x.sup.y + (1-x.sup.y)*coord.L[[1 + x.sup.y]]))[1:le[x.sup.y+1]], 
                               pt1_y = (((rep(rg.L[[2 - x.sup.y]][1],le[x.sup.y+1]) * (1-x.sup.y)) + x.sup.y*coord.L[[1 + x.sup.y]]))[1:le[x.sup.y+1]],
                               u1_x = x.sup.y, u1_y = 1-x.sup.y, d1=rg.L[[2 - x.sup.y]][2]-rg.L[[2 - x.sup.y]][1],
                               pt2_x = l.data[[j]]$pt[,1], pt2_y = l.data[[j]]$p[,2], eps =eps),
                    nrow = le[x.sup.y+1])
        if (length(M)>1){
          m.coord <- match(coords[plan.z.flag, x.sup.y+1],coord.L[[1 + x.sup.y]])
          for (col.idx in 1:(ncol(M)/2)){
            keep  <- keep | ((M[m.coord, col.idx*2 -1]<=coords[plan.z.flag, 2 - x.sup.y]) & (coords[plan.z.flag, 2 - x.sup.y] <= M[m.coord, col.idx*2]))
          }
        }
        
        if (within){ value <- (level.idx %% 2) == 0
        }else {value <- TRUE}
        Vb$vol3D.data [Ridx[plan.z.flag, ]][keep] <- value #Vb$vol3D.data [Ridx[plan.z.flag, ]] | keep
      }
      if (verbose) pb$tick()
    }
  }
  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  
  if (alias=="") return (Vb)
  return(.set.ref.obj(Vb,list(struct)))
}
