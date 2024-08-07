#' Transform the grid of a volume class object into the grid of another
#' @description The \code{vol.regrid} function transforms the grid of a volume 
#' according to the grid of another.
#' @param vol "volume" class object to regrid.
#' @param back.vol "volume" class object whose grid will be used for regriding. 
#' Its \code{$ref.pseudo} must exist in the \code{T.MAT} list.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}. If \code{T.MAT = NULL}, 
#' \code{back.vol$ref.pseudo} must be equal to \code{vol$ref.pseudo}.
#' @param interpolate Boolean, default to \code{TRUE}. 
#' @param method method of interpolation, defaulted to 'NN' ie 'Nearest Neighbor'. See Details. 
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be that of \code{vol}.
#' @param verbose Boolean. if \code{TRUE} (default) a progress bar is displayed.
#' @details The interpolation method is chosen from:
#' \itemize{
#' \item \code{'NN'}: the value of a voxel is calculated from its nearest adajcent neighbors.
#' \item \code{'Av'}: the value of a voxel is the weighted average of 
#' the voxels contained in a box, whose sides are automatically calculated from 
#' the \code{back.bol$dxyz} steps. 
#' }
#' @return Returns a copy of \code{vol}, in which grid is that of \code{back.vol}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c ("mr", "rtdose"), 
#'                              dxyz = rep (step, 3), beam.nb = 4)
#' MR <- patient$mr[[1]]
#' D <- patient$rtdose[[1]]
#'
#' # change grid
#' D.on.MR <- vol.regrid (vol = D, back.vol = MR, interpolate = TRUE,
#'                        T.MAT = patient$T.MAT, alias = "",
#'                        description = NULL, verbose = FALSE)
#'
#' # maximum dose location
#' max.dose.in.MR <- get.xyz.from.index (which.max (D.on.MR$vol3D.data), D.on.MR)
#' display.plane (bottom = MR, view.coord = max.dose.in.MR[3],
#'                top= D.on.MR, bottom.col = grey.colors(255, start = 0, end = 1),
#'                bottom.breaks = seq (0, 500, length.out = 256),
#'                bg = "#00ff00", interpolate = FALSE)

#' @import progress
#' @importFrom methods is
#' @export
#' 
vol.regrid <- function (vol, back.vol, T.MAT = NULL, interpolate = TRUE, 
                        method = c("NN","Av"), alias = "", 
                        description=NULL, verbose = TRUE){
  
  if (is.null(vol)) return (NULL)
  if (!is (vol, "volume")) stop ( "vol should be a volume class object.")
  if (!is (back.vol, "volume")) stop ("back.vol should be a volume class object.")
  if (is.null (vol$vol3D.data)) stop("vol should contain vol3D.data. Load voxel data.")
  if (is.null (back.vol$vol3D.data)) stop("back.vol should contain vol3D.data. Load voxel data.")
  if (!is.logical(interpolate)) stop ("interpolate should be TRUE or FALSE.")
  if (!is.null(T.MAT)){
    if (!is (T.MAT, "t.mat")) stop ("T.MAT should be a t.mat class object or NULL.")
  }

  method <- method[1]
  send <- FALSE
  if (is.null(description)) description<- vol$description
  
  if (grid.equal (vol, back.vol)) {
    new.vol <- vol.copy (vol, alias = alias, modality = vol$modality, description =description)
    send <- TRUE
  } else {
    M.Tref <- get.rigid.M (T.MAT, back.vol$ref.pseudo, vol$ref.pseudo)
    if (is.null (M.Tref)) 
      stop ("vol$ref.pseudo and back.vol$ref.pseudo are different. Specify T.MAT")

    
    new.vol <- vol.copy (back.vol, alias = alias, modality = vol$modality, 
                         description =description)
    new.vol$patient <- vol$patient
    new.vol$patient.name <- vol$patient.name
    new.vol$patient.bd <- vol$patient.bd
    new.vol$patient.sex <- vol$patient.sex
    new.vol$acq.date <- vol$acq.date
    new.vol$study.date <- vol$study.date
    new.vol$creation.date <- vol$creation.date
    new.vol$study.time <- vol$study.time
    new.vol$number <- vol$number
    new.vol$unit <- vol$unit
    
    new.vol$vol3D.data[] <- NA
  }
  if (alias != "")  new.vol <- .set.ref.obj(new.vol, list(vol,back.vol), add=FALSE)
  
  if (send) return (new.vol)
  if (verbose) pb <- progress_bar$new(format = " processing [:bar] :percent",
                                      total = 2, width= 60)
  
  # if (is.na(vol$min.pixel) | is.na(vol$max.pixel)) return (new.vol)
  # if (abs (det(vol$xyz.from.ijk))<1e-6) {
  #   warning ("non invertible matrix, check arguments.")
  #   return (NULL)
  # }
  # M.old.from.new <- solve (vol$xyz.from.ijk) %*% M.Tref %*% new.vol$xyz.from.ijk
  # if (abs (det(M.old.from.new))<1e-6) {
  #   warning ("non invertible matrix, check arguments.")
  #   return (NULL)
  # }
  # 
  
  M.xyz.from.new <- M.Tref %*% new.vol$xyz.from.ijk
  #2D
  idx.c <- which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
  idx.r <-  which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  
  cube.idx <- vol$cube.idx
  cube.idx[1,c(1,4,5,8)] <- cube.idx[1,c(1,4,5,8)] - 0.5
  cube.idx[1,c(2,3,6,7)] <- cube.idx[1,c(2,3,6,7)] + 0.5
  cube.idx[2,c(1,2,5,6)] <- cube.idx[2,c(1,2,5,6)] - 0.5
  cube.idx[2,c(3,4,7,8)] <- cube.idx[2,c(3,4,7,8)] + 0.5
  cube.idx[3,1:4] <- cube.idx[3,1:4] - 0.5
  cube.idx[3,5:8] <- cube.idx[3,5:8] + 0.5
  if (length(idx.c)>0) {
    if(abs(vol$xyz.from.ijk[idx.r,4]-M.xyz.from.new[idx.r,4])>1e-6) return(new.vol)
    u <- vol$xyz.from.ijk 
    u[idx.r,idx.c]<- 1
    M.old.from.new <- solve(u)  %*% M.xyz.from.new
    
    M.xyz.from.new[which(apply(abs(M.xyz.from.new[1:3,1:3]),2,sum)==0) ,
                   which(apply(abs(M.xyz.from.new[1:3,1:3]),1,sum)==0)]<- 1
    new.cube.idx <- solve(M.xyz.from.new )%*% vol$xyz.from.ijk %*% cube.idx
    new.cube.idx[abs(new.cube.idx) < 1.0e-9] <- 0
  
  } else {
    M.old.from.new <- solve (vol$xyz.from.ijk) %*% M.Tref %*% new.vol$xyz.from.ijk
    new.cube.idx <- solve (M.old.from.new) %*% cube.idx
    new.cube.idx[abs(new.cube.idx) < 1.0e-9] <- 0
  }
  
  
  # new.cube.idx <- solve (M.old.from.new) %*% vol$cube.idx
  rg.i <- floor (min(new.cube.idx[1,])):ceiling(max(new.cube.idx[1,]))
  rg.j <- floor (min(new.cube.idx[2,])):ceiling(max(new.cube.idx[2,]))
  rg.k <- floor (min(new.cube.idx[3,])):ceiling(max(new.cube.idx[3,]))
  rg.i <- rg.i[!is.na (match (rg.i,(1:new.vol$n.ijk[1])-1))]
  rg.j <- rg.j[!is.na (match (rg.j,(1:new.vol$n.ijk[2])-1))]
  back.rgk.loc <- match (rg.k, new.vol$k.idx)
  fna <- is.na(back.rgk.loc)
  rg.k <- rg.k[!fna]
  
  ijk.selection <- as.matrix (expand.grid (rg.i, rg.j, rg.k, 1))
  
  back.rgk.loc <-  back.rgk.loc[!fna]
  ijk.Rselection <-  as.matrix (expand.grid (rg.i+1, rg.j+1, back.rgk.loc))
  
  
  k.idx<- match(0:max(vol$k.idx),vol$k.idx)
  k.loc <- k.idx-1
  fna <- is.na(k.idx)
  k.idx[!fna] <- vol$k.idx
  k.loc[fna] <- max(vol$k.idx)+1
  k.idx[fna] <- max(vol$k.idx)+1
  
  ijk <-matrix((ijk.selection %*% t(M.old.from.new))[ ,1:3, drop = FALSE],ncol=3)
  ijk[abs(ijk)<1e-9] <- 0
  if (verbose) pb$tick()
  vol3D <- as.numeric(vol$vol3D.data)
  vol3D[is.na(vol3D)] <- NaN
  s_ijk <- c(1,1,1)
  if (method == "Av"){
    s_ijk <-apply(abs(rbind(c(1,0,0,0),c(0,1,0,0), c(0,0,1,0)) %*% t(M.old.from.new))[,1:3, drop =FALSE],1,max)
  }
  new.vol$vol3D.data[ijk.Rselection] <- round(.getvaluefromijkC (vol3D = vol3D,
                                                           interpolate = interpolate,
                                                           i = as.numeric(ijk[ ,1]),
                                                           j = as.numeric(ijk[ ,2]),
                                                           k = as.numeric(ijk[ ,3]),
                                                           k_idx = k.idx,
                                                           k_loc = k.loc, n_ijk=vol$n.ijk,s_ijk = s_ijk),9)
  
  if (verbose) pb$tick()
  new.vol$vol3D.data[is.nan(new.vol$vol3D.data)] <- NA
  # new.vol$vol3D.data[abs(new.vol$vol3D.data)<=1e-9] <- 0
  if (new.vol$modality == "binary") {
    new.vol$vol3D.data <- new.vol$vol3D.data>=0.5
  } 
  if (any(!is.na(new.vol$vol3D.data))){
    new.vol$min.pixel <- min ( new.vol$vol3D.data, na.rm=TRUE)
    new.vol$max.pixel <- max ( new.vol$vol3D.data, na.rm=TRUE)
  } else {
    new.vol$min.pixel <- NA
    new.vol$max.pixel <- NA
  }
  return (new.vol)

}
