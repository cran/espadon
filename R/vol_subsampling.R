#' Subsampling a volume
#' @description The \code{vol.subsampling} function sub-samples the grid of a "volume" class 
#' object.
#' @param vol "volume" class object.
#' @param fact.ijk Strictly positive integer, or a vector of 3 strictly positive integers.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be \code{paste ("subsampling" ,vol$description)}.
#' @return Returns a "volume" class object, in which 3D volume grid is 
#' subsampled: the voxel size is multiplied by \code{fact.ijk} and the center location 
#' of the volume is invariant.
#' @seealso \link[espadon]{vol.oversampling}.
#' @examples
#' vol <- vol.create(n.ijk = c(10,10,1),dxyz = c(2,2,2), ref.pseudo = "ref1", 
#'                   modality ="test", pt000 = c(0,0,0))
#' vol$vol3D.data[] <- array(1:prod(vol$n.ijk), dim = vol$n.ijk)
#' vol$max.pixel <- prod(vol$n.ijk)
#' vol$min.pixel <- 1
#' mid <- as.numeric (apply (get.extreme.pt (vol), 1, mean))
#
#' vol_ss <- vol.subsampling (vol, fact.ijk= 2)
#' mid_ss <- as.numeric (apply (get.extreme.pt (vol_ss), 1, mean))
#
#' display.plane(vol,interpolate = FALSE, view.coord = mid[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid[1], mid[2], pch=16, col="red")
#' display.plane(vol_ss,interpolate = FALSE, view.coord = mid_ss[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid_ss[1], mid_ss[2], pch=16, col="red")
#' @export  

vol.subsampling <- function(vol, fact.ijk= 2, interpolate = TRUE, alias = "", description = NULL){
  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  if (!((length(fact.ijk)==1 | length(fact.ijk)==3)) | 
      !all(abs(as.integer(fact.ijk)) == fact.ijk) | 
      any(fact.ijk < 0)){
    stop ("fact.ijk must be a strictly positive integer or a vector of 3 strictly positive integers")
  }
  if (length(fact.ijk)==1) fact.ijk <- rep(fact.ijk,3)
  t.mat <- ref.cutplane.add(vol, ref.cutplane = "rcp" )
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo="rcp", T.MAT =  t.mat)
  new.dxyz <- vol_$dxyz
  new.n.ijk <- vol_$n.ijk
  f <- new.n.ijk<2
  fact <- fact.ijk
  fact[f] <- 1
  new.dxyz <-  new.dxyz*fact
  true.n.ijk <- c(vol$n.ijk[1:2], vol$k.idx[vol$n.ijk[3]]-vol$k.idx[1]+1)
  new.n.ijk <- ceiling(true.n.ijk/fact)
  
  center.idx <-c(true.n.ijk/2 - 0.5, 0)
  
  if (is.null(description))  description <- paste("subsampling", vol$description)
  back.vol <- vol.create(n.ijk = new.n.ijk, dxyz=new.dxyz,
                         mid.pt = as.numeric(center.idx %*% t(vol_$xyz.from.ijk))[1:3], 
                         ref.pseudo = "rcp", modality = vol_$modality, alias = alias, description = description)
  ijkt <- as.matrix(expand.grid(0:(back.vol$n.ijk[1]-1),0:(back.vol$n.ijk[2]-1),0:(back.vol$n.ijk[3]-1),0))
  xyz <-  ijkt %*% t(back.vol$xyz.from.ijk)
  
  back.vol$vol3D.data[ijkt[,1:3]+1] <- get.value.from.xyz(xyz[, 1:3],vol_,interpolate =  interpolate)
  if (any(!is.na(back.vol$vol3D.data))) {
    back.vol$min.pixel <- min(back.vol$vol3D.data, na.rm = TRUE)
    back.vol$max.pixel <- max(back.vol$vol3D.data, na.rm = TRUE)
  } else {
    back.vol$min.pixel <- NA
    back.vol$max.pixel <- NA
  }
  vol.in.new.ref (back.vol,new.ref.pseudo=vol$ref.pseudo, T.MAT =  t.mat,
                  alias = alias,
                  description = description)   
}