#' Oversampling a volume
#' @description The \code{vol.oversampling} function oversamples the grid of a "volume" class 
#' object.
#' @param vol "volume" class object.
#' @param fact.ijk Strictly positive integer, or a vector of 3 strictly positive integers.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be \code{paste ("oversampling" ,vol$description)}.
#' @return Returns a "volume" class object, in which 3D volume grid is 
#' oversampled: the voxel size is divided by \code{fact.ijk}.
#' @seealso \link[espadon]{vol.subsampling}.
#' @examples
#' vol <- vol.create(n.ijk = c(10,10,1),dxyz = c(2,2,2), ref.pseudo = "ref1", 
#'                   modality ="test", pt000 = c(0,0,0))
#' vol$vol3D.data[] <- array(1:prod(vol$n.ijk), dim = vol$n.ijk)
#' vol$max.pixel <- prod(vol$n.ijk)
#' vol$min.pixel <- 1
#' mid <- as.numeric (apply (get.extreme.pt (vol), 1, mean))
#'  
#' vol_os <- vol.oversampling (vol, fact.ijk= c(2,2,1))
#' mid_os <- as.numeric (apply (get.extreme.pt (vol_os), 1, mean))
#' 
#' display.plane(vol,interpolate = FALSE, view.coord = mid[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid[1], mid[2], pch=16, col="red")
#' display.plane(vol_os,interpolate = FALSE, view.coord = mid_os[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid_os[1], mid_os[2], pch=16, col="red")
#' @export  

vol.oversampling <- function(vol, fact.ijk = 2, alias = "",interpolate = TRUE, description = NULL){
  

  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
 
  
  if (!((length(fact.ijk)==1 | length(fact.ijk)==3)) | 
      !all(abs(as.integer(fact.ijk)) == fact.ijk) | 
      any(fact.ijk <= 0)){
    stop ("fact.ijk must be a strictly positive integer or a vector of 3 strictly positive integers")
  }
 
    
  if (length(fact.ijk)==1) fact.ijk <- rep (fact.ijk, 3)
  if (all(fact.ijk==c(1,1,1))) return (vol)
  
  
  t.mat <- ref.cutplane.add(vol, ref.cutplane = "rcp" )
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo="rcp", T.MAT =  t.mat)
  new.dxyz <- vol_$dxyz
  true.n.ijk <- c(vol$n.ijk[1:2], vol$k.idx[vol$n.ijk[3]]-vol$k.idx[1]+1)
  new.n.ijk <- true.n.ijk
  
  # f <- new.n.ijk>1
  f <- rep(TRUE,3)
  new.dxyz[f] <-  new.dxyz[f]/fact.ijk[f]
  new.n.ijk[f] <- new.n.ijk[f] * fact.ijk[f] #(new.n.ijk[f]-1)*fact.ijk[f] + 1
  #back.vol <- vol.create(n.ijk = new.n.ijk, dxyz=new.dxyz, pt000 = vol_$xyz0[1,], ref.pseudo = "rcp")
  back.vol <- vol.create(n.ijk = new.n.ijk, dxyz=new.dxyz, 
                         pt000 = vol_$xyz0[1,]+(new.dxyz-vol_$dxyz)/2, ref.pseudo = "rcp")
  
  if (is.null(description))  description <- paste("oversampling", vol$description)
  vol_ <- vol.regrid(vol_,back.vol, alias = "", interpolate=interpolate, description = description)
  vol_$object.alias <- vol$object.alias
  vol_$object.info <- vol$object.info
  
  if (vol$modality == "binary") vol_$vol3D.data <- vol_$vol3D.data >= 0.5
  if (any(!is.na(vol_$vol3D.data))) {
    vol_$min.pixel <- min(vol_$vol3D.data, na.rm = TRUE)
    vol_$max.pixel <- max(vol_$vol3D.data, na.rm = TRUE)
  } else {
    vol_$min.pixel <- NA
    vol_$max.pixel <- NA
  }
  
  vol.in.new.ref (vol_,new.ref.pseudo=vol$ref.pseudo, T.MAT =  t.mat,
                  alias = alias,
                  description = description)                  
}