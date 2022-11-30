#' Oversampling a volume
#' @description The \code{vol.oversampling} function oversamples the grid of a "volume" class 
#' object.
#' @param vol "volume" class object.
#' @param fact.ijk Strictly positive integer, or a vector of 3 strictly positive integers.
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
#' vol_os <- vol.oversampling (vol, fact.ijk= 2)
#' mid_os <- as.numeric (apply (get.extreme.pt (vol_os), 1, mean))
#' 
#' display.plane(vol,interpolate = FALSE, view.coord = mid[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid[1], mid[2], pch=16, col="red")
#' display.plane(vol_os,interpolate = FALSE, view.coord = mid_os[3], 
#'               abs.rng = c(-5,25), ord.rng = c(-5,25), bg="green")
#' points (mid_os[1], mid_os[2], pch=16, col="red")
#' @export  

vol.oversampling <- function(vol, fact.ijk = 2, alias = "", description = NULL){
  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  
  if (!((length(fact.ijk)==1 | length(fact.ijk)==3)) | 
      !all(abs(as.integer(fact.ijk)) == fact.ijk) | 
      any(fact.ijk < 0)){
    stop ("fact.ijk must be a strictly positive integer or a vector of 3 strictly positive integers")
  }
  
  if (length(fact.ijk)==1) fact.ijk <- rep (fact.ijk, 3)
  
  t.mat <- ref.cutplane.add(vol, ref.cutplane = "rcp" )
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo="rcp", T.MAT =  t.mat)
  new.dxyz <- vol_$dxyz
  true.n.ijk <- c(vol$n.ijk[1:2], vol$k.idx[vol$n.ijk[3]]-vol$k.idx[1]+1)
  new.n.ijk <- true.n.ijk
  
  f <- new.n.ijk>1
  new.dxyz[f] <-  new.dxyz[f]/fact.ijk[f]
  new.n.ijk[f] <- (new.n.ijk[f]-1)*fact.ijk[f] + 1
  back.vol <- vol.create(n.ijk = new.n.ijk, dxyz=new.dxyz, pt000 = vol_$xyz0[1,], ref.pseudo = "rcp")
  if (is.null(description))  description <- paste("oversampling", vol$description)
  vol_ <- vol.regrid(vol_,back.vol, alias = vol$object.alias, 
                     description = description)
  
  vol.in.new.ref (vol_,new.ref.pseudo=vol$ref.pseudo, T.MAT =  t.mat,
                  alias = alias,
                  description = description)                  
}