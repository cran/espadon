#' Conversion of the indices of a point into ijk vector
#' @description The \code{get.ijk.from.index} function converts the voxel indices 
#' of \code{vol$vol3D.data} (for example, obtained with the function \code{which}) 
#' into a vector or matrix of DICOM indices i, j, k.
#' @param idx Index, or matrix of voxel indices of the array \code{vol$vol3D.data}.
#' @param vol "volume" class object.
#' @return Returns an i, j, k column matrix of the DICOM indices of the points 
#' of \code{vol$vol3D.data}.
#' @seealso \link[espadon]{get.value.from.ijk}, \link[espadon]{display.kplane}
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#'
#' # voxels location where the dose is greater than 99.9% Dmax
#' Dmax <- max (D$vol3D.data, na.rm = TRUE) # D$max.pixel
#' get.ijk.from.index (which (D$vol3D.data >= 0.999 * Dmax), D)
#' # or
#' get.ijk.from.index (which (D$vol3D.data >= 0.999 * Dmax, arr.ind = TRUE), D)
#'
#' ijk <- as.numeric (get.ijk.from.index (which.max (D$vol3D.data), D))
#' display.kplane (D, k = ijk[3]) 

#' @export
get.ijk.from.index <- function (idx, vol){
  
  return(.get.ijkt.from.index (idx, vol)[,1:3])
  
}
