#' Conversion of the indices of a point, into xyz coordinate vector in the patient's 
#' frame of reference
#' @description The \code{get.xyz.from.index} function converts the indices of a 
#' voxel of \code{vol$vol3D.data} (for example, obtained with the function 
#' \code{which}) into a vector or matrix of x, y, z coordinates in the patient's 
#' frame of reference.
#' @param idx Index, or matrix of voxel indices in the array \code{vol$vol3D.data}. 
#' The first index of the array is 1.
#' @param vol "volume" class object.
#' @return Returns a column-matrix of coordinates in the patient's reference frame, 
#' corresponding to the indices \code{idx}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for better
#' # result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#'
#' # voxels location where the dose is greater than 99.9% Dmax
#' Dmax <- max (D$vol3D.data, na.rm = TRUE) # D$max.pixel
#' get.xyz.from.index (which (D$vol3D.data >= 0.99 * Dmax), D)
#' # or
#' get.xyz.from.index (which (D$vol3D.data >= 0.99 * Dmax, arr.ind = TRUE), D)

#' @export
get.xyz.from.index <- function (idx, vol){
  if (length(idx)==0) return(NULL)
  idx_ <- .get.ijkt.from.index (idx, vol)
  pt <- matrix ((idx_ %*%  t(vol$xyz.from.ijk))[,1:3], ncol=3)
  colnames (pt) <- c("x","y","z")
  return (pt)
}