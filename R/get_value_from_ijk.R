#' Value of the volume at a selection of DICOM indices
#' @description The \code{get.value.from.ijk} function calculates the value of 
#' a "volume" class object at DICOM indices i, j, k, whether they are 
#' integers or not.
#' @param ijk Vector or 3-column matrix of DICOM indices.
#' @param vol "volume" class object.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @param s.ijk Vector of 3 positive numbers greater than or equal to 1,
#' representing the new voxel ijk-sizes in which averaging is calculated.
#' @return Returns a vector of the values of the volume at the requested DICOM 
#' indices.
#' @seealso \link[espadon]{get.ijk.from.index}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#' # isodose
#' Dmax <- max (D$vol3D.data, na.rm = TRUE)
#' Dmax
#' idx <- which (D$vol3D.data >= (Dmax -1) & D$vol3D.data <= (Dmax - 0.2))
#' ijk <- get.ijk.from.index (idx, D)
#' get.value.from.ijk (ijk, vol = D, interpolate = FALSE)

#' @export
#' @importFrom methods is
get.value.from.ijk <- function (ijk, vol, interpolate = TRUE, s.ijk = c(1,1,1))  {
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  s_ijk <- abs(s.ijk)
  # s_ijk [s_ijk < 1] <- 1
  if (length(s_ijk) < 3) s_ijk <- c(s_ijk, rep(s_ijk[length(s_ijk) ],2))[1:3]
  # if (length(s.ijk)!=3 | !all(s_ijk[1:min(length(s.ijk),3)]==s.ijk[1:min(length(s.ijk),3)])) 
  #   warning("s.ijk is set to c(", paste(s_ijk, collapse = ", "),")")
  toC3M <- function (vect) {return(matrix(vect,ncol=3))}
  ijk <- toC3M (ijk)

  k.idx<- match(0:max(vol$k.idx),vol$k.idx)
  k.loc <- k.idx-1
  fna <- is.na(k.idx)
  k.idx[!fna] <- vol$k.idx
  k.loc[fna] <- max(vol$k.idx)+1
  k.idx[fna] <- max(vol$k.idx)+1
  vol3D <- as.numeric(vol$vol3D.data)
  # vol3D[is.na(vol3D)] <- NaN
  value <- .getvaluefromijkC (vol3D = vol3D,
                              interpolate = interpolate,
                              i = as.numeric(ijk[ ,1]),
                              j = as.numeric(ijk[ ,2]),
                              k = as.numeric(ijk[ ,3]),
                              k_idx = k.idx,
                              k_loc = k.loc, n_ijk=vol$n.ijk,
                              s_ijk = s_ijk)
  
  # value [is.nan(value)] <- NA
  return(value)
}