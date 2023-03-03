#' Voxel values on a selection of points
#' @description The \code{get.value.from.xyz} function calculates the voxel values at
#' the x, y, z coordinate points in the chosen frame of reference.
#' @param xyz Vector of length 3, corresponding to the x, y, z
#' coordinates (in mm) of a point in \code{xyz.ref.pseudo} frame of
#' reference, or 3-column matrix  or dataframe of x, y, z coordinates of several points.
#' @param vol "volume" class object.
#' @param xyz.ref.pseudo \code{ref.pseudo} in which the \code{xyz}
#' coordinate points are given. This \code{ref.pseudo} must exist in the \code{T.MAT} 
#' list. If \code{ref.pseudo} is \code{NULL} then the point with coordinates xyz 
#' is considered to be in the patient frame of reference \code{vol$ref.pseudo}.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm}, 
#' \link[espadon]{load.patient.from.dicom} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{xyz.ref.pseudo} must be equal to \code{vol$ref.pseudo} 
#' or \code{NULL}.
#' @param interpolate Boolean, default to \code{FALSE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @param verbose Boolean, default to \code{FALSE}. If \code{verbose = TRUE}, then 
#' the xyz coordinates are printed.
#' @return Returns a vector of the voxel values at the requested coordinates.
#' @seealso \link[espadon]{get.xyz.from.index}
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#' get.value.from.xyz (xyz = matrix (c (0, 0, 0, 10, 10, 10), 
#'                     ncol = 3, byrow = TRUE), vol = D)
#'
#' # isodose
#' Dmax <- max (D$vol3D.data, na.rm = TRUE)
#' idx <- which (D$vol3D.data >= (Dmax -1) & D$vol3D.data <= (Dmax - 0.3))
#' pt <- get.xyz.from.index (idx, D)
#' get.value.from.xyz (pt, vol = D, interpolate = FALSE, verbose = TRUE)

#' @export
#' @importFrom methods is
get.value.from.xyz <- function (xyz, vol, xyz.ref.pseudo = NULL, T.MAT = NULL,
                                interpolate = TRUE, verbose = FALSE){
  
  if (is.data.frame(xyz)) {
    xyz <- as.matrix (xyz)
  } else if (!is.matrix(xyz)){
    if (length(xyz)!=3) stop("xyz must be a vector of length 3 or 3-column matrix.")
    xyz <- matrix (xyz, ncol = 3) 
  }

  if (ncol(xyz)!=3) stop("xyz must be a vector of length 3 or 3-column matrix.")
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if (is.null (xyz.ref.pseudo)) xyz.ref.pseudo <- vol$ref.pseudo
  
  M.Tref <- get.rigid.M (T.MAT,src.ref=xyz.ref.pseudo, dest.ref=vol$ref.pseudo)
  if (is.null (M.Tref)) {
    warning ("vol$ref.pseudo and xyz.ref.pseudo are different.")
    return (NULL)
  }
  pt_ <- cbind (xyz,1)
  ijk <-(pt_  %*% t(M.Tref) %*% t(solve(vol$xyz.from.ijk)))[,1:3]
  
  if(verbose){
    colnames(pt_)  <- c("x","y","z","t")
    rownames(pt_)  <- 1:nrow(pt_)
    cat ("Value at\n")
    print(pt_[ , 1:3])
    cat ("\n")
  }
  
  return(get.value.from.ijk (ijk, vol, interpolate))
}