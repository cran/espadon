#' Indices relating to the coordinates of the points
#' @description The \code{get.ijk.from.xyz} function calculates the i, j, k DICOM
#' indices of the points given in the patient's reference frame.
#' @param xyz Vector of length 3, corresponding to the x, y, z coordinates (in mm) of a point in
#' the patient's frame of reference, or 3-column matrix of x, y, z coordinates of several points.
#' @param vol "volume" class object.
#' @param verbose Boolean, default to \code{FALSE}. If \code{verbose = TRUE}, then the \code{xyz} coordinates
#' are printed.
#' @return Returns a vector or a matrix of the i, j, k DICOM indices of the x, y, z coordinate points
#' in the patient's frame of reference.
#' @note The vector or matrix is made up of real numbers. It is up to the user to make the indices as integer.
#' @note The indices of the first voxel \code{vol} are 0, 0, 0.
#' WARNING: As i,j,k are DICOM indices, they are not directly related to array indices.
#' To get the value of the \code{vol$vol3D.data}, use the function
#' \link[espadon]{get.value.from.ijk}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "ct", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' get.ijk.from.xyz (xyz = CT$xyz0[1,], vol = CT, verbose = TRUE)
#' get.ijk.from.xyz (xyz = c (1,1,1), vol = CT, verbose = TRUE)
#'
#' index <- get.ijk.from.xyz (xyz = c (1,1,1), vol = CT)
#' floor (index)
#'
#' index <- get.ijk.from.xyz (xyz = matrix (c (0,0,0,1,1,1), ncol = 3, byrow = TRUE), 
#'                            vol = CT)
#' floor (index)

#' @export
#' @importFrom methods is

get.ijk.from.xyz <- function (xyz = matrix(c(0,0,0),ncol=3), vol, verbose = FALSE) {
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  pt_ <- cbind (matrix (xyz, ncol = 3),1)
  ijk <-(pt_  %*%  t(solve(vol$xyz.from.ijk)))[,1:3]
  ijk <-matrix(ijk, ncol=3)
  if(verbose){
    colnames(pt_)  <- c("x","y","z","t")
    rownames(pt_)  <- 1:nrow(pt_)
    cat ("index of\n")
    print(pt_[ , 1:3])
    cat ("\n")
  }
  colnames(ijk)  <- c("i","j","k")
  rownames(ijk)  <- 1:nrow(ijk)
  return (ijk)
  
}