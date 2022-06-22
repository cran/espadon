#' Coordinates of the extreme points
#' @description The \code{get.extreme.pt} function returns the x, y, z coordinates
#' of the 2 extreme voxels of the rectangular parallelepiped, containing 
#' the volume \code{vol}. These coordinates are given in the \code{ref.pseudo} 
#' frame of reference .
#' @param vol "volume" class object.
#' @param ref.pseudo Pseudonym of the frame of reference in which you want the 
#' coordinates.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{ref.pseudo} must be equal to \code{vol$ref.pseudo}.
#' @return Returns a dataframe of min and max columns, and x, y and z rows, 
#' representing the coordinates of the 2 extreme voxels of the rectangular 
#' parallelepiped.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = "ct", roi.name = "", dxyz = c (5, 5, 5))
#' CT <- patient$ct[[1]]
#'
#' # xyz extreme coordinate
#' get.extreme.pt (CT)

#' @export
#' @importFrom methods is
get.extreme.pt <- function (vol, ref.pseudo = vol$ref.pseudo, T.MAT = NULL) {
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  M <- vol$xyz.from.ijk
  if (ref.pseudo != vol$ref.pseudo) {
    M_ <- get.rigid.M(T.MAT,vol$ref.pseudo,ref.pseudo)
    if (is.null(M_)){
      warning ("different ref.pseudo. Enter T.MAT")
      return (NULL)
    }
    M <- M_ %*% M
  }
  ext_ <- (M %*% vol$cube.idx)[1:3, ]
  ext <- data.frame (min= apply (ext_,1,min),max= apply (ext_,1,max))
  
  row.names (ext) <- c("x", "y", "z")
  ext
}