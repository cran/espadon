#' Image value along an axis
#' @description The \code{get.line} function calculates the value of the points 
#' of a volume \code{vol} along an axis in any direction.
#' @param vol "volume" class object.
#' @param origin Vector of x, y, z coordinates belonging to the line to extract. 
#' If \code{interpolate = FALSE}, these coordinates are replaced by the coordinates 
#' of the voxel closest to \code{origin}.
#' @param orientation Directing vector of the line in the \code{vol} frame of 
#' reference. This vector is internally normalized.
#' @param grid Vector, representing the curvilinear coordinates on the line to extract.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @return Returns a dataframe, composed of the columns \code{$x}, \code{$y}, \code{$z},
#' representing the coordinates of the points where the values are taken in 
#' \code{vol} volume, the column \code{$s} representing the curvilinear abcissa, 
#' and the column \code{$value} representing values along \code{$s}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#'
#' # Dose at maximum dose
#' origin <- get.xyz.from.index (which.max (D$vol3D.data), D)
#' display.plane (bottom = D, view.coord =  origin[3], 
#'                bg = "#0000ff")
#'
#' # Dose profile at x = origin[1] and z = origin[3].
#' l <- get.line (vol = D, origin = origin,
#'                orientation = c (0, 1, 0), interpolate = FALSE)
#' plot (l$y, l$value, type = "l")
#' grid ()
#'
#' # Dose profile at y = origin[2] and z = origin[3].
#' l <- get.line (D, origin = origin,
#'                orientation = c (1, 0, 0), interpolate = FALSE)
#' plot (l$s, l$value, type = "l")
#' grid ()

#' @export
#' @importFrom methods is
get.line <- function (vol, origin = c (0, 0, 0),
                      orientation = c(1, 0, 0),
                      grid = seq (-100, 100, 1), interpolate = TRUE){
  
  if (!is (vol, "volume")) 
    stop ("vol should be a volume class object.")

  if(is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
 
  if ((length(grid) ==0)| !is.numeric(grid)) 
    stop ("grid should be a numeric vector.")
 
  v1 <-orientation/ as.numeric(sqrt(orientation%*%orientation))
  M.grid <- matrix(rep(grid,3),byrow=F,ncol=3, dimnames = list(NULL,c("x","y","z")))
  pt <-  sweep(sweep(M.grid,2,v1,"*"),2,origin,"+")
  # pt <- as.matrix(do.call(rbind.data.frame, lapply(grid, function(ind) ind*v1+ origin)))
  # colnames(pt) <- c("x","y","z")
  ijk <- round(get.ijk.from.xyz(pt, vol),6)
  df <-as.data.frame(cbind (pt, s = grid, value= get.value.from.ijk(ijk, vol=vol,interpolate = interpolate)))
  
  return (df)
  
}