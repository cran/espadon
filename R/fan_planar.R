#' Creation of pyramid fan object passing through pixels of a plane.
#' @description The \code{fan.planar} function creates a "fan" class object 
#' containing, among others, the coordinates of the unit director vectors of 
#' the rays of rectangular pyramid fan. Rays are passing through all pixels of a
#' plane, which represent the pyramid basis. It is for instance useful to compute 
#' rt-image.
#' @param vol "volume" class object.
#' @param k Positive number specifying the plane index that the rays of the fan must cross.
#' By default, k is the central plane.
#' @param origin Numeric vector, giving the xyz coordinates of the fan origin. 
#' By default \code{c (0, 0, 0)}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' @return Returns a "fan" class object (see \link[espadon]{espadon.class} 
#' for class definitions) containing, among others, 
#' \itemize{
#' \item \code{$xyz} : a matrix of 3 columns giving the xyz coordinates of the fan rays. 
#' \item \code{$local.coords} : a list of the ijkt DICOM coordinates of the crossed plane,
#' and the transfer matrix to xyz.from.ijk to compute xyz coordinates in \code{$ref.pseudo}.
#' }
#' @seealso \link[espadon]{fan.sphere}, \link[espadon]{fan.beam}, 
#' \link[espadon]{fan.to.voxel}
#' @export
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct"), dxyz = rep (step, 3))
#' fan <- fan.planar (patient$ct[[1]], origin = patient$ct[[1]]$xyz0[1,])
#' head (fan$xyz)
#' if (interactive()) {
#'   rgl::open3d ()
#'   points3d (fan$xyz)
#' }
fan.planar  <- function(vol, k = vol$k.idx[ceiling(length(vol$k.idx)/2)], 
                        origin = c(0,0,0),
                        alias = "", 
                        description = "planar fan"){
  ijkt <- as.matrix(expand.grid(0:(vol$n.ijk[1]-1), 0:(vol$n.ijk[2]-1), k,1))
  pt <-(ijkt %*% t(vol$xyz.from.ijk))[, 1:3] 
  colnames (pt) <- c("x","y","z")
  
  pt <- sweep(pt, 2, origin)
  dir <- apply(pt,2,mean) 
  dir <- dir/.fnorm(dir)
  
  ptn <- sweep(pt,1,.fnorm(pt),"/")
  
  u <- vector.product ( vol$orientation[4:6],dir)
  u <- u/.fnorm(u)
  v <- vector.product (dir, u )
  
  orientation <- c(u,v)
  origin <- matrix(as.numeric(origin),ncol=3,byrow = TRUE)
  
  fan_ <- list()
  fan_$patient <- ""
  fan_$patient.name <- ""
  fan_$patient.bd <- ""
  fan_$patient.sex <- ""
  
  fan_$file.basename <- ""
  fan_$file.dirname <- ""
  fan_$object.name <- alias
  fan_$object.alias <- alias
  
  fan_$frame.of.reference <- vol$frame.of.reference
  fan_$ref.pseudo <- vol$ref.pseudo
  fan_$modality <- "fan"
  fan_$description <- description
  fan_$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  fan_$origin <- matrix(origin,ncol=3,byrow = TRUE)
  fan_$direction <- matrix(dir,ncol=3,byrow = TRUE)
  fan_$orientation <- orientation
  
  
  fan_$xyz <- ptn
  
  fan_$local <- list()
  fan_$local [["ijkt"]]  <- ijkt
  fan_$local [["xyz.from.ijk"]]  <- vol$xyz.from.ijk
  
  
  class(fan_) <- "fan"
  return(fan_)
}