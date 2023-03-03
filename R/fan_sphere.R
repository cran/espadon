#' Creation of spherical fan object.
#' @description The \code{fan.sphere} function  creates a "fan" class object 
#' containing, among others, the coordinates of the unit director vectors of 
#' the rays of a spherical fan.
#' @param angle Positive number specifying the angle (or mean angle in case of 
#' \code{method = "random"}) between two nearest vectors.
#' @param method Requested method of fan calculation from among 'regular' and 
#' 'random'. By default, \code{method = regular}. See details. 
#' @param origin Numeric vector, giving the xyz coordinates of the fan origin. 
#' By default \code{c (0, 0, 0)}.
#' @param ref.pseudo Character string, frame of reference pseudonym of the 
#' created object.
#' @param frame.of.reference Character string, frame of reference of the 
#' created object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' @return Returns a "fan" class object (see \link[espadon]{espadon.class} 
#' for class definitions) containing, among others, 
#' \itemize{
#' \item \code{$xyz} : a matrix of 3 columns giving the xyz coordinates of the fan rays. 
#' \item \code{$local} : a matrix of 2 columns indicating the polar angle 
#' \code{theta} (rad) and the azimuthal angle \code{phi} (rad) of each ray are added.
#' }
#' @seealso \link[espadon]{fan.beam}, \link[espadon]{fan.planar}, \link[espadon]{fan.to.voxel}

#' @details The "regular" and "random" method are explained by \emph{Deserno} \strong{\[1\]}.
#' \itemize{
#' \item If \code{method = "regular"}, the returned vectors composing \code{$xyz} matrix 
#' are regularly equidistributed  at the specified angle.   
#' \item If \code{method = "random"}, the returned vectors composing \code{$xyz} matrix 
#' are randomly equidistributed at the specified angle.
#' }
#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{Deserno}{espadon}
#' @export
#' @examples
#' regular.fan <- fan.sphere (angle = 30)
#' head (regular.fan$xyz)
#' random.fan <- fan.sphere (angle = 30, method = "random")
#' head (random.fan$xyz)
#' library (rgl)
#' open3d ()
#' points3d (regular.fan$xyz)
#' open3d ()
#' points3d (random.fan$xyz)

fan.sphere <- function (angle = 1, method = c("regular","random"),
                        origin = c(0,0,0),
                        ref.pseudo = "", 
                        frame.of.reference = "",
                        alias = "", 
                        description = "spherical fan"){
  if (angle==0) stop("angle must not be 0")
  method <- method[1]
  
  
  fan_ <- list()
  fan_$patient <- ""
  fan_$patient.name <- ""
  fan_$patient.bd <- ""
  fan_$patient.sex <- ""
  
  fan_$file.basename <- ""
  fan_$file.dirname <- ""
  fan_$object.name <- alias
  fan_$object.alias <- alias
  
  fan_$frame.of.reference <- frame.of.reference
  fan_$ref.pseudo <- ref.pseudo
  fan_$modality <- "fan"
  fan_$description <- description
  fan_$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  fan_$origin <-  matrix(as.numeric(origin),ncol=3,byrow = TRUE)
  fan_$direction <- matrix(c(0,0,1),ncol=3,byrow = TRUE)
  fan_$orientation <- c(1,0,0,0,1,0)
 
  
  
  class(fan_) <- "fan"
  
  if (method == "regular") {
    M <- matrix(.fansphereC(abs(angle)),ncol=5, byrow=TRUE, 
                dimnames = list(NULL, c("x","y","z", "theta", "phi")))
    fan_$xyz <- M[,1:3]
    fan_$local <- as.matrix(M[,4:5])
    
          
  } else if (method == "random"){
    a = 4*pi/((angle*pi/180.0)^2);
    MC.nb = ceiling(a)
    
    z <- runif (MC.nb,-1,1)
    phi <- runif (MC.nb,0, 2*pi)
    x = sqrt(1-z^2)*cos(phi)
    y = sqrt(1-z^2)*sin(phi)
    M <- as.matrix (data.frame(x=x, y = y, z = z, theta = acos(z), phi = phi))
    fan_$xyz <- M[,1:3]
    fan_$local <- as.matrix(M[,4:5])
  } 
  return(fan_)
}