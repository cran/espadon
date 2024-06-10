#' Creation of pyramid fan object with constant angle step.
#' @description The \code{fan.beam} function creates a "fan" class object 
#' containing, among others, the coordinates of the unit director vectors of 
#' the rays of rectangular pyramid fan. Rays are uniformly distributed by angle.
#' @param alpha Positive number specifying the half-angle of the conical beam.
#' @param dalpha Positive number specifying the step of the angle between the 
#' rays of the cone beam. 
#' @param orientation Vector orientation of the pyramid base composed by the 
#' 2 orthonormal vectors coordinates.
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
#' \item \code{$local} : a matrix of 2 columns indicating the deflection angle (in rad) in 
#' the main directions defined by \code{orientation}.
#' }
#' @seealso \link[espadon]{fan.planar}, \link[espadon]{fan.sphere}, \link[espadon]{fan.to.voxel}.
#' @export
#' @examples
#' fan <- fan.beam (alpha = 30, dalpha = 1)
#' head (fan$xyz)
#' library (rgl)
#' open3d ()
#' points3d (fan$xyz)

fan.beam  <- function(alpha, dalpha,  orientation = c(0,0,1,1,0,0), 
                      origin = c(0,0,0),
                      ref.pseudo = "ref1", 
                      frame.of.reference = "",
                      alias = "", 
                      description = "beam fan"){
  
  step <- seq(dalpha,alpha,dalpha)
  dl <- tan(c(-rev(step),0, step) * pi / 180)
  pt <- as.matrix(expand.grid(dl,dl,1))
  colnames(pt) <- c("alpha.x","alpha.y","z")
  ptn <- cbind(sweep(pt,1,.fnorm(pt),"/"),0)
  colnames(ptn) <- c("x","y","z","t")
  # ptn <- unique(ptn)
  
  TM <- solve(.ref.create (orientation))
 
  
  rownames(ptn) <- NULL
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
  
  fan_$origin <- matrix(as.numeric(origin),ncol=3,byrow = TRUE)
  fan_$direction <- matrix(as.numeric(TM[1:3,3]),ncol=3,byrow = TRUE)
  fan_$orientation <- orientation
  
  
  fan_$xyz <- matrix((ptn %*% t(TM))[,1:3],ncol=3, dimnames = list(NULL,c("x","y","z")))
  # 
  # ptn <- as.data.frame(ptn)
  # ptn$theta<-  acos(ptn$z) 
  # n.xy <- .fnorm(as.matrix(ptn[,1:2]), dim=2)
  # f.xy <-  n.xy!=0
  # ptn$phi  <- rep(0,nrow(ptn))
  # ptn$phi[f.xy] <- acos (ptn$x[f.xy] / n.xy[f.xy])
  # ptn$phi[ptn$y<0]<- 2*pi - ptn$phi[ptn$y<0]
  fan_$local <- as.matrix(pt[,1:2])
  class(fan_) <- "fan"
  return(fan_)
}

.fnorm <- function(pt, dim=3){
  pt <- matrix(pt, ncol=dim)
  apply(pt,1,function(v) sqrt(sum(v^2)))
}
