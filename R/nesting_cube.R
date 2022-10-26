#' Restriction of a volume to a rectangular parallelepiped
#' @description The \code{nesting.cube} function restricts or increases
#' a volume to the rectangular parallelepiped defined by its 2 extreme vertices.

#' @param vol "volume" class object.
#' @param pt.min minimum x, y, z coordinates of the rectangular parallelepiped vertex.
#' @param pt.max maximum x, y, z coordinates of the rectangular parallelepiped vertex.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. If
#' the \code{description = NULL} (default value), it will be set to \code{vol$description}.
#' @return Returns a "volume" class object, in which 3D volume is restricted 
#' or increased to be circumscribed to the requested rectangular parallelepiped. 
#' If the created volume exceeds the initial volume, new voxels are set to \code{NA}.
#' @seealso \link[espadon]{add.margin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = "ct", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' # Calculation of the new CT restricted to the parallelepiped reduced by 10 mm.
#' pt.CT <- get.extreme.pt (CT) # extreme points of CT
#' new.pt.CT <- pt.CT + matrix (rep (c (+ 12, -12), 3), ncol = 2, byrow = TRUE)
#' new.CT <- nesting.cube (CT, new.pt.CT$min, new.pt.CT$max, alias = "new CT")
#' \dontrun{  
#' # check for change
#' display.3D.stack (CT)
#' display.3D.stack (new.CT, line.col="red")
#' }
#' @export
#' @importFrom methods is

nesting.cube <- function (vol, pt.min, pt.max, alias = "", description = NULL){
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  struct.cube <- vol$cube.idx
  struct.cube[1, c(1,4,5,8)] <- pt.min[1]
  struct.cube[1, c(2,3,6,7)] <- pt.max[1]
  struct.cube[2, c(1,2,5,6)] <- pt.min[2]
  struct.cube[2, c(3,4,7,8)] <- pt.max[2]
  struct.cube[3, 1:4] <- pt.min[3]
  struct.cube[3, 5:8] <- pt.max[3]
  
  struct.in.vol.cube.idx <- round(solve (vol$xyz.from.ijk) %*% struct.cube,6)
  rg.i <- c(floor(min(struct.in.vol.cube.idx[1,])):ceiling(max(struct.in.vol.cube.idx[1,])))
  rg.j <- c(floor(min(struct.in.vol.cube.idx[2,])):ceiling(max(struct.in.vol.cube.idx[2,])))
  rg.k <- c(floor(min(struct.in.vol.cube.idx[3,])):ceiling(max(struct.in.vol.cube.idx[3,])))
  pt000 <- c(rg.i[1],rg.j[1],rg.k[1], 1) %*% t(vol$xyz.from.ijk)
  
  flag.i <- rg.i>=0 & rg.i<vol$n.ijk[1]
  flag.j <- rg.j>=0 & rg.j<vol$n.ijk[2]
  flag.k <- !is.na (match(rg.k, vol$k.idx))
  
  if (!is.null(vol$local.gridx)) {
    vol$local.gridx <- NULL
    vol$local.gridy <- NULL
  }
  
  new.vol <- vol.copy (vol = vol, alias = alias, description = description)
  
  new.vol$n.ijk <- c (length( rg.i), length (rg.j), length (rg.k))
  new.vol$xyz.from.ijk[ ,4] <- pt000
  new.vol$k.idx <- 0:(new.vol$n.ijk[3]-1)
  new.vol$xyz0  <- matrix ((as.matrix (expand.grid (0, 0, new.vol$k.idx,1))%*% t(new.vol$xyz.from.ijk))[ ,1:3],ncol=3)
  
  new.vol$cube.idx <- matrix ( c(0, 0, 0, 1, new.vol$n.ijk[1]-1, 0, 0, 1,
                                 new.vol$n.ijk[1]-1, new.vol$n.ijk[2]-1, 0, 1, 0, new.vol$n.ijk[2]-1, 0, 1,
                                 0, 0, new.vol$n.ijk[3]-1, 1, new.vol$n.ijk[1]-1, 0, new.vol$n.ijk[3]-1, 1,
                                 new.vol$n.ijk[1]-1, new.vol$n.ijk[2]-1, new.vol$n.ijk[3]-1, 1, 0, new.vol$n.ijk[2]-1, new.vol$n.ijk[3]-1, 1),
                               nrow=4, byrow= FALSE)
  
  new.vol$vol3D.data <- array(NA, dim=new.vol$n.ijk)
  
  rg.i.to.complete <- (1:new.vol$n.ijk[1])[flag.i]
  rg.j.to.complete <- (1:new.vol$n.ijk[2])[flag.j]
  rg.k.to.complete <- (1:new.vol$n.ijk[3])[flag.k]
  
  new.vol$vol3D.data[rg.i.to.complete, rg.j.to.complete, rg.k.to.complete] <- vol$vol3D.data[rg.i[flag.i] + 1,
                                                                                             rg.j[flag.j] + 1,
                                                                                             rg.k[flag.k] + 1]
  #new.vol$vol3D.data[new.vol$vol3D.data<=1e-6] <- 0
  if (any(!is.na(new.vol$vol3D.data))){
    new.vol$min.pixel <- min(new.vol$vol3D.data, na.rm = TRUE)
    new.vol$max.pixel <- max(new.vol$vol3D.data, na.rm = TRUE)
  } else {
    new.vol$min.pixel <- NA
    new.vol$max.pixel <- NA
  }
  return(new.vol)
  
}