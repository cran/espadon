################################################################################
#' Binary volume clustering
#' \loadmathjax
#' @description The \code{bin.clustering} function groups and labels TRUE voxels 
#' that have a 6-connectivity (i.e. sharing a common side).
#' @param vol "volume" class object, of \code{"binary"} modality
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol$object.alias,"clustering")}
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), of \code{"cluster"} modality. This object contains the 
#' \code{$cluster.info} field, detailing the label and volumes in \mjeqn{cm^{3}}{ascii} 
#' of the different clusters. Note that the label "0" is used for the background.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "ct", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' # generation of a binary volume
#' b <- bin.from.vol(CT, min = -80, max = 20)
#' 
#' # Display of the n = 3 largest volumes
#' n <- 3
#' cluster.b <- bin.clustering (b)
#' 
#' col <- c ("#00000000", rainbow (n))
#' breaks <- seq (0, n, length.out = n+2)

#' display.plane (CT, top = b, main = "Before clustering",
#'                view.coord = 50, top.col = col, top.breaks = breaks, 
#'                interpolate = FALSE)
#' display.plane (CT, top = cluster.b, main = "After clustering", 
#'                view.coord = 50, top.col = col, top.breaks = breaks, 
#'                interpolate = FALSE)

#' @export
#' @importFrom methods is
bin.clustering <- function (vol, alias="", description = NULL ) {
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if ((vol$modality!="binary")) {
    warning ("vol must be of binary modality.")
    return (NULL)
  }
  
  new.vol <- vol.copy (vol,alias = alias, modality = "cluster", description = description)
    
  na.voxel <- is.na(vol$vol3D.data)
  nb.na <- sum (na.voxel)
  vol$vol3D.data[na.voxel] <- FALSE
  vol3D <- .labelbrowser(as.vector(vol$vol3D.data), vol$n.ijk)
  vol3D[vol3D==prod(vol$n.ijk)] <- NA
  
  na.voxel_ <- is.na(vol3D)
  nb.na_ <- sum (na.voxel_)
  color <- vol3D [!na.voxel_]
  color.le <- table(color)
  ord <- order(color.le, decreasing = TRUE)
  color <- match(color,as.numeric (names(color.le))[ord])
  vol3D [!na.voxel_]  <- color
  color.le <- color.le[ord]
  label <- c("NA","bg",as.character(1: length(color.le)))
  value <-  c(NA, 0, 1:length(color.le))             
  volume.cc <- c (nb.na, nb.na_-nb.na, color.le)*prod(vol$dxyz)/1000
  vol3D[na.voxel_] <- 0
  vol3D[na.voxel] <- NA
  new.vol$vol3D.data <- vol3D
  new.vol$max.pixel=max(new.vol$vol3D.data,na.rm =T)
  new.vol$min.pixel=min(new.vol$vol3D.data,na.rm =T)
  new.vol$vol3D.data <- array(new.vol$vol3D.data, new.vol$n.ijk)
  new.vol$cluster.info <- data.frame(label= label , value=value,volume.cc=volume.cc)
  return(new.vol)
}