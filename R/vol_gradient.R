#' Gradient of a volume
#' @description The \code{vol.gradient} function calculates the 3D gradient of a 
#' "volume" class object
#' @param vol "volume" class object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol$object.alias, "gradient")}.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), with the same grid and modality 
#' as \code{vol}, gradient of \code{vol}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 3
#' pat <- toy.load.patient (modality = c ("ct"), dxyz = rep (step, 3))
#' CT.gradient <- vol.gradient (pat$ct[[1]])
#' 
#' display.plane (CT.gradient, view.type = "sagi", view.coord = 61, 
#'                interpolate = FALSE)

#' @export
vol.gradient <- function (vol,  alias = "", description = NULL) {
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol1$vol3D.data.")
    return (NULL)
  }
  if (is.null(description)) description <-  paste (vol$object.alias, "gradient")
  Gr <- vol.copy (vol, alias = alias, description = description)
  
  
  Grx <- Gry <- Grz <- vol
  
  n <- vol$n.ijk[1]
  Grx$vol3D.data[2:(n - 1),,] <- vol$vol3D.data[3:n,,] - vol$vol3D.data[1:(n - 2),,]
  Grx$vol3D.data <- Grx$vol3D.data / vol$dxyz[1]
  n <- vol$n.ijk[2]
  Gry$vol3D.data[,2:(n - 1),] <- vol$vol3D.data[,3:n,] - vol$vol3D.data[,1:(n - 2),]
  Gry$vol3D.data <- Gry$vol3D.data / vol$dxyz[2]
  n <- vol$n.ijk[3]
  Grz$vol3D.data[,,2:(n - 1)] <- vol$vol3D.data[,,3:n] - vol$vol3D.data[,,1:(n - 2)]
  Grz$vol3D.data <- Grz$vol3D.data / vol$dxyz[3]
  
  Gr$vol3D.data <- sqrt (Grx$vol3D.data^2 + Gry$vol3D.data^2 + Grz$vol3D.data^2)
  Gr$vol3D.data[c(1, Gr$n.ijk[1]),,] <- 0
  Gr$vol3D.data[,c(1, Gr$n.ijk[2]),] <- 0
  Gr$vol3D.data[,,c(1, Gr$n.ijk[3])] <- 0
  
  if (any(!is.na(Gr$vol3D.data))) {
    Gr$min.pixel <- min(Gr$vol3D.data, na.rm = TRUE)
    Gr$max.pixel <- max(Gr$vol3D.data, na.rm = TRUE)
  } else {
    Gr$min.pixel <- NA
    Gr$max.pixel <- NA
  }
  return (Gr)
  
}

