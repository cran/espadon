#' Median filter on a volume
#' @description The \code{vol.median} function applies a 26-connectivity median 
#' filter on all the voxels of a "volume" class object.
#' @param vol "volume" class object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol$object.alias, "median")}.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), with the same grid and modality 
#' as \code{vol}, in which voxels are filtered by a 26-connectivity median filter.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 3
#' pat <- toy.load.patient (modality = c ("ct"), dxyz = rep (step, 3))
#' CT.median <- vol.median (pat$ct[[1]])
#' 
#' display.plane (CT.median, view.type = "sagi", view.coord = 61, 
#'                interpolate = FALSE)
#' @export
vol.median <- function (vol,  alias = "", description = NULL) {
    
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  if (is.null(description)) description <-  paste (vol$object.alias, "median")
  Md <- vol.copy (vol, alias = alias, description = description)
  

  ball_ijk <- as.matrix(expand.grid((-1):1,(-1):1,(-1):1))
  Md$vol3D.data[] <- .medianfilterC(as.numeric(vol$vol3D.data),vol$n.ijk, 
                                    0:(prod(vol$n.ijk)-1),
                                    as.numeric(ball_ijk[,1]), as.numeric(ball_ijk[,2]), 
                                    as.numeric(ball_ijk[,3]))
  

  if (any(!is.na(Md$vol3D.data))) {
    Md$min.pixel <- min(Md$vol3D.data, na.rm = TRUE)
    Md$max.pixel <- max(Md$vol3D.data, na.rm = TRUE)
  } else {
    Md$min.pixel <- NA
    Md$max.pixel <- NA
  }
  return(Md)
}


