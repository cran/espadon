################################################################################
#' Inversion of a binary
#' @description The \code{bin.inversion} function creates a "volume" class object,
#' of "binary" modality, representing the inverse (logical \code{NOT}) of another binary object.
#' @param vol "volume" class object, of "binary" modality
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste ("!", vol$object.alias, sep = "")}.
#' @return Returns a "volume" class object of "binary" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol}, inverse of \code{vol}.

#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' bin.patient <- bin.from.roi (CT, struct = S, roi.name = c ("patient"),
#'                              alias = "patient")
#' inverse.patient <- bin.inversion (bin.patient, alias = "inv (patient)")
#'
#' display.plane(CT, top = inverse.patient, interpolate = FALSE)

#' @export
#' @importFrom methods is
bin.inversion <- function (vol, alias = "", description = NULL) {
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if ((vol$modality!="binary")) {
    warning ("vol must be of binary modality.")
    return (NULL)
  }
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (is.null(description)) description <- paste ("!",vol$object.alias, sep="")
  Vb <- vol.copy (vol, alias = alias, modality = "binary",
                  description = description)
  
  Vb$vol3D.data <- !vol$vol3D.data
  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  return(Vb)
  
}