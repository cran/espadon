#' Creation of a binary volume according to the voxel values of a volume
#' @description The \code{bin.from.vol} function creates a "volume" class 
#' object, of "binary" modality, in which the voxels fulfilling a condition on 
#' their value are selected.
#' @param vol "volume" class object.
#' @param min Minimum value of the selected voxel. Default to \code{-Inf}.
#' @param max Maximum value of the selected voxel. Default to \code{+Inf}.
#' @param in.selection Boolean, defaults to \code{TRUE}. If \code{in.selection = FALSE},
#'  the selected pixels are those whose value is not between \code{min} and \code{max}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (min, vol$object.alias, max, sep = "<=")} or if 
#' \code{in.selection = FALSE}, \code{paste ("!(", description, ")")}.

#' @return Returns a "volume" class object of \code{"binary"} modality, 
#' with the same grid as \code{vol}, in which the selected voxels
#' (i.e. set to TRUE) are those satisfying the following conditions :
#' \itemize{
#' \item If \code{in.selection = TRUE}, then \code{min <= vol$vol3D.data <= max}.
#' \item If \code{in.selection = FALSE}, then \code{vol$vol3D.data < min} or 
#' \code{max < vol$vol3D.data}
#' }
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 3
#' patient <- toy.load.patient (modality = "ct", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' bin.bone <- bin.from.vol (CT, min = 300, max = 3000, alias = "bone")
#' display.plane (CT, top = bin.bone, interpolate = FALSE)

#' @export
#' @importFrom methods is
bin.from.vol <- function (vol, min=-Inf, max=Inf, in.selection = TRUE,
                          alias="", description = NULL) {
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (is.null(description)) {
    description <-  paste (min, vol$object.alias, max,sep= " <= ")
    if (!in.selection) description <- paste ("!(", description, ")")
  }
  Vb <- vol.copy (vol, alias = alias, modality="binary", description=description)
  Vb$min.pixel <- FALSE
  Vb$max.pixel <- TRUE
  
  na.vox <-  is.na(Vb$vol3D.data)
  Vb$vol3D.data[!na.vox] <- (Vb$vol3D.data[!na.vox] >= min) &  (Vb$vol3D.data[!na.vox] <= max)
  
  if (!in.selection) Vb$vol3D.data <- !Vb$vol3D.data

  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- min(Vb$vol3D.data, na.rm=T)
  Vb$max.pixel <- max(Vb$vol3D.data, na.rm=T)
  Vb$vol3D.data  <- Vb$vol3D.data==T
  Vb$min.pixel <- Vb$min.pixel==T
  Vb$max.pixel <- Vb$max.pixel==T
  return (Vb)
  
}