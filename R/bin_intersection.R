################################################################################
#' Intersection of two binaries
#' @description The \code{bin.intersection} function creates a "volume" 
#' class object, of "binary" modality, representing the intersection (logical 
#' \code{AND}) of two binary objects.
#' @param vol1,vol2 "volume" class objects, of "binary" modality.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol1$object.alias, "&", vol2$object.alias)}.
#' @return Returns a "volume" class object of "binary" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol1} and \code{vol2}, intersection of \code{vol1} and \code{vol2}.
#' @export
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 3
#' patient <- toy.load.patient (modality = c("mr", "rtstruct"), 
#'                              roi.name = c("brain", "labyrinth processing unit"), 
#'                              dxyz = rep (step, 3))
#' MR <- patient$mr[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' z.brain <- S$roi.info$Gz[S$roi.info$roi.pseudo == "brain"]
#' 
#' # Try to discriminate the processing unit with binary selections
#' bin.brain <- bin.from.roi (MR, struct = S, roi.name = "brain",
#'                            alias = "brain", T.MAT = patient$T.MAT)
#' bin.pu.density <- bin.from.vol (MR, min = 160)   
#'    
#' display.plane (MR, top = bin.pu.density, display.ref = S$ref.pseudo,
#'                view.coord = z.brain, T.MAT = patient$T.MAT, 
#'                interpolate = FALSE, main = "before brain intersection")                      
#' bin.pu <- bin.intersection (vol1 = bin.pu.density, vol2 = bin.brain, 
#'                             alias = "processing unit")
#' display.plane (MR, top = bin.pu, display.ref = S$ref.pseudo,
#'                view.coord = z.brain, T.MAT = patient$T.MAT, 
#'                interpolate = FALSE, main = "after brain intersection")


#' @importFrom methods is
bin.intersection <- function (vol1, vol2, alias = "", description = NULL) {
  
  
  if (!is (vol1, "volume") | !is (vol2, "volume")){
    warning ("vol1 or vol2 should be volume class objects.")
    return (NULL)
  }
  
  if ((vol1$modality!="binary") | (vol2$modality!="binary")) {
    warning ("both volumes must be of binary modality.")
    return (NULL)
  }
  
  
  if(is.null(vol1$vol3D.data)){
    warning ("empty vol1$vol3D.data.")
    return (NULL)
  }
  if(is.null(vol2$vol3D.data)){
    warning ("empty vol2$vol3D.data.")
    return (NULL)
  }
  
  #verifier que les volumes ont le même support
  if (!grid.equal (vol1, vol2)) {
    warning ("both volumes must share the same grid.")
    return (NULL)
  }
  if (is.null(description)) description <-  paste (vol1$object.alias, "&", vol2$object.alias)
  Vb <- vol.copy (vol1, alias = alias, modality = "binary",
                  description = description)
  
  Vb$vol3D.data <- vol1$vol3D.data & vol2$vol3D.data
  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  
  return(Vb)
  
}