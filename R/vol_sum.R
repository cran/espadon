#' Sum of 2 volumes
#' @description The \code{vol.sum} function adds two "volume" class objects 
#' of the same grid and of the same modality.
#' @param vol1,vol2 "volume" class objects. The 
#' 2 volumes must have the same modality, and the same grid (i.e. share the same 
#' position of the voxels).
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol1$object.alias, "+", vol2$object.alias)}.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), with the same grid and modality 
#' as \code{vol1} and \code{vol2}, sum of \code{vol1} and \code{vol2}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' pat<- toy.load.patient (modality = c ( "rtdose"), dxyz = rep (step, 3),
#'                         beam.nb = 3)
#' 
#' # Double dose
#' D <- vol.sum (pat$rtdose[[1]], pat$rtdose[[1]])
#' pat$rtdose[[1]]$max.pixel
#' D$max.pixel

#' @export
#' @importFrom methods is
vol.sum<- function (vol1, vol2, alias = "", description = NULL) {
  
  if (!is (vol1, "volume")) {
    warning ("vol1 should be a volume class object.")
    return (NULL)
  }
  if (!is (vol2, "volume")) {
    warning ("vol2 should be a volume class object.")
    return (NULL)
  }
  
  if (is.null(vol1)|is.null(vol2)) {
    warning ("vol1 or vol2 must not be null.")
    return (NULL)
  }
  
  if (vol1$modality!=vol2$modality) {
    warning ("both volumes must have same modality.")
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
    warning ("both volumes must share same grid.")
    return (NULL)
  }
  if (is.null(description)) description <-  paste (vol1$object.alias, "+", vol2$object.alias)
  V <- vol.copy (vol1, alias = alias, description = description)
  
  V$vol3D.data <- vol1$vol3D.data + vol2$vol3D.data
  
  if (any(!is.na(V$vol3D.data))){
    V$min.pixel <- min ( V$vol3D.data, na.rm=TRUE)
    V$max.pixel <- max ( V$vol3D.data, na.rm=TRUE)
  } else {
    V$min.pixel <- NA
    V$max.pixel <- NA
  }
  
  
  return(V)
}