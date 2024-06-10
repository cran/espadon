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
  
  if ((!is (vol1, "volume") & !is.null(vol1)) | (!is (vol2, "volume") & !is.null(vol1))) 
    stop ("vol1 or vol2 should be volume class objects.")
  if (!is.null(vol1)){
    if (is.null(vol1$vol3D.data)) stop ("empty vol1$vol3D.data.")
  }
  if (!is.null(vol2)){
    if (is.null(vol2$vol3D.data)) stop ("empty vol2$vol3D.data.")
  }
  
  if (is.null(vol1)) return (vol2)
  if (is.null(vol2)) return (vol1)
  
  if (vol1$modality!=vol2$modality) warning ("both volumes should have same modality.")
   
  #verifier que les volumes ont le même support
  if (!grid.equal (vol1, vol2)) stop ("both volumes must share same grid.")

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
  
  if (alias=="") return(V)
  return(.set.ref.obj(V,list(vol2)))
}