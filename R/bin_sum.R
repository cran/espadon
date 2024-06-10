################################################################################
#' Sum of two binaries
#' @description The \code{bin.sum} function creates a "volume" class object 
#' of "binary" modality, representing the sum (logical \code{OR}) of two binary 
#' objects.
#' @param vol1,vol2 "volume" class objects of "binary" modality.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol1$object.alias, "+", vol2$object.alias)}.
#' @return Returns a "volume" class object of "binary" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol1} and \code{vol2}, sum of \code{vol1} and \code{vol2}.

#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "eye",
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' z.leye<- S$roi.info$Gz[S$roi.info$roi.pseudo == "lefteye"]
#'
#'
#' # 'left eye' et 'right eye' binaries
#' bin.left.eye <- bin.from.roi (CT, struct = S, roi.sname = "lefteye",
#'                               alias = "left eye", verbose = FALSE)
#' bin.right.eye <- bin.from.roi (CT, struct = S, roi.name = "righteye",
#'                                alias = "right eye", verbose = FALSE)
#' bin.eyes <- bin.sum (bin.left.eye, bin.right.eye, alias = "eyes")
#'
#' display.plane (CT, top = bin.eyes, struct = S, roi.sname = "eye",
#'                view.coord = z.leye, legend.shift = -90 ,
#'                interpolate = FALSE)

#' @export
#' @importFrom methods is
bin.sum <- function (vol1, vol2 , alias = "", description = NULL) {
  
  if ((!is (vol1, "volume") & !is.null(vol1)) | (!is (vol2, "volume") & !is.null(vol1))) 
    stop ("vol1 or vol2 should be volume class objects.")
  if (!is.null(vol1)){
    if (vol1$modality!="binary") stop ("vol1 must be of binary modality.")
    if (is.null(vol1$vol3D.data)) stop ("empty vol1$vol3D.data.")
  }
  if (!is.null(vol2)){
    if (vol2$modality!="binary") stop ("vol2 must be of binary modality.")
    if (is.null(vol2$vol3D.data)) stop ("empty vol2$vol3D.data.")
  }
  if (is.null(vol1)) return (vol2)
  if (is.null(vol2)) return (vol1)

  #verifier que les volumes ont le même support
  if (!grid.equal (vol1, vol2)) stop ("both volumes must share the same grid.")
 
  if (is.null(description)) description <-  paste (vol1$object.alias, "+", vol2$object.alias)
  Vb <- vol.copy (vol1, alias = alias, modality = "binary",
                  description = description)
  
  Vb$vol3D.data <- vol1$vol3D.data | vol2$vol3D.data
  # Vb$vol3D.data[is.na(Vb$vol3D.data)] <- FALSE
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  
  if (alias=="") return(Vb)
  return(.set.ref.obj(Vb,list(vol2)))
}
