#' Volume class object according to binary selection
#' @description The \code{vol.from.bin} function selects a part of a "volume"
#' class object of "binary" modality which has the same grid. It is especially
#' useful to restrict voxel data in region of interest.
#' @param vol "volume" class object, containing data to restrict.
#' @param sel.bin "volume" class object, of "binary" modality.
#' \code{vol} and \code{sel.bin} must have the same grid.
#' @param alias Character string, \code{$alias} of the created object
#' @param description Character string, describing the created object.
#' If \code{description = NULL} (default value),
#' it will be set to \code{paste (vol$object.alias, "from", sel.bin$object.alias)}
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), in which non-selected voxels have the value \code{NA}, 
#' and selected voxels have the original value of \code{vol}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' patient <- toy.load.patient (modality = c ("ct", "rtstruct"), 
#'                              roi.name = "brain", dxyz = c (4, 4, 4))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' # select the brain in the volume
#' bin.brain <- bin.from.roi (vol = CT, struct = S, roi.name = "brain",
#'                            verbose = FALSE)
#' vol.brain <- vol.from.bin (CT, bin.brain)
#'# display at the center of gravity of the brain Gz
#' Gz <- S$roi.info [grep("^brain", S$roi.info$roi.pseudo),]$Gz
#' display.plane (bottom = vol.brain, view.coord = Gz, struct = S,
#'                roi.sname = "brain", bg = "#00ff00", interpolate = FALSE)

#' @export
#' @importFrom methods is
vol.from.bin <- function (vol, sel.bin, alias = "", description = NULL){
  
  if (is.null (vol)) return (NULL)
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  
  if (!is (sel.bin, "volume"))
    stop ("sel.bin should be a volume class object.")
  if (sel.bin$modality!="binary") stop ("sel.bin must be modality binary.")
  
  #verifier que les volumes ont le même support
  if (!grid.equal (vol, sel.bin)) stop ("vol and sel.bin must share same grid.")

  if (is.null(description)) description = paste (vol$object.alias, "from", sel.bin$object.alias)
  V <- vol.copy (vol, alias = alias, description = description)
  
  V$vol3D.data[] <- NA
  keep <- which(sel.bin$vol3D.data==TRUE)
  
  V$vol3D.data[keep] <- vol$vol3D.data[keep]
  if (any(!is.na(V$vol3D.data))){
    V$min.pixel <- min ( V$vol3D.data, na.rm=TRUE)
    V$max.pixel <- max ( V$vol3D.data, na.rm=TRUE)
  } else {
    V$min.pixel <- NA
    V$max.pixel <- NA
  }
  
  if (alias=="") return(V)
  return(.set.ref.obj(V, list(sel.bin)))
}
