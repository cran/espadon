#' Change of frame of reference of a volume
#' @description The \code{vol.in.new.ref} function allows you to change the 
#' frame of reference of a volume.
#' @param vol "volume" class object.
#' @param new.ref.pseudo pseudonym of the frame of reference in which the volume 
#' should be located. This \code{new.ref.pseudo} must exist in the \code{T.MAT} 
#' list.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be that of the \code{vol} volume.
#' @return Returns "volume" class object in the new frame of reference 
#' \code{new.ref.pseudo}.
#' @seealso \link[espadon]{struct.in.new.ref}
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = c (4, 4, 4))
#' CT <- patient$ct[[1]]
#' CT.in.new.ref <- vol.in.new.ref (CT, patient$mr[[1]]$ref.pseudo, patient$T.MAT)

#' @export
#' @importFrom methods is

vol.in.new.ref <- function (vol, new.ref.pseudo, T.MAT, alias="",description=NULL) {
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if (new.ref.pseudo!=vol$ref.pseudo){
    if (!is (T.MAT, "t.mat")) {
      warning ("T.MAT should be a t.mat class object.")
      return (NULL)
    }}
  
  if (is.null(vol$vol3D.data))  {
    message ("vol$vol3D.data is NULL.")
  }
  M <- get.rigid.M (T.MAT, vol$ref.pseudo, new.ref.pseudo)
  if (is.null (M)) return (NULL)
  V <- vol.copy (vol, alias = alias,
                 description = description)
  V$ref.pseudo <- new.ref.pseudo
  V$frame.of.reference <- T.MAT$ref.info[T.MAT$ref.info$ref.pseudo==new.ref.pseudo, ]$ref
  V$xyz.from.ijk <- M %*% vol$xyz.from.ijk
  V$xyz.from.ijk[abs(V$xyz.from.ijk) < 1.0e-6] <- 0
  V$patient.xyz0 <- matrix ((as.matrix (expand.grid (0, 0, vol$k.idx, 1)) %*% t(V$xyz.from.ijk))[ ,1:3],ncol=3)
  V$patient.orientation <- c(V$xyz.from.ijk[1:3, 1]/vol$dxyz[1], V$xyz.from.ijk[1:3, 2]/vol$dxyz[2])
  return (V)
}