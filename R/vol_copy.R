######################################################################################
#' Creating a volume from another one
#' @description The \code{vol.copy} function creates a "volume" class object, 
#' with the same grid as the \code{vol} volume object.
#' @param vol "volume" class object, template of the created object.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param modality Character string, modality of the created volume. If 
#' \code{modality = NULL}, then the created object will have the modality of 
#' \code{vol}.
#' @param description Character string, description of the returned object. If 
#' \code{descritption = NULL}, then the created object will have the description 
#' of \code{vol}.
#' @param number number of the returned volume. If \code{number = NULL}, then 
#' the returned object will have the number of \code{vol}.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), with the same grid as \code{vol}, in which \code{$vol3D.data} 
#' is initialized to \code{NA}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' patient <- toy.load.patient (modality = "ct", roi.name = "",dxyz = c (4, 4, 4))
#' CT <- patient$ct[[1]]
#'
#' # creating a volume
#' vol.from.CT <- vol.copy (CT, alias = "ct reference")
#' str (vol.from.CT)

#' @export
#' @importFrom methods is
vol.copy <- function (vol, alias = "", modality = NULL, description = NULL, 
                      number = NULL)  {
  if (!is(vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  vol_ <- vol
  if (!is.null (modality)) vol$modality <- modality
  if (!is.null (description)) vol$description <- description
  
  vol$file.basename <- ""
  vol$file.dirname <- ""
  vol$object.name <- alias
  vol$object.alias <- alias
  vol$ref.object.alias <- NULL
  vol$object.info <- NULL
  vol$ref.object.info <- NULL
  vol$creation.date <- format(Sys.Date(), "%Y%m%d")
  if (!is.null (number)) vol$number <- number
  if (!is.null(vol$error)) vol$error <- NULL
  
  if(alias =="") return (vol)
  return(.set.ref.obj(vol,list(vol_),add=FALSE))
}
