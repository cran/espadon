#' Error volume
#' @description The function \code{vol.error} provides the error volume between 2 volumes.
#' @param vol,vol.ref "volume" class objects.
#' @param T.MAT "t.mat" class object to link the reference frames of \code{vol} 
#' and \code{vol.ref}. \code{T.MAT} can be created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{vol$ref.pseudo} must be equal to \code{vol.ref$ref.pseudo}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to "error".
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), with the same grid as \code{vol.ref}, and representing 
#' the error between \code{vol} and \code{vol.ref}.
#' @seealso \link[espadon]{vol.abserror}
#' @examples
#' # loading of toy-patient objects (decrease dxyz)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "sct","rtstruct"), 
#'                              roi.name = c("eye", "brain","gizzard"),
#'                              dxyz = rep (step, 3))
#' 
#' patient$ct[[1]]$description
#' patient$ct[[2]]$description
#' # Creation of the absolute error volume between ct and synthetic ct
#' vE <- vol.error (patient$ct[[2]], patient$ct[[1]], T.MAT = patient$T.MAT)
#' 
#' # Display
#' palette_vE <- colorRampPalette(c("#0000FF","#FFFFFF","#FF0000")) (100)
#' breaks_vE <- seq(floor (vE$min.pixel), ceiling (vE$max.pixel), 
#'                   length.out = 101)
#'                   
#' layout (mat = matrix(c(rep(1,6),2,2), ncol=4))                                
#' plot (vE, view.coord =61, view.type = "trans", 
#'       col = palette_vE, breaks = breaks_vE)
#' display.palette(palette_vE, breaks = breaks_vE, 
#'                 cex.axis = 1.2, main = vE$unit)
#' par(mfrow=c(1,1))
#' @export
#' @importFrom methods is
vol.error <- function (vol, vol.ref, T.MAT = NULL, alias = "", description = NULL) {
  if (is.null(description)) description <- "error"
  error.vol<- vol.regrid(vol,vol.ref,T.MAT=T.MAT,alias =alias, description= description,
                         verbose = FALSE)
  error.vol$modality <- "error"
  error.vol$vol3D.data <- error.vol$vol3D.data - vol.ref$vol3D.data
  error.vol$max.pixel <- max(error.vol$vol3D.data, na.rm = TRUE)
  error.vol$min.pixel <- min(error.vol$vol3D.data, na.rm = TRUE)
  return(error.vol)
}