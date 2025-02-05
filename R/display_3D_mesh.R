#' 3D display of a mesh
#' @description The \code{display.3D.mesh} function performs a 3D display of a mesh.
#' @param mesh "mesh" class object, created by the \link[espadon]{mesh.from.bin} 
#' function. See \link[espadon]{espadon.class} for class definitions.
#' @param display.ref Character string. Pseudonym of the frame of reference used for display.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} or
#' \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{mesh} must be 
#' displayed in \code{display.ref = mesh$ref.pseudo}.
#' @param ... Additional arguments passed to \link[rgl]{shade3d} as \code{color}, \code{specular}, 
#' \code{alpha}...
#' @return Returns a display of \code{mesh} in the current \pkg{RGL} window if it exists, 
#' in a new window otherwise.
#' @seealso \link[espadon]{mesh.from.bin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "",
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' # creation of the patient mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "patient", verbose = FALSE)
#' mesh.patient <- mesh.from.bin (bin, alias = "patient", verbose = FALSE)
#'
#' # display of the patient mesh, with transparency
#' rgl::open3d()
#' display.3D.mesh (mesh.patient, color = "burlywood2", specular = "#404040")


#' @export
#' @importFrom rgl shade3d
#' @importFrom methods is

display.3D.mesh <- function (mesh, display.ref = mesh$ref.pseudo, T.MAT = NULL, ...) {  
    
  if (!is (mesh, "mesh")) stop ("mesh should be a mesh class object")
  M <- get.rigid.M (T.MAT, mesh$ref.pseudo, display.ref)
  if (is.null (M)) stop ("cannot display anything in selected frame of reference.")
  mesh$mesh$vb <- M %*% mesh$mesh$vb
  
  shade3d (mesh$mesh, ...)
  

}