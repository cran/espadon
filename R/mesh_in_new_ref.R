#' Change of frame of reference of a mesh
#' @description The \code{mesh.in.new.ref} function allows you to change the frame 
#' of reference of a mesh.
#' @param mesh "volume" class object.
#' @param new.ref.pseudo pseudonym of the frame of reference in which the mesh 
#' should be located. This \code{new.ref.pseudo} must exist in the \code{T.MAT} list.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be that of the \code{mesh}.
#' @return Returns "mesh" class object in the new frame of reference \code{new.ref.pseudo}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct", "mr"), 
#'                              roi.name = "", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' #creation of the patient mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "patient", verbose = FALSE)
#' mesh.patient <- mesh.from.bin (bin, alias = "patient", verbose = FALSE)
#' 
#' # mesh in the MR frame of reference
#' new.mesh <- mesh.in.new.ref (mesh.patient, patient$mr[[1]]$ref.pseudo, 
#'                              T.MAT = patient$T.MAT)
#' 
#' str (new.mesh, max.level = 2)                             

#' @export
#' @importFrom methods is

mesh.in.new.ref <- function (mesh, new.ref.pseudo, T.MAT = NULL, alias="",description=NULL) {
  
  if (is.null(mesh)) return(NULL)
  if (!is (mesh, "mesh")) stop ("mesh should be a mesh class object.")
  
  if (new.ref.pseudo!=mesh$ref.pseudo){
    if (!is (T.MAT, "t.mat")) stop ("T.MAT should be a t.mat class object.")
  }
  
  M <- get.rigid.M (T.MAT, mesh$ref.pseudo, new.ref.pseudo)
  if (is.null (M)) stop("no transfer matrix between mesh$ref.pseudo and new.ref.pseudo")
  
  mesh_ <- list(object.alias=mesh$object.alias,object.info=mesh$object.info)
  class(mesh_) <- "mesh"
  
  mesh$ref.pseudo <- new.ref.pseudo
  mesh$frame.of.reference <- T.MAT$ref.info[T.MAT$ref.info$ref.pseudo==new.ref.pseudo, ]$ref
  if (!is.null(description)) mesh$description  <- description
  mesh$file.basename <- ""
  mesh$file.dirname <- ""
  mesh$object.name <- alias
  mesh$object.alias <- alias
  mesh$object.info <- NULL
  mesh$ref.object.alias <- NULL
  mesh$ref.object.info <- NULL
  
  mesh$mesh$vb <- M %*% mesh$mesh$vb
  if (!is.null(mesh$mesh$normals)) {
    mesh$mesh$normals[4,] <- 0
    mesh$mesh$normals <- M %*% mesh$mesh$normals
    mesh$mesh$normals[4,] <- 1
  }
  # return (mesh)
  if (alias=="") return(mesh)
  return(.set.ref.obj(mesh,list(mesh_)))
}