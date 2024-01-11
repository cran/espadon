#' Creation of a mesh according to a binary volume
#' @description The \code{mesh.from.bin} function creates a \code{mesh} class 
#' object from a volume object of "binary" modality.
#' @param bin "volume" class object of "binary" modality.
#' @param alias Character string, \code{$alias} of the mesh defining 
#' the \code{$alias} of the created object.
#' @param tol Tolerance in mm, applied for mesh simplification. See \link[Rvcg]{vcgClean}. 
#' The default value, equal to half the smallest voxel edge, limits meshing errors.
#' @param smooth.iteration Number of iterations applied in the smoothing algorithm. 
#' See \link[Rvcg]{vcgSmooth}. 
#' @param smooth.type character: select smoothing algorithm. Available are "taubin", 
#' "laplace", "HClaplace", "fujiLaplace", "angWeight" (and any sensible abbreviations). 
#' By default, set to "taubin". See \link[Rvcg]{vcgSmooth}.
#' @param smooth.lambda numeric: parameter for Taubin smooth. See \link[Rvcg]{vcgSmooth}.
#' @param smooth.mu numeric: parameter for Taubin smooth. See \link[Rvcg]{vcgSmooth}.
#' @param smooth.delta numeric: parameter for Scale dependent laplacian smoothing 
#' (see reference below).and maximum allowed angle (in radians) for deviation 
#' between normals Laplacian (surface preserving). See \link[Rvcg]{vcgSmooth}.
#' @param verbose Boolean, by default set to \code{FALSE}. Allows you to inhibit comments.
#' @seealso \link[Rvcg]{vcgSmooth} 
#' @return Returns a "mesh" class object. This is a list including the
#' following 6 elements:
#' \itemize{
#' \item \code{$patient}: set to \code{bin$patient}.
#' \item \code{$patient.bd}: set to \code{bin$patient.bd}.
#' \item \code{$patient.name}: set to \code{bin$patient.name}.
#' \item \code{$patient.sex}: set to \code{bin$patient.sex}.
#' \item \code{$file.basename}: set to "".
#' \item \code{$file.dirname}: set to "".
#' \item \code{$object.name}: set to "".
#' \item \code{$object.alias}: set to the \code{alias} argument of the function.
#' \item \code{$frame.of.reference}: set to \code{bin$frame.of.reference}.
#' \item \code{$ref.pseudo} : set to \code{bin$ref.pseudo}.
#' \item \code{$modality} : set to \code{"mesh"}.
#' \item \code{$description}: By default, set to \code{paste (bin$object.alias, "mesh")}.
#' \item \code{$creation.date}: set to \code{Sys.Date}.
#' \item \code{$nb.faces}: set to the number of faces of the mesh.
#' \item \code{$mesh}: list of 3 elements defining the mesh:
#' \tabular{rl}{
#' \tab - \code{$vb}: array made up of the generalized coordinates (x, y, z, 1) of the vertices of the triangles.\cr 
#'  \tab There are as many columns as there are vertices.\cr
#' \tab - \code{$it}: array of the 3 indices of the vertices forming a triangle, arranged by column.\cr
#'  \tab There are as many columns as there are triangles in the mesh.\cr
#' \tab - \code{$normals}: array made up of the generalized coordinates (x, y, z, 1) of the normal vectors of each triangle.\cr
#'  \tab There are as many columns as there are vertices.\cr
#' }
#' }
#' @note To compute the mesh, all NA voxels of the binary volume \code{bin} are 
#' set to FALSE. If all voxels are equal to FALSE, the function returns the code \code{NULL}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "",
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' # creation of the patient mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "patient")
#' mesh.patient <- mesh.from.bin (bin, alias = "patient", verbose = FALSE)
#' str (mesh.patient)

#' @importFrom Rvcg vcgClean vcgIsolated vcgSmooth
#' @importFrom misc3d contour3d
#' @importFrom rgl tmesh3d
#' @importFrom methods is
#' @export
mesh.from.bin <- function (bin, alias="", tol = min(bin$dxyz)/2, 
                           smooth.iteration = 10, 
                           smooth.type = c("taubin", "laplace", "HClaplace", 
                                           "fujiLaplace", "angWeight",
                                           "surfPreserveLaplace"),
                           smooth.lambda = 0.5,
                           smooth.mu = -0.53,
                           smooth.delta = 0.1,
                           verbose = FALSE) {
  
  
  if (!is (bin, "volume")) {
    warning ("bin should be a volume class object.")
    return (NULL)
  }
  if ((bin$modality!="binary")) {
    warning ("bin must be modality binary.")
    return (NULL)
  }
  
  vol3D <- bin$vol3D.data
  vol3D[,,1] <- FALSE
  vol3D[,,dim (vol3D)[3]] <- FALSE
  vol3D[is.na(vol3D)] <- FALSE
  
  if (all(vol3D==FALSE)) return (NULL)
  mesh_ <- list ()
  
  mesh_$patient <- bin$patient
  mesh_$patient.name <- bin$patient.name
  mesh_$patient.bd <- bin$patient.bd
  mesh_$patient.sex <- bin$patient.sex
  
  mesh_$file.basename <- ""
  mesh_$file.dirname <- ""
  mesh_$object.name <- alias
  mesh_$object.alias <- alias
  mesh_$frame.of.reference <- bin$frame.of.reference
  mesh_$ref.pseudo <- bin$ref.pseudo
  mesh_$modality <- "mesh"
  mesh_$description <- paste (bin$object.alias, "mesh")
  mesh_$creation.date <- format(Sys.Date(), "%Y%m%d")

  mesh_$nb.faces <- 0
  
 
  
  mesh <- contour3d (vol3D,
                     x= 1:bin$n.ijk[1] -1, y= 1:bin$n.ijk[2] -1, z= bin$k.idx,
                     .5, fill=TRUE, smooth=TRUE, draw=FALSE)
  
  
  vertices <- matrix (as.vector (t (cbind (mesh$v1, mesh$v2, mesh$v3))), ncol=3, byrow=TRUE)
  
  vertices <- cbind (vertices, rep (1, nrow(vertices)))
  N.triangles <- nrow (vertices) / 3
  indices <- cbind (3 * (1:N.triangles) - 2, 3 * (1:N.triangles) - 1, 3 * (1:N.triangles))
  mesh <- tmesh3d (t (vertices), t (indices))
  
  if (verbose) cat ("mesh.from.bin INFO : Basic mesh cleanings\n")
  mesh <- vcgIsolated(mesh, silent = !verbose)
  mesh <- vcgClean(mesh, sel=c(1:7), tol = tol, iterate = TRUE, silent = !verbose)

  if (smooth.iteration > 0){
    if (verbose) cat ("mesh.from.bin INFO : mesh smoothing\n")
    smooth.type <- smooth.type[1]
    mesh <- vcgSmooth(mesh, iteration = smooth.iteration, type=smooth.type,
                      lambda = smooth.lambda, mu = smooth.mu,
                      delta = smooth.delta)
  }

  mesh$vb <- bin$xyz.from.ijk %*% mesh$vb
  mesh$remvert <- NULL
  mesh_$mesh <- mesh
  mesh_$nb.faces <- dim(mesh$it)[2]
  class(mesh_) <- "mesh"
  if (alias=="") return (mesh_)
  return(.set.ref.obj(mesh_,list(bin)))
  
}