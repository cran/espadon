#' Voxel value at a given depth of a mesh
#' @description The \code{get.value.from.mesh} function is used to retrieve
#' the values of an object of class "volume" at the desired depth of a surface
#' described by a mesh. If the mesh corresponds to the "patient" contour,
#' the zero depth is the skin, positive depths enter the patient and negative depths exit to the outside.
#' @param mesh espadon "mesh" class object, or rgl/Rvcg "mesh3d" class object. "mesh3d" class object
#' shall an additional field "ref.pseudo" specifying the mesh frame of reference.
#' @param vol "volume" class object.
#' @param method string specifying the desired method for retrieving measurements in \code{vol}.
#' by default "point". Other methods exist, for example "disk" or "sphere. See details.
#' @param depth Numeric value, representing the depth, relative to the surface 
#' of the mesh, at which values are retrieved. 0 corresponds to the surface, 
#' positive values enter the volume used to define the mesh and negative values leave it.
#' @param radius Positive number, defining the radius of the disk or sphere, 
#' depending on the desired method.
#' @param spacing spacing of the measurement points on the disk or sphere.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.
#' @param FUN function to be applied to reduce the data ("disk" or "sphere" method) 
#' to a single value.
#' Default, median value.
#' @param ... Additional arguments passed to \code{FUN} if needed.
#' @return Returns a vector of values measured at the requested depth, with the 
#' desired method, filtered by \code{FUN}, at each vertex of the mesh.
#' @details The \code{get.value.from.mesh} function works at each vertex of the mesh.
#' It moves along the normal at that point to the desired depth.
#' \itemize{
#' \item When the method is "point", it simply retrieves the value of the volume 
#' \code{vol} specified at that point.
#' \item When the method is "disk", the values are retrieved on the disk orthogonal 
#' to the normal,with radius \code{radius}. 
#' \item When the method is "sphere", the values are retrieved inside the sphere 
#' of radius \code{radius}. 
#' }
#' For "disk" or "sphere", the measurement points are
#' spaced by \code{spacing}. For \code{radius=5} and \code{spacing=1}, "disk" and "sphere" perform
#' 78 and 523 measurements respectively.
#' In both cases, the measured values must be reduced to a single result using the 
#' \code{FUN} function. By default, espadon uses the median, but it can be 
#' provided with more complex functions to filter the data efficiently (see example below).
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "",
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' # creation of the patient mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "patient", 
#'                      verbose = FALSE)
#' mesh.patient <- mesh.from.bin (bin, alias = "patient", verbose = FALSE)
#' 
#' # density value on the skin contour, extracted from CT
#' density <- get.value.from.mesh (mesh.patient, CT ,depth = 0) 
#' 
#' if (interactive()){
#'   # Display of mesh, with RVV pal
#'   density[density < -1000] <- -1000
#'   density[density > 1000] <- 1000
#'   col <- pal.RVV(255)[cut (density, seq (-1000, 1000, length.out = 256), 
#'                            include.lowest=TRUE)]

#'   rgl::open3d ()
#'   display.3D.mesh (mesh.patient, col = col)
#' }
#' @importFrom methods is
#' @importFrom rgl addNormals
#' @export
get.value.from.mesh <- function (mesh, vol,
                                 method = c("point","disk", "sphere"),
                                 depth = 0, radius = 5, spacing = 1,
                                 T.MAT = NULL,
                                 FUN = median, ...) {
  
  if (!is (mesh, "mesh") &  !is (mesh, "mesh3d")) 
    stop ("mesh should be a mesh or rgl::mesh3d class object.")
  
  if (!is (vol, "volume")) 
    stop ("vol should be a volume class object.")
  
  
  if (is.null (mesh$ref.pseudo)) {
    ref.pseudo <- vol$ref.pseudo
  } else {
    ref.pseudo <- mesh$ref.pseudo
  }
  
  
  transfert.M <- get.rigid.M (T.MAT, src.ref=vol$ref.pseudo, dest.ref = ref.pseudo)
  
  if (is.null (transfert.M)){
    if (vol$ref.pseudo!=ref.pseudo) 
      stop ("Different ref.pseudo. Enter T.MAT.")
  } 
  
  if (is (mesh, "mesh")) mesh <- mesh$mesh
  
  if (is.null (mesh$normals)) mesh <- addNormals(mesh)
  
  radius <- abs(radius[1])
  spacing <- abs(spacing[1])
  method <- method[1]
  
  if (is.na(match(method, c("point","disk", "sphere")))){
    warning (("method is not a point or disk or sphere. It is set to 'point' value.")) 
  }
  
  if (method == "point") {
    xyz0 <- matrix (c (0, 0, 0), ncol=3)
  } else if (method == "disk") {
    r.grid <- c (-rev (seq (0, radius, by=spacing)[-1]), seq (0, radius, by=spacing))
    xyz0 <- as.matrix (expand.grid(r.grid, r.grid, 0))
    
  } else if (method == "sphere") {
    r.grid <- c (-rev (seq (0, radius, by=spacing)[-1]), seq (0, radius, by=spacing))
    xyz0 <- as.matrix (expand.grid(r.grid, r.grid, r.grid))
  }
  xyz0 <- xyz0[apply (xyz0, 1, function (v) sum (v^2) <= radius^2), ]
  
  if (!is.matrix (xyz0)) xyz0 <- matrix (xyz0, nrow=1)
  
  Mn0 <- rbind (mesh$vb[1:3, ], mesh$n[1:3, ])
  
  uv0 <- apply (Mn0, 2, function (Mn) {
    n <- Mn[4:6]
    idx <- which.max(abs (n))[1]
    dummy <- rep (0, 3)
    dummy[idx] <- n[idx]
    
    u <- vector.product (n, dummy)
    u <- u / sqrt (sum (u^2))
    v <- vector.product (n, u)
    v <- v / sqrt (sum (v^2))
    c (u, v)
  })
  
  Mnuv0 <- rbind (Mn0, uv0)
  
  xyz <- as.numeric (apply (Mnuv0, 2, function (Mnuv) {
    Mt <- matrix (c (Mnuv[7:12], Mnuv[4:6]), nrow=3, byrow=TRUE)
    dum <- t (xyz0 %*% Mt) + Mnuv[1:3] - depth * Mnuv[4:6]
    return (as.numeric (dum))
  }))
  
  xyz <- matrix (xyz, ncol=3, byrow=TRUE)
  
  v <- get.value.from.xyz (xyz, vol, xyz.ref.pseudo = ref.pseudo, T.MAT = T.MAT)
  
  V <- matrix (v, nrow=nrow (xyz0))
  
  v <- apply (V, 2, FUN, ...)
  
  return (v)
  
}