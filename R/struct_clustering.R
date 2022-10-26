#' Clustering volumes by RoI
#' \loadmathjax
#' @description The \code{struct.clustering} function creates a new volume 
#' in which voxels are clustered and labeled by region of interest defined in an 
#' rt-struct.
#' @param vol "volume" class object.
#' @param struct "struct" class object.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} object.
#' By default \code{roi.idx = NULL}. See Details.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{struct$ref.pseudo} must be equal to 
#' \code{vol$ref.pseudo}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (struct$object.alias,"clustering")}
#' @param verbose Boolean. if \code{TRUE} (default), the RoI studied are listed.
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are
#' all set to \code{NULL}, all RoI are selected.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), of \code{"cluster"} modality. This object contains the 
#' \code{$cluster.info} field, detailing the label and volumes in \mjeqn{cm^{3}}{ascii} 
#' of the different clusters. Note that the label \code{NA} or value 0 is used for the voxels 
#' which are not contained in any RoI (air for instance).
#' @seealso \link[espadon]{get.roi.connection}
#' @export
#'
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c ("mr", "rtstruct"),  
#'                              dxyz = rep (step, 3))
#' MR <- patient$mr[[1]]
#' S <- patient$rtstruct[[1]]
#' cluster.vol <- struct.clustering (MR, S, T.MAT = patient$T.MAT, verbose = FALSE)
#' head (cluster.vol$cluster.info)
#' 
#' # Display
#' n = nrow(cluster.vol$cluster.info)
#' col = c ("#00000000", rainbow (n))
#' breaks <- seq (0, n, length.out = n+2)
#' 
#' display.plane (cluster.vol, main = "RoI clustering", view.coord = 0, 
#'                bottom.col = col, bottom.breaks = breaks, interpolate = FALSE)
#' 
struct.clustering <- function (vol, struct, roi.name = NULL, roi.sname = NULL, 
                               roi.idx = NULL, T.MAT = NULL, alias = "", 
                               description = NULL, verbose = TRUE) {
  
  
  roi.idx <- select.names (names = struct$roi.info$roi.pseudo, roi.name = roi.name, 
                           roi.sname=roi.sname, roi.idx=roi.idx)
  
  if (length (roi.idx) == 0) {
    warning ("no roi selected.")
    return (NULL)
  }
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if (!is (struct, "struct")) {
    warning ("struct should be a struct class object.")
    return (NULL)
  }
  

  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if(is.null(struct$roi.data)){
    warning ("empty roi.data.")
    return (NULL)
  }
  
  if (is.null (T.MAT)){
    if (vol$ref.pseudo!=struct$ref.pseudo) {
      warning ("different ref.pseudo. Enter T.MAT")
      return (NULL)
    }
  }
  
  o <- order (struct$roi.info$vol[roi.idx], decreasing = TRUE)
  
  if (is.null(description)) description <- paste (struct$object.alias,"clustering")
  b <- vol.copy(vol, alias = alias, modality = "cluster", 
                description = description, number = NULL)
  b$modality <- "cluster"
  b$vol3D.data <- array ("0", dim(b$vol3D.data))
  
  for (i in roi.idx[o]) {
    if (verbose) cat (paste(struct$roi.info$roi.pseudo[i],"\n"))
    b_ <- bin.from.roi (vol, struct, roi.idx = i, T.MAT = T.MAT)
    selected.idx <- which (b_$vol3D.data)
    if (length(selected.idx))
      b$vol3D.data[selected.idx] <- paste (b$vol3D.data[selected.idx], i, sep="|")
  }
  
  zz.l <- unique (as.vector (b$vol3D.data))
  
  zz.n <- sapply (zz.l, function (s) {
    s_ <- as.numeric (strsplit (s, "[|]")[[1]])
    if (length (s_) != 1) {
      return (paste (struct$roi.info$roi.pseudo[s_[-1]], sep="|", collapse = "|"))
    }
    return ("NA")
  })
  names (zz.n) <- NULL
  
  zz <- match (b$vol3D.data, zz.l) - 1
  
  b$vol3D.data <- array (zz, dim = dim(b$vol3D.data))
  b$max.pixel=max(b$vol3D.data, na.rm =T)
  b$min.pixel=min(b$vol3D.data, na.rm =T)
   
  zz.t <- (hist(zz, breaks = (b$min.pixel:(b$max.pixel+1))-0.5, plot = FALSE)$counts) *
    prod (b$dxyz) / 1000.0
  
  b$cluster.info <- data.frame (label = zz.n, value = (1:length (zz.n)) - 1, 
                                volume.cc = as.numeric (zz.t))
  
  return (b)
}