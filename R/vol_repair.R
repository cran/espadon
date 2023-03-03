#' repairing missing planes of volumes
#' @description The \code{vol.repair} function repairs missing planes in volumes.
#' @param vol "volume" class object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object. If
#' \code{description = NULL} (default value), it will be set to
#' \code{paste (vol$object.alias, "repair")}.
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class}
#' for class definitions), with no missing plane, if \code{vol} is to be repaired. 
#' Returns \code{vol} otherwise.
#' @details Missing planes at download can generate errors or unpredictible results 
#' in espadon processing. The \code{vol.repair} function detects such missing 
#' planes and recreates their value by interpolation.
#' @examples
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "mr", "rtstruct", "rtdose"),
#'                              roi.name  = "",
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' CT <- patient$ct[[1]]
#'
#' # this function removes a plane in a volume in order to simulate
#' # a dicom transfer issue
#' remove.plane <- function (vol, k) {
#'   idx <- which (vol$k.idx == k)
#'   vol$n.ijk[3] <- vol$n.ijk[3] - 1
#'   vol$xyz0 <- vol$xyz0[-idx, ]
#'   vol$k.idx <- vol$k.idx[-idx]
#'   vol$missing.k.idx <- TRUE
#'   vol$vol3D.data <- vol$vol3D.data[, , -idx]
#'   return (vol)
#' }
#'
#' # Creation of CT.damaged without the 29th slice.
#' CT.damaged<- remove.plane (CT, 29)
#' CT.fix <- vol.repair (CT.damaged)
#' 
#' # Display 
#' par (mfrow=c(3, 3))
#' for (k in 28:30) {
#' display.kplane (CT, k, main = paste("CT @ k =",k),interpolate = FALSE)
#' display.kplane (CT.damaged, k, main = "damaged CT",interpolate = FALSE)
#' display.kplane (CT.fix, k, main = "fixed CT", interpolate = FALSE)
#' }
#' @export
vol.repair <- function (vol, alias = "", description = NULL) {
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (!vol$missing.k.idx) return(vol)
  
  if (is.null(description)) description <-  paste (vol$object.alias, "repair")
  rest <- vol.copy (vol, alias = alias, description = description)
  
  rest$k.idx <- vol$k.idx[1]:vol$k.idx[length (vol$k.idx)]
  rest$missing.k.idx <- FALSE
  rest$n.ijk[3] <- length (rest$k.idx)
  
  rest$xyz0 <- matrix((as.matrix(expand.grid(0, 0, rest$k.idx,
                                             1)) %*% t(rest$xyz.from.ijk))[, 1:3], ncol = 3)
  
  
  rest$vol3D.data <- array(NA, dim=rest$n.ijk)
  m  <- match (rest$k.idx, vol$k.idx)
  fna <- is.na(m)
  missing.k <- rest$k.idx[fna]
  
  rest$vol3D.data[, ,which(!fna)] <- vol$vol3D.data
  # missing.k <- rest$k.idx[which (!rest$k.idx %in% vol$k.idx)]
  
  
  for (k in missing.k) {
    
    dummy <- vol
    dummy$n.ijk[3] <- 2
    p.up.idx <- which(dummy$k.idx > k)[1]
    p.dw.idx <- p.up.idx-1
    
    dummy$k.idx <- c(0, 1)
    dummy$vol3D.data <- vol$vol3D.data[, , c(p.dw.idx, p.up.idx)]
    
    # get.value.from.ijk(c(256, 256, 0.5), dummy)
    
    k.inter <- (k - vol$k.idx[p.dw.idx]) / (vol$k.idx[p.up.idx] - vol$k.idx[p.dw.idx])
    
    ijk.to.restore <- as.matrix (expand.grid (1:vol$n.ijk[1] - 1, 1:vol$n.ijk[2] -1, k.inter))
    img <- matrix (get.value.from.ijk (ijk.to.restore, dummy, interpolate = TRUE), nrow=dim(dummy$vol3D.data)[1])
    
    rest$vol3D.data[, , which (rest$k.idx == k)] <- img
  }
  
  if (vol$modality == "binary") rest$vol3D.data <- rest$vol3D.data >= 0.5
  
  if (any(!is.na(rest$vol3D.data))) {
    rest$min.pixel <- min(rest$vol3D.data, na.rm = TRUE)
    rest$max.pixel <- max(rest$vol3D.data, na.rm = TRUE)
  } else {
    rest$min.pixel <- NA
    rest$max.pixel <- NA
  }
  return (rest)
}