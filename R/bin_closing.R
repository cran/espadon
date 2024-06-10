################################################################################
#' Binary volume closing
#' @description The \code{bin.closing} function performs a morphological 
#' operation of closing, using a sphere, on a "volume" class object of "binary" modality.
#' Closing is useful for :
#' \itemize{
#' \item filling holes that are smaller than the \code{radius},
#' \item merging two shapes close to each other.
#' }
#' @param vol "volume" class object, of "binary" modality
#' @param radius Positive number, in millimeters. By default, radius = 10.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol$object.alias, "closing r =", radius)}.
#' @return Returns a "volume" class object of "binary" modality 
#' (see \link[espadon]{espadon.class} for class definitions),  with the same grid 
#' as \code{vol}, in which \code{$vol3D.data} has been transformed  by the closing 
#' operation.
#' @note Closing can be time consuming, try to reduce the \code{binary}
#' volume to the strict minimum, before any operations.
#' @seealso \link[espadon]{bin.dilation}, \link[espadon]{bin.erosion}, 
#' \link[espadon]{bin.opening}, \link[espadon]{add.margin}, \link[espadon]{nesting.cube}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "mr", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' MR <- patient$mr[[1]]
#'
#' # generation of a binary volume
#' b <- bin.from.vol(MR, min = 15, max = 30)
#' 
#' b.closing <- bin.closing (b, radius = step)
#' display.plane (bottom = MR, top = b, main = "Before closing", 
#'                view.coord = -20, interpolate = FALSE)
#' display.plane (bottom = MR, top = b.closing, main = "After closing", 
#'                view.coord = -20, interpolate = FALSE)
#' @export
#' @importFrom stats fft
#' @importFrom methods is
bin.closing <- function (vol, radius=10, alias = "", description = NULL) {
  
  if (is.null(vol)) return (NULL)
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if ((vol$modality!="binary")) stop ("vol must be of binary modality.")
  if(is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")

  radius <- rep(radius,3)[1:3]
  orientation <- matrix(vol$orientation,ncol=2, byrow=FALSE)
  M.cut.from.refpseudo <- solve(as.matrix(cbind(orientation[,1], orientation[,2],
                                                vector.product(orientation[,1], orientation[,2])),
                                          dimnames = list(NULL,NULL)))
  radius <- abs(as.vector(radius  %*%t(M.cut.from.refpseudo)))
  if (all(radius -vol$dxyz<0)) {
    warning ("At least one radius in the xyz direction must be larger than vol$dxyz.")
    return (vol)
  }
  
  if (!is.null(vol$missing.k.idx)) {if(vol$missing.k.idx) message ("missing k.idx, unpredictable result.")}
  
  if (is.null(description)) description <-  paste (vol$object.alias, "closing r =",  paste(unique(radius), collapse=" "))
  
  
  Vb <- vol.copy (vol, alias = alias, modality = "binary", description = description)
  
  na.pt <- is.na(Vb$vol3D.data)
  
  #augmentation provisoire eventuelle du volume
  delta <- ceiling(radius/abs(Vb$dxyz))
  
  n.i <- c (1:delta[1],Vb$n.ijk[1]-(1:delta[1]) + 1)
  n.i <- n.i[!is.na(match(n.i,1:Vb$n.ijk[1]))]
  n.j <- c (1:delta[2],Vb$n.ijk[2]-(1:delta[2]) + 1)
  n.j <- n.j[!is.na(match(n.j,1:Vb$n.ijk[2]))]
  n.k <- c (1:delta[3],Vb$n.ijk[3]-(1:delta[3]) + 1)
  n.k <- n.k[!is.na(match(n.k,Vb$k.idx+1))]
  
  contour.pts <- unique(rbind(expand.grid(i = n.i, 
                                          j = 1:Vb$n.ijk[2], 
                                          k = 1 :Vb$n.ijk[3]),
                              expand.grid(i = 1:Vb$n.ijk[1], 
                                          j = n.j, 
                                          k = 1 :Vb$n.ijk[3]),
                              expand.grid(i = 1:Vb$n.ijk[1], 
                                          j = 1:Vb$n.ijk[2], 
                                          k = n.k)))
  add.contour <- any(Vb$vol3D.data[as.matrix(contour.pts)], na.rm=TRUE)
  if (!add.contour) delta <- c(0,0,0)
  add.post <- !(Vb$n.ijk %% 2)
  if (add.contour | any(add.post)) Vb <- .vol.border.tuning (Vb, pre.nijk = +delta, 
                                                             post.nijk = delta+as.numeric(add.post))
  ####
  
  Vb$vol3D.data [is.na(Vb$vol3D.data)] <- FALSE
  
  kernel <- .kernel (Vb,radius)
  fft_kernel <- fft(kernel)
  
  expansion <- Re (fft (fft (Vb$vol3D.data) * fft_kernel, inverse=TRUE))
  Vb$vol3D.data <- expansion <= 0.5
  
  if (add.contour) {
    n.i <- c (1:delta[1],Vb$n.ijk[1]-(1:delta[1]) + 1)
    n.i <- n.i[!is.na(match(n.i,1:Vb$n.ijk[1]))]
    n.j <- c (1:delta[2],Vb$n.ijk[2]-(1:delta[2]) + 1)
    n.j <- n.j[!is.na(match(n.j,1:Vb$n.ijk[2]))]
    n.k <- c (1:delta[3],Vb$n.ijk[3]-(1:delta[3]) + 1)
    n.k <- n.k[!is.na(match(n.k,Vb$k.idx+1))]
    
    contour.pts <- unique(rbind(expand.grid(i = n.i, 
                                            j = 1:Vb$n.ijk[2], 
                                            k = 1 :Vb$n.ijk[3]),
                                expand.grid(i = 1:Vb$n.ijk[1], 
                                            j = n.j, 
                                            k = 1 :Vb$n.ijk[3]),
                                expand.grid(i = 1:Vb$n.ijk[1], 
                                            j = 1:Vb$n.ijk[2], 
                                            k = n.k)))
    Vb$vol3D.data [as.matrix(contour.pts)] <- FALSE
  }
  
  
  
  expansion <- Re (fft (fft (Vb$vol3D.data) * fft_kernel, inverse=TRUE))
  Vb$vol3D.data <- expansion <= 0.5
  
  if (add.contour| any(add.post)) Vb <- .vol.border.tuning (Vb, pre.nijk = -delta,
                                                            post.nijk = -delta-as.numeric(add.post))
  
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  Vb$vol3D.data [na.pt] <- NA
  
  return (Vb)
}