################################################################################
#' Binary volume opening
#' @description The \code{bin.opening} function performs a morphological operation 
#' of opening, using a sphere, on a "volume" class object of "binary" modality.
#' Opening is useful for :
#' \itemize{
#' \item removing volumes that are smaller than the \code{radius},
#' \item smoothing shapes.
#' }
#' @param vol "volume" class object, of "binary" modality.
#' @param radius Positive number, in millimeters. By default, radius = 10.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to 
#' \code{paste (vol$object.alias, "opening r =", radius)}.
#' @return Returns a "volume" class object of "binary" modality (see 
#' \link[espadon]{espadon.class} for class definitions), with the same grid as 
#' \code{vol}, in which \code{$vol3D.data} has been transformed by the opening 
#' operation.
#' @note Opening can be time consuming, try to reduce the \code{binary}
#' volume to the strict minimum, before any operations.
#' @seealso \link[espadon]{bin.dilation}, \link[espadon]{bin.erosion}, 
#' \link[espadon]{bin.closing}, \link[espadon]{add.margin}, \link[espadon]{nesting.cube}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "mr", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' MR <- patient$mr[[1]]
#'
#' # generation of a binary volume
#' b <- bin.from.vol(MR, min = 15,max = 30)
#' 
#' b.opening <- bin.opening (b, radius = step)
#' display.plane (bottom = MR, top = b, main = "Before opening", 
#'                view.coord = -20, interpolate = FALSE)
#' display.plane (bottom = MR, top = b.opening, main = "After opening", 
#'                view.coord = -20, interpolate = FALSE)
#' @export
#' @importFrom stats fft
#' @importFrom methods is
bin.opening <- function (vol, radius=10, alias = "", description = NULL) {
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if ((vol$modality!="binary")) {
    warning (" vol must be of binary modality.")
    return (NULL)
  }
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (radius < max(vol$dxyz)) {
    warning ("radius must be larger than max(vol$dxyz).")
    return (vol)
  }
  
  if (!is.null(vol$missing.k.idx)) {if (vol$missing.k.idx) message ("missing k.idx, unpredictable result.")}
  
  if (is.null(description)) description <-  paste (vol$object.alias, "opening r =", radius)
  
  Vb <- bin.inversion(vol, alias = alias, description = description)
  
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
  #####
  
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
  Vb$vol3D.data <- expansion > 0.5
  if (add.contour| any(add.post)) Vb <- .vol.border.tuning (Vb, pre.nijk = -delta,
																	 post.nijk = -delta-as.numeric(add.post))
  
  Vb$min.pixel <- all(Vb$vol3D.data)
  Vb$max.pixel <- any(Vb$vol3D.data)
  Vb$vol3D.data [na.pt] <- NA
  
  return (Vb)
}
  