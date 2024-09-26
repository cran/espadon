################################################################################
#' Binary volume dilation
#' @description The \code{bin.dilation} function enlarges a 
#' "volume" class object, of \code{"binary"} modality, by means of 
#' convolution with a sphere.
#' Dilation is useful for :
#' \itemize{
#' \item filling holes that are smaller than the \code{radius},
#' \item enlarging capes,
#' \item filling narrow channels,
#' \item merging two shapes close to each other.
#' }
#' @param vol "volume" class object, of \code{"binary"} modality
#' @param radius Positive number, or xyz-vector of 3 positive numbers.By default, 
#' radius = 10.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL}
#' (default value), it will be set to \code{paste (vol$object.alias, "dilataion r =", radius)}.
#' @return Returns a "volume" class object of \code{"binary"} modality (see 
#' \link[espadon]{espadon.class} for class definitions), with 
#' the same grid as \code{vol}, in which the selected volume has been enlarged 
#' by the \code{radius}.
#' @seealso \link[espadon]{bin.erosion}, \link[espadon]{bin.opening}, 
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
#' b.dilation <- bin.dilation (b, radius = step)
#' display.plane (bottom = MR, top = b, main = "Before dilation", 
#'                view.coord = -20, interpolate = FALSE)
#' display.plane (bottom = MR, top = b.dilation, main = "After dilation", 
#'                view.coord = -20,interpolate = FALSE)

#' @export
#' @importFrom stats fft
#' @importFrom methods is
bin.dilation <- function (vol, radius=10, alias = "", description = NULL) {

  if (is.null(vol)) return (NULL)
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if ((vol$modality!="binary")) stop ("vol must be of binary modality.")
  if(is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  
  radius <- rep(radius,3)[1:3]
  M.cut.from.refpseudo <- solve(as.matrix(cbind(vol$orientation[1:3], vol$orientation[4:6],
                                                vector.product(vol$orientation[1:3], vol$orientation[4:6])),
                                          dimnames = list(NULL,NULL)))
  radius <- abs(as.vector(radius  %*%t(M.cut.from.refpseudo)))

  if (all(radius -vol$dxyz<0)) {
    warning ("At least one radius in the xyz direction must be larger than vol$dxyz.")
    return (vol)
  }
  
  if (!is.null(vol$missing.k.idx)) {if (vol$missing.k.idx) message ("missing k.idx, unpredictable result.")}
  
  if (is.null(description)) description <-  paste (vol$object.alias, "dilation r =", paste(unique(radius), collapse=" "))
  
  Vb <- vol.copy (vol, alias = alias, modality = "binary", description = description)
  na.pt <- is.na(Vb$vol3D.data)
  
  #augmentation provisoire eventuelle du volume
  delta <- ceiling(radius/abs(Vb$dxyz)) + 1
  
  ext.ijk <- apply(get.ijk.from.index(which(vol$vol3D.data>0), vol), 2, range) + 
    rbind(-delta,delta)
  
  pre.nijk <- - ext.ijk[1,]
  post.nijk <- ext.ijk[2,] - vol$n.ijk + 1 + as.numeric(!((ext.ijk[2,]-ext.ijk[1,] + 1) %% 2))
  
  add.ijk <-any(c(pre.nijk,post.nijk)!=0)
  if (add.ijk) Vb <- .vol.border.tuning (Vb, pre.nijk = pre.nijk, post.nijk = post.nijk)
  
    #####
  
  Vb$vol3D.data [is.na(Vb$vol3D.data)] <- FALSE
  
  kernel <- .kernel (Vb,radius)
  fft_kernel <- fft(kernel)
  
  expansion <- Re (fft (fft (Vb$vol3D.data) * fft_kernel, inverse=TRUE))
  Vb$vol3D.data <- expansion > 0.5

  if (add.ijk) Vb <- .vol.border.tuning (Vb, pre.nijk = -pre.nijk, post.nijk = -post.nijk) 
  
  Vb$vol3D.data [na.pt & !Vb$vol3D.data] <- NA
  Vb$min.pixel <- all(Vb$vol3D.data, na.rm = TRUE)
  Vb$max.pixel <- any(Vb$vol3D.data, na.rm = TRUE)
  
  return (Vb)
}