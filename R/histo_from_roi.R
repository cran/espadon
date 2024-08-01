#' Histogram according to a RoI
#' @description The \code{histo.from.roi} function calculates the histogram of
#' the volume voxels belonging to a RoI.
#' @param vol "volume" class object
#' @param struct "struct" class object.
#' @param roi.name Exact name of a RoI in \code{struct} object. By default 
#' \code{roi.name = NULL}. See Details.
#' @param roi.sname Name or part of name of a RoI in \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Value of the index of a RoI that belong to the \code{struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.patient.from.dicom} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{struct$ref.pseudo} must be equal to 
#' \code{vol$ref.pseudo}.
#' @param breaks Vector giving the breakpoints between histogram cells. If 
#' \code{breaks = NULL}, the chosen breakpoints are those used by the 
#' \link[graphics]{hist} function by default. If \code{breaks} are specified, 
#' outside values of \code{vol$vol3D.data} are not taken into account.
#' @param MC If different from \code{NULL} (default value), number of calculations 
#' that will be performed, by Monte-Carlo, by randomly moving the chosen RoI 
#' over a random distance, generated according to a normal distribution with
#' mean translation defined by \code{offset} and
#' standard deviation \code{sd}.
#' @param sd Vector representing the standard deviation of distances in the 3 
#' directions x, y and z.
#' @param offset Vector representing the translation of the RoI in the 3 
#' directions x, y and z.
#' @param over.sampling.factor Strictly positive integer, or a vector of 3 strictly 
#' positive integers, default to 1. Defined to oversample grids of \code{vol}.
#' Oversampling can be very time consuming.
#' @param alias Character string, \code{$alias} of the created object
#' @param description Character string, describing the the created object. If 
#' the \code{description = NULL} (default value), it will be set to 
#' \code{struct$roi.info$roi.pseudo[roi.idx]}

#' @return Returns "histo" class object. This is a list including:
#' \itemize{
#' \item \code{$alias}: alias of the histo object.
#' \item \code{$description}: description of the histo object.
#' \item \code{$breaks}: vector breakpoints
#' \item \code{$mids}: vector of cell centers.
#' \item \code{$mids.unit}: Character string, representing the unit of the abcissa
#' of the histogram. For instance, "Gy", when \code{vol} is a rtdose.
#' \item \code{counts}: count of voxels whose value is included in the limits 
#' defined by \code{$breaks}.
#' \item \code{dV_dx}: differential histogram, expressed in cm3 by voxel units, 
#' at each \code{$mids}.
#' \item \code{MC.counts}: array of \code{MC} rows. Each row \code{i} represents 
#' the histogram of the voxels
#' contained in the RoI, whose points have been shifted by \code{$MC.dxyz[i,]}.
#' \item \code{MC.dV_dx}: array of \code{MC} rows. Each row \code{i} represents 
#' the differential histogram
#' of the voxels contained in the RoI, the points of which have been shifted by 
#' \code{$MC.dxyz[i,]}.
#' \item \code{MC.dxyz}: array of \code{MC} rows, representing the offset applied 
#' to the RoI.
#' }
#' @details \code{roi.name}, \code{roi.sname}, and \code{roi.idx} must select
#' only one RoI.
#' @note Using Monte-Carlo can be time consuming for large RoI.
#' @note If you only want the result just for a translation, use the arguments 
#' \code{MC = 1}, \code{sd = 0} and \code{offset =} desired translation vector.
#' @seealso \link[espadon]{histo.vol}, \link[espadon]{histo.from.bin},
#' \link[espadon]{display.histo}, \link[espadon]{display.dV_dx}
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for better
#' # result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), 
#'                              roi.name = "gizzard", dxyz = rep (step, 3), 
#'                              beam.nb = 3)
#'
#' # Calculation of the histogram
#' H <- histo.from.roi (patient$rtdose[[1]], patient$rtstruct[[1]], 
#'                      roi.name = "gizzard", 
#'                      breaks = seq (0, 60, by = 1))
#' str (H)

#' @export

#' @importFrom stats rnorm sd
#' @importFrom methods is
histo.from.roi  <- function (vol, struct, roi.name = NULL, roi.sname = NULL, 
                             roi.idx = NULL, T.MAT = NULL,
                             breaks = NULL,
                             MC = NULL, sd = c (1, 1, 1), offset = c (0, 0, 0), 
                             over.sampling.factor = 1,
                             alias = "", description = NULL) {
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if (!is (struct, "struct")){
    warning ("struct should be a struct class object.")
    return (NULL)
  }
  
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (length (roi.idx) != 1) {
    warning ("multiple names or no names forbidden.")
    return (NULL)
  }
  
  if (!((length(over.sampling.factor)==1 | length(over.sampling.factor)==3)) | 
      !all(abs(as.integer(over.sampling.factor)) == over.sampling.factor) | 
      any(over.sampling.factor < 0)){
    warning ("over.sampling.factor should be an integer >0 or vector of length 3 of integers >0")
    return (NULL)
  }
  
  if (length(over.sampling.factor)==1) over.sampling.factor=rep(over.sampling.factor,3)
  
  vol.alias <- vol$object.alias
  vol.info <- vol$object.info
  vol <- nesting.roi(vol, struct, roi.idx=roi.idx, T.MAT=T.MAT, vol.restrict = TRUE, 
                     xyz.margin = abs(offset)+5*abs(sd)+ abs(vol$dxyz))
  
  if (all(over.sampling.factor==c(1,1,1))) {
    vol_o <- vol
  } else {
    vol_o <- vol.oversampling (vol, fact.ijk = over.sampling.factor)
  }

  vol_o$object.alias <- vol.alias
  vol_o$object.info <- vol.info
  bin.vol <- bin.from.roi (vol_o, struct, roi.idx=roi.idx, T.MAT=T.MAT, modality ="weight", verbose = FALSE)
  if (is.null(bin.vol)){
    warning ("no RoI in volume.")
    return (NULL)
  }
  bin.vol$object.alias <- struct$object.alias
  bin.vol$object.info <- struct$object.info
  if (is.null(description)) 
    description <- paste(vol$object.alias, "histogram in ",struct$roi.info$roi.pseudo[roi.idx]) 

  H <- histo.from.bin (vol_o, bin.vol, breaks=breaks, alias = alias,
                       description = description)
  
  if (is.null (MC)) return (H)
  
  MC <- round (MC, 0)
  if  (MC <=0) return (H)
  
  H$nb.MC <- MC
  
  contour_ <- struct
  H$MC.counts <- matrix(0, nrow=MC, ncol=length(H$mids))
  H$MC.dV_dx <- matrix(0, nrow=MC, ncol=length(H$mids))
  H$MC.dxyz <- matrix(0, nrow=MC, ncol=3)
  colnames (H$MC.dxyz) <- c ("dx", "dy", "dz")
  for  (MC.idx in 1:MC){
    H$MC.dxyz [MC.idx, ] <- rnorm (3, offset, sd)
    contour_$roi.data[[roi.idx]] <- lapply (struct$roi.data[[roi.idx]],
                                            function (cont) {
                                              cont$pt <-  data.frame (t (t (cont$pt) + H$MC.dxyz [MC.idx, ]))
                                              return (cont)
                                            })
    bin.vol <- bin.from.roi (vol_o, struct=contour_, roi.idx=roi.idx, T.MAT=T.MAT, modality = "weight", verbose = FALSE)
    h<- histo.from.bin (vol_o, bin.vol, breaks=H$breaks)
    H$MC.counts[MC.idx, ] <- h$counts
    H$MC.dV_dx[MC.idx, ] <- h$dV_dx
  }
  
  class(H) <- "histo"
  return (H)
}