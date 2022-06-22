#' Histogram of a volume
#' \loadmathjax
#' @description The \code{histo.vol} function calculates the voxel values  
#' histogram of "volume" class object.
#' @param vol "volume" class object.
#' @param breaks Vector giving the breakpoints between histogram cells. If 
#' \code{breaks = NULL},
#' the chosen breakpoints are those used by the \link[graphics]{hist} function 
#' by default. If \code{breaks} are specified, outside values of \code{vol$vol3D.data} 
#' are not taken into account.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. If the 
#' \code{description = NULL} (default value), it will be set to \code{vol$description}.


#' @return Returns a "histo" class object. This is a list including:
#' \itemize{
#' \item \code{$patient}: set to \code{vol$patient}.
#' \item \code{$patient.bd}: set to \code{vol$patient.bd}.
#' \item \code{$patient.sex}: set to \code{vol$patient.sex}.
#' \item \code{$file.basename}: set to "".
#' \item \code{$file.dirname}: set to "".
#' \item \code{$object.name}: set to "".
#' \item \code{$object.alias}: alias of the histo object.
#' \item \code{$frame.of.reference}: set to \code{vol$frame.of.reference}.
#' \item \code{$ref.pseudo} : set to \code{vol$ref.pseudo}.
#' \item \code{$modality} : set to \code{"histo"}.
#' \item \code{$description}: description of the histo object.
#' \item \code{$creation.date}: set to \code{Sys.Date}.
#' \item \code{$nb.MC}: set to 0.
#' \item \code{$breaks}: vector breakpoints
#' \item \code{$mids}: vector of cell centers.
#' \item \code{$mids.unit}: Character string, representing the unit of the abcissa
#' of the histogram. For instance, "Gy", when \code{vol} is a rtdose.
#' \item \code{counts}: count of voxels whose value is included in the limits 
#' defined by \code{$breaks}.
#' \item \code{dV_dx}: differential histogram, expressed in \mjeqn{cm^3}{ascii} by voxel units, 
#' at each \code{$mids}.
#' }
#' @seealso \link[espadon]{histo.from.roi}, \link[espadon]{histo.from.bin},
#' \link[espadon]{display.histo}, \link[espadon]{display.dV_dx}
#' @examples
#' # loading of toy-patient objects
#' step <- 3
#' patient <- toy.load.patient (modality = "ct", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' # histogram and display
#' H <- histo.vol (CT, breaks = seq (3, ceiling (CT$max.pixel), 1), 
#'                 alias = "CT_hist")
#' str (H)
# display.dV_dx (H, type = "h", lwd = 2)
# display.histo (H, log = "y", lwd = 2)

#' @export
#' @importFrom graphics hist
#' @importFrom methods is
histo.vol  <- function (vol, breaks = NULL, alias = "", description = NULL){

  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (is.null(breaks)) {
    H <- hist(vol$vol3D.data, plot=FALSE)
  } else {
    keep <- (vol$vol3D.data >= breaks[1]) & (vol$vol3D.data < rev(breaks)[1])
    H <- hist(vol$vol3D.data[keep], breaks=breaks, include.lowest=TRUE, plot=FALSE)
  }
  if (is.null(description)) description <- vol$description
  r <- list ()
  
  r$patient <- vol$patient
  r$patient.bd <- vol$patient.bd
  r$patient.sex <- vol$patient.sex
  
  r$file.basename <- ""
  r$file.dirname <- ""
  r$object.name <- ""
  r$object.alias <- alias
  
  r$frame.of.reference <- vol$frame.of.reference
  r$ref.pseudo <- vol$ref.pseudo
  r$modality <- "histo"
  r$description <- description
  
  # r$acq.date <- ""
  # r$study.date <- ""
  r$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  r$nb.MC <- 0
  r$step <- H$breaks[2]-H$breaks[1]
  r$breaks=H$breaks
  r$mids=H$mids
  r$mids.unit=vol$unit
  r$counts=H$counts
  r$dV_dx=H$counts /(H$breaks[2:length(H$breaks)]-H$breaks[1:(length(H$breaks)-1)]) * abs(prod (vol$dxyz))/1000
  
  class(r) <- "histo"
  return(r)
}
