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
#' @param weight espadon object of the volume class, where \code{weight$vol3D.data} 
#' represents the weight of \code{vol$vol3D.data}.
#' \code{description = NULL} (default value), it will be set to \code{vol$description}.


#' @return Returns a "histo" class object. This is a list including:
#' \itemize{
#' \item \code{$patient}: set to \code{vol$patient}.
#' \item \code{$patient.name}: set to \code{vol$patient.name}.
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
histo.vol  <- function (vol, breaks = NULL, alias = "", description = NULL, weight = NULL){

  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  # 
  # if (is.null(breaks)) {
  #   H <- hist(vol$vol3D.data, plot=FALSE)
  # } else {
  #   keep <- (vol$vol3D.data >= breaks[1]) & (vol$vol3D.data < rev(breaks)[1])
  #   H <- hist(vol$vol3D.data[keep], breaks=breaks, include.lowest=TRUE, plot=FALSE)
  # }
  
  
  
  keep <- !is.na(vol$vol3D.data)
  vol3D <- vol$vol3D.data[keep]
  if (!is.null(weight)) {weight3D <- weight$vol3D.data[keep]} 
  else {weight3D = rep(1, length(vol3D))}
  
  if (is.null(breaks)) {
    breaks <- seq(vol$min.pixel,vol$max.pixel, length.out = 30)
    step  <- ceiling(breaks[2]-breaks[1])
    breaks <- floor(vol$min.pixel/step)*step + (0:29) * step
    breaks <- breaks[breaks - 2*step <vol$max.pixel]
    tab <- data.frame(D=cut(vol3D, breaks = breaks, include.lowest=TRUE),
                      weight = weight3D)
  } else {
    keep <- (vol3D>= breaks[1]) & (vol3D< rev(breaks)[1])
    tab <- data.frame(D=cut(vol3D[keep], breaks = breaks, include.lowest=TRUE),
                      weight = weight3D[keep])
  }
  
  diff_step <- diff(breaks)
  byC <- as.numeric(by(tab,tab$D,function(v) sum(v$weight)))
  byC[is.na(byC)] <- 0
  tot <- sum(byC * diff_step)
  
  
  if (is.null(description)) description <- vol$description
  r <- list ()
  
  r$patient <- vol$patient
  r$patient.name <- vol$patient.name
  r$patient.bd <- vol$patient.bd
  r$patient.sex <- vol$patient.sex
  
  r$file.basename <- ""
  r$file.dirname <- ""
  r$object.name <- alias
  r$object.alias <- alias
  
  r$frame.of.reference <- vol$frame.of.reference
  r$ref.pseudo <- vol$ref.pseudo
  r$modality <- "histo"
  r$description <- description
  
  # r$acq.date <- ""
  # r$study.date <- ""
  r$creation.date <- format(Sys.Date(), "%Y%m%d")


  # r$dV_dx=H$counts /diff(H$breaks) * abs(prod (vol$dxyz))/1000
  r$nb.MC <- 0
  r$step <- step
  r$breaks <- breaks
  r$mids <- breaks[1:(length(breaks)-1)] + diff_step/2
  r$mids.unit <- vol$unit
  r$counts <- byC
  r$dV_dx <- byC /diff_step * abs(prod (vol$dxyz))/1000
  
  class(r) <- "histo"
  if (alias=="") return(r)
  return(.set.ref.obj(r,list(vol)))
}
