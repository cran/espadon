#' Cumulative Dose Volume Histogram
#' @description The \code{histo.DVH} function calculates, for each dose, 
#' the volume receiving at least this dose.
#' @param histo "histo" class object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' If the \code{description = NULL} (default value), it will be set to 
#' \code{histo$description}.

#' @return Returns a "dvh" class object. This is a list including:
#' \itemize{ 
#' \item \code{$patient}: set to \code{histo$patient}.
#' \item \code{$patient.bd}: set to \code{histo$patient.bd}.
#' \item \code{$patient.sex}: set to \code{histo$patient.sex}.
#' \item \code{$file.basename}: set to "".
#' \item \code{$file.dirname}: set to "".
#' \item \code{$object.name}: set to "".
#' \item \code{$object.alias}: alias of the dvh object..
#' \item \code{$frame.of.reference}: set to \code{histo$frame.of.reference}.
#' \item \code{$ref.pseudo} : set to \code{histo$ref.pseudo}.
#' \item \code{$modality} : set to \code{"dvh"}.
#' \item \code{$description}: description of  the dvh object. By default, 
#' set to \code{histo$description}.
#' \item \code{$creation.date}: set to \code{Sys.Date}.
#' \item \code{$nb.MC}: set to \code{histo$nb.MC}.
#' \item \code{$breaks}: vector breakpoints.
#' \item \code{$mids}: vector of cell centers.
#' \item \code{$mids.unit}: Character string, representing the unit of the abcissa
#' of the DVH. For instance, "Gy", when \code{vol} is a rtdose.
#' \item \code{$vol}: cumulative volume receiving at least the doses defined by \code{$mids}.
#' \item \code{$pcv}: percentage of the total volume receiving at least the doses defined by \code{$mids}.
#' \item \code{$MC.vol}: cumulative volume associated with \code{histo$MC.dV_dx}, if it exists.
#' \item \code{$MC.pcv}: percentage of the total volume associated with \code{histo$MC.dV_dx}, if it exists.
#' \item \code{$MC.dxyz}: set to \code{histo$MC.dxyz}, if it exists.
#' }

#' @seealso \link[espadon]{histo.from.roi}, \link[espadon]{histo.from.bin},
#' \link[espadon]{histo.vol}, \link[espadon]{display.DVH}, \link[espadon]{display.DVH.pc}
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), 
#'                              roi.name = "gizzard", dxyz = rep (step, 3), 
#'                              beam.nb = 3)
#'
#' # Calculation of the histogram
#' H <- histo.from.roi (patient$rtdose[[1]], patient$rtstruct[[1]], 
#'                      roi.name = "gizzard", 
#'                      breaks = seq (0, 60, by = 1))
#'
#' # DVH
#' DVH <- histo.DVH (H)
#' str (DVH)

#' @export
#' @importFrom methods is

histo.DVH  <- function (histo, alias = "", 
                        description = histo$description) {
  
  if (!is (histo, "histo")) {
    warning ("histo should be a histo class object.")
    return (NULL)
  }
  
  dbin <- (histo$breaks[2:length(histo$breaks)]-histo$breaks[1:(length(histo$breaks)-1)])

  dvh <- list ()
  
  dvh$patient <- histo$patient
  dvh$patient.bd <- histo$patient.bd
  dvh$patient.sex <- histo$patient.sex
  
  dvh$file.basename <- ""
  dvh$file.dirname <- ""
  dvh$object.name <- ""
  dvh$object.alias <- alias
  
  dvh$frame.of.reference <- histo$frame.of.reference
  dvh$ref.pseudo <- histo$ref.pseudo
  dvh$modality <- "dvh"
  dvh$description <- description
  dvh$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  dvh$nb.MC <- histo$nb.MC
  dvh$step <- histo$step
  dvh$breaks <-histo$breaks
  dvh$mids <- histo$mids
  dvh$mids.unit <- histo$mids.unit
  dvh$vol <- rev (cumsum (rev (histo$dV_dx * dbin)))

  
  
  dvh$pcv <-  100 * dvh$vol/ dvh$vol[1]
  
  if (!is.null (histo$MC.counts)) {
    dvh$MC.vol <- t(apply (histo$MC.dV_dx,1, function (dV_dx) rev (cumsum (rev (dV_dx*dbin)))))
    dvh$MC.pcv <- t(apply (dvh$MC.vol, 1, function (vol) 100 * vol/ vol[1]))
    dvh$MC.dxyz <- histo$MC.dxyz
  }
  class(dvh) <- "dvh"
  return (dvh)
}