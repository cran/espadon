#' Histogram according to a binary
#' @description The \code{histo.from.bin} function computes the voxels histogram 
#' of the selection defined by the binary object \code{sel.bin} of a "volume" 
#' class object.
#' @param vol "volume" class object
#' @param sel.bin "volume" class object, of \code{binary} modality
#' @param breaks Vector giving the breakpoints between histogram cells. If 
#' \code{breaks = NULL}, the chosen breakpoints are those used by the 
#' \link[graphics]{hist} function by default. If \code{breaks} are specified, 
#' outside values of \code{vol$vol3D.data} are not taken into account.
#' @param alias Character string, \code{$alias} of the created object
#' @param description Character string, describing the the created object. If the 
#' \code{description = NULL} (default value),it will be set to \code{vol$description}
#' @return Returns a "histo" class object. See \link[espadon]{histo.vol}.
#' @seealso \link[espadon]{histo.from.roi}, \link[espadon]{histo.vol},
#' \link[espadon]{display.histo}, \link[espadon]{display.dV_dx}

#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct","rtstruct"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#' bin.patient <- bin.from.roi (patient$ct[[1]], struct = patient$rtstruct[[1]],
#'                              roi.name = "patient", verbose = FALSE)
#' # ct histogram in patient  volume
#' H <- histo.from.bin (patient$ct[[1]], sel.bin = bin.patient, breaks = NULL, 
#'                      alias = "patient_hist")
#' str(H)
#' 
#' \dontrun{    
#' # Skin dose histogram 
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), roi.name = "", 
#'                              dxyz = c (2, 2, 2), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' # Creation of the skin contour of 3 mm
#' bin.patient <- bin.from.roi (D, struct = S, roi.name = "patient", 
#'                              alias = "patient", verbose = FALSE)
#' inverse.patient <- bin.inversion (bin.patient, alias = "inv (patient)")
#' expansion <- bin.dilation (inverse.patient, radius = 3, 
#'                            alias = "inv (patient) + 3")
#' contour.3mm <- bin.intersection (bin.patient, expansion, 
#'                                  alias = "contour 3 mm")
#'
#' # Dose histogram in this volume
#' H <- histo.from.bin (D, sel.bin = contour.3mm, breaks = NULL, 
#'                      alias = "Skin dose")
#' str(H)
#' }

#' @export
#' @importFrom methods is


histo.from.bin  <- function (vol, sel.bin, breaks = NULL, alias = "", 
                             description = NULL) {
  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (!is (sel.bin, "volume")) stop ("sel.bin should be a volume class object.")

  if ((sel.bin$modality!="binary" & sel.bin$modality!="weight")) stop ("sel.bin must be modality binary or weight.")

  
  vol.sel <- vol.from.bin (vol, sel.bin, alias="dum")
  vol.sel$object.alias <- vol.sel$ref.object.alias
  vol.sel$object.info <- vol.sel$ref.object.info
  weight <- NULL
  if (sel.bin$modality =="weight") weight  <- sel.bin
  return (histo.vol (vol.sel, breaks = breaks, alias=alias,
                      description = description,weight=weight))
}