#' Display of a DVH
#' \loadmathjax
#' @description The \code{display.DVH} function displays the
#' Dose Volume Histogram of a "dvh" class object. Y-units are \mjeqn{cm^3}{ascii}.
#' @param dvh "dvh" class object.
#' @param add Boolean indicating whether to display the background image.
#' @param xgrid Boolean indicating the display of the x grid.
#' @param ygrid Boolean indicating the display of the y grid.
#' @param MC.plot Boolean. If \code{MC.plot = TRUE}, then \code{display.DVH} displays,
#' if they exist, the quantile zones (Prob = 0, .025, .25, .5, .75, .975, 1)
#' of MC DVH variations.
#' @param MC.col Character string, a valid palette with 4 colours corresponding 
#' to 100%, 95%, 50% and median of MC data.
#' @param ... Additional arguments xlab, ylab, xlim, ylim, main, type, col, lwd, lty and log
#' managed by the \link[base]{plot} function.
#' @return Returns a plot of the cumulative histogram included in \code{dvh}, 
#' with its median, and the quantile areas (0%-100%), (2.5%-97.5%) and (25%-75%) 
#' of the \code{dvh$vol} variations, if they exist.
#' @seealso \link[espadon]{display.DVH.pc}
#' @examples
#' # DVH without MCMC
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5 
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), roi.name = "",
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' H <- histo.from.roi (patient$rtdose[[1]], patient$rtstruct[[1]], roi.name = "ptv", 
#'                      breaks = seq (0, 60, by = 2))
#' DVH <- histo.DVH (H)
#' display.DVH (DVH)
#' 
#' \dontrun{
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), 
#'                              roi.name = "gizzard",
#'                              dxyz = c (2, 2, 2), beam.nb = 3)
#' 
#' # Calculation of the histogram
#' H <- histo.from.roi (patient$rtdose[[1]], patient$rtstruct[[1]], 
#'                      roi.name = "gizzard", 
#'                      breaks = seq (0, 60, by = 1), MC = 100)
#' 
#' # DVH
#' DVH <- histo.DVH (H)
#' display.DVH (DVH, MC.plot = TRUE, ylim = c (0, 40))
#' }
#' @export
#' @importFrom methods is
display.DVH<- function (dvh, add=FALSE, xgrid = TRUE, ygrid = TRUE,
                        MC.plot=FALSE, MC.col = grey.colors(4, rev = TRUE), ...) {
  
  if (!is (dvh, "dvh"))
    stop ("dvh should be a dvh class object.")

 
  args <- list(...)
  y <- dvh$vol
  MC.y <- dvh$MC.vol
  if ( is.null(args[['ylab']])) {
    new.ylab <- expression ("cm"^3)
  } else {
    new.ylab <- args[['ylab']]
  }
  .display.histo (histo=dvh,y=y, new.ylab=new.ylab,
                  MC.y=MC.y, MC.plot = MC.plot, MC.col = MC.col, 
				  xgrid = xgrid, ygrid = ygrid, add=add,...)
  
}