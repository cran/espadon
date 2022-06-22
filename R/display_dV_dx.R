#' Display of the volume density of a histogram
#' \loadmathjax
#' @description The \code{display.dV_dx} function displays the volume density
#' of a "histo" class object. Y-units are \mjeqn{cm^3.Gy^{-1}}{ascii}.
#' @param histo "histo" class object. See \link[espadon]{espadon.class} for class definitions.
#' @param add Boolean indicating whether to display the background image.
#' @param xgrid Boolean indicating the display of the x grid.
#' @param ygrid Boolean indicating the display of the y grid.
#' @param MC.plot Boolean. If \code{MC.plot = TRUE}, then \code{display.dV_dx} displays,
#' if they exist, the quantile zones (Prob = {0, .025, .25, .5, .75, .975, 1})
#' of variations in volume density.
#' @param MC.col Character string, a valid palette with 4 colours corresponding 
#' to 100%, 95%, 50% and median of MC data.
#' @param ... Arguments xlab, ylab, xlim, ylim, main, type, col, lwd, lty and log
#' managed by the \link[base]{plot} function.
#' @return Returns a plot of the differential histogram included in \code{histo}, 
#' with its median, and the quantile areas (0%-100%), (2.5%-97.5%) and (25%-75%) 
#' of the \code{histo$dv_dx} variations, if they exist.
#' @seealso \link[espadon]{display.histo}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"),
#'                              roi.name = "gizzard", dxyz = rep (step, 3), 
#'                              beam.nb = 3)
#' 
#' # Calculation of the differential histogram
#' H <- histo.from.roi (patient$rtdose[[1]], patient$rtstruct[[1]], 
#'                      roi.name = "gizzard", breaks = seq (0, 60, by = 2))
#' display.dV_dx (H, lwd = 2, col = '#00ff00', ylim = c (0,10))

#' @export
#' @importFrom methods is
display.dV_dx <- function (histo, add = FALSE, xgrid = TRUE, ygrid = TRUE,
                           MC.plot=FALSE, MC.col = grey.colors(4, rev = TRUE), ...) {
  
  if (!is (histo, "histo")){
    stop ("histo should be a histo class object.")
  }
  
  args <- list(...)
  y <- histo$dV_dx
  MC.y <- histo$MC.dV_dx
  if ( is.null(args[['ylab']])) {
    new.ylab <- bquote("cm"^3 ~ "/" ~ .(histo$mids.unit))
  } else {
    new.ylab <- args[['ylab']]
  }
  .display.histo (histo=histo,y=y,new.ylab = new.ylab,
                  MC.y=MC.y, MC.plot = MC.plot, 
				   MC.col = MC.col,
				   xgrid = xgrid, ygrid = ygrid, add=add, ...)
  
}