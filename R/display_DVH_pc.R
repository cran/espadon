#' Display of a cumulative DVH in percent of total volume
#' @description The \code{display.DVH.pc} function displays the Dose Volume
#' Histogram of "dvh" class object. Y-units are percents of total volume.
#' @param dvh "dvh" class object. See \link[espadon]{espadon.class} for class definitions.
#' @param add Boolean indicating whether to display the background image.
#' @param xgrid Boolean indicating the display of the x grid.
#' @param ygrid Boolean indicating the display of the y grid.
#' @param MC.plot Boolean. If \code{MC.plot = TRUE}, then \code{display.DVH.pc}
#' displays, if they exist, the quantile zones (Prob = 0, .025, .25, .5, .75, .975, 1)
#' of MC DVH variations.
#' @param MC.col Character string, a valid palette with 4 colours corresponding 
#' to 100%, 95%, 50% and median of MC data.
#' @param ... Arguments xlab, ylab, xlim, ylim, main, type, col, lwd, lty and log
#' managed by the \link[base]{plot} function.
#' @return Returns a plot in percent of total volume of the cumulative histogram 
#' included in \code{dvh}, with its median, and the quantile areas (0%-100%), 
#' (2.5%-97.5%) and (25%-75%) of the \code{dvh$pcv} variations, if they exist.

#' @seealso \link[espadon]{display.DVH}
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
#' display.DVH.pc (DVH)

#' @export
#' @importFrom methods is
display.DVH.pc <- function (dvh, add = FALSE, xgrid = TRUE, ygrid = TRUE,
                            MC.plot = FALSE, MC.col = grey.colors (4, rev = TRUE), 
                            ...) {
  
  if (!is (dvh, "dvh")) stop ("dvh should be a dvh class object.")

  
  args <- list(...)
  y <- dvh$pcv
  MC.y <- dvh$MC.pcv
  if ( is.null(args[['ylab']])) {
    new.ylab <-"%"
  } else {
    new.ylab <- args[['ylab']]
  }
  ylab <- ifelse (is.null(args[['ylab']]),"%", args[['ylab']])
  .display.histo (histo=dvh,y=y,new.ylab=new.ylab,MC.y=MC.y, MC.plot = MC.plot, MC.col = MC.col,
                  xgrid = xgrid, ygrid = ygrid, add=add, ...)
}