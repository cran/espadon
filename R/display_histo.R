#' Display of the counts of a histogram
#' @description The \code{display.histo} function displays the counts of
#' "histo" class object.
#' @param histo "histo" class object.
#' @param add Boolean indicating whether to display the background image.
#' @param xgrid Boolean indicating the display of the x grid.
#' @param ygrid Boolean indicating the display of the y grid.
#' @param MC.plot Boolean. If \code{MC.plot = TRUE}, then \code{display.histo} displays,
#' if they exist, the quantile zones (Prob = 0, .05, .25, .5, .75, .95, 1)
#' of variations in counts.
#' @param MC.col Character string, a valid palette with 4 colours corresponding to 
#' 100%, 95%, 50% and median of MC data.
#' @param ... Additional arguments xlab, ylab, xlim, ylim, main, type, col, lwd, lty and log
#' managed by the \link[base]{plot} function.
#' @seealso \link[espadon]{display.dV_dx}.
#' @return Returns a plot of the counts included in \code{histo}, with its median, 
#' and the quantile areas (0%-100%), (2.5%-97.5%) and (25%-75%) of the \code{histo$counts} 
#' variations, if they exist.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 3
#' patient <- toy.load.patient (modality = "ct", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' # histogram and display
#' H <- histo.vol (CT, breaks = seq (3, ceiling (CT$max.pixel), 1), 
#'                 alias = "CT_hist")
#' display.histo (H, log = "y", lwd = 2)

#' @export
#' @importFrom methods is
display.histo <- function (histo, add = FALSE, xgrid = TRUE, ygrid = TRUE,
                          MC.plot=FALSE, MC.col = grey.colors (4, rev = TRUE), 
                          ...) {
  
  if (!is (histo, "histo")) 
    stop ("display.histo ERROR : histo should be a histo class object.")
  
  
  
  args <- list(...)
  y <- histo$counts
  MC.y <- histo$MC.counts
  if ( is.null(args[['ylab']])) {
    new.ylab <-"counts"
  } else {
    new.ylab <- args[['ylab']]
  }
  .display.histo (histo=histo,y=y,new.ylab=new.ylab,MC.y=MC.y, MC.plot = MC.plot,
                  MC.col = MC.col,
                  xgrid = xgrid, ygrid = ygrid, add=add, ...)
}