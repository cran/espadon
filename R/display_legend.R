#' Display of the RoI legend
#' @description The \code{display.legend} function displays in an image the list 
#' of requested RoI and their associated color.
#' @param struct "struct" class object.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the 
#' \code{struct} object. By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @param lwd Line thickness, defaults to 1
#' @param cex Font size, default to 1.
#' @param displayed.roi.name Vector. If different from \code{NULL}, it represents
#' the replacement names of selected RoI if needed.  
#' @param bg color of the background.
#' @param text.col color of the legend text.
#' @details \code{roi.name}, \code{roi.sname}, and \code{roi.idx}
#' indicates the RoI to display. If all three are set to NULL, all RoI are selected.
#' @return  Returns  display of the RoI names and their associated color in the 
#' active graphics window.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = c("rtstruct"), dxyz = c (5, 5, 5))
#' S <- patient$rtstruct[[1]]
#'
#' display.legend (struct = S, roi.idx = 2:10, lwd = 2)

#' @export
#' @importFrom grDevices dev.new
#' @importFrom methods is
display.legend <- function (struct = NULL, roi.name = NULL,
                            roi.sname = NULL, roi.idx = NULL, lwd = 1, cex = 1, 
                            displayed.roi.name = NULL, bg = "black", 
                            text.col = "white") {
  
  if (!is (struct, "struct")) 
    stop ("struct should be a struct class object.")

  list.roi.idx <- NULL
  
  if (!is.null (struct)) list.roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (is.null (list.roi.idx)) return (NULL)
  
  c.name <- struct$roi.info$roi.pseudo[list.roi.idx]
  if (!is.null(displayed.roi.name))
    if(length(displayed.roi.name)==length (struct$roi.info$roi.pseudo[list.roi.idx])) c.name <- displayed.roi.name
  
  ncol=ceiling(length(list.roi.idx)*cex/22.5)
  op <- par(no.readonly = TRUE)

  par(mar=c(0,0,0,0))
  plot(c(0,0),type="n", xlim=c(0,1), ylim=c(0,1),xaxt="n", yaxt="n", xlab="", ylab="")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =bg)
  
  leg <- legend("topleft",legend =c.name , col =struct$roi.info$color[list.roi.idx],seg.len=1,
                ncol=ncol ,lwd=lwd, bty="o", cex=cex, text.col=text.col,bg=bg)

  on.exit (par(op))
  
}