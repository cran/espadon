#' Rainbow palette
#' @description The \code{pal.rainbow} function produces a color palette adapted 
#' to dose repesentation.
#' @param n Integer, number of colors to be in the palette
#' @param alpha Vector representing the opacity, in the range of 0 (transparent) 
#' to 1 (opaque). If \code{alpha = NULL}, all colors are opaque.
#' @return Returns a color-labeled vector of size \code{n}.
#' @examples
#' pal <- pal.rainbow (255)
#'
#' image (x = seq (0, 70, length.out = 255), y = 1,
#'        z = matrix (seq (0, 70, length.out = 255), ncol = 1),
#'        col = pal,
#'        main = "Rainbow colors")

#' @export

pal.rainbow <- function (n, alpha = seq (0.8, 0, length.out=n)) {
  if (n<2) stop("n must be >1")
  if (is.null(alpha)) alpha =rep(1,n)
  pal <- rainbow (n, s = seq (1, 0, length.out = n),
           start = 0, end = 4/6,
           alpha = alpha,
           rev = TRUE)
  attr(pal,"label") = "rainbow"
  return(pal)
  
}