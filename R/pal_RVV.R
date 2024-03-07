#' Conversion of Hounsfied Units to Realistic Volume Vizualization colors
#' @description The \code{RVV.pal} function produces a color palette where 
#' Hounsfield Units in the range -1000 HU to 1000 HU are converted into 
#' realistically colorized virtual anatomy (for use with CT), developped by 
#' \emph{J.C. Silverstein and al} \[1\]
#' @param n Integer, number of colors to be in the palette
#' @param alpha Vector representing the opacity, in the range of 0 (transparent) 
#' to 1 (opaque). If \code{alpha = NULL} (default), all colors are opaque, and 
#' no alpha channel is added to the colors.
#' @param min.col,max.col respectively the color below -1000HU (by default, black,
#'  i.e. \code{"#000000"}) and above +1000HU (by default, white, i.e. \code{"#FFFFFF"})
#' @return Returns a color-labeled vector of size \code{n}.

#' @importFrom Rdpack reprompt
#' @references 
#' \[1\] \insertRef{SILVERSTEIN2008927}{espadon}
#' @examples
#' pal <- pal.RVV (256)
#'
#' image (x = seq (-1000, 1000, length.out = 1024), y = 1,
#'        z = matrix (seq (-1100, 1100, length.out = 1024), ncol = 1),
#'        col = pal,
#'        main = "Realistic Volume Vizualization colors")

#' @export
#' @importFrom grDevices col2rgb hsv rgb2hsv rgb
#' @importFrom stats approx

pal.RVV <- function (n, alpha = NULL, min.col="#000000", max.col="#FFFFFF") {
  
  adjust.col <- function (Y, col) {
    
    get.Y <- function (rgb) {
      return (0.212656 * rgb[1] + 0.715158 * rgb[2] + 0.0721856 * rgb[3])
    }
    
    Y.col <- get.Y (col)
    
    if (Y.col > Y) {
      repeat {
        hsv.col <- rgb2hsv (col, maxColorValue = 255)
        hsv.col[3] <- hsv.col[3] - .01
        if (hsv.col[3] < 0) hsv.col[3] <- 0
        col <- col2rgb (hsv (hsv.col[1], hsv.col[2], hsv.col[3]))
        Y.col <- get.Y (col)
        if (Y.col <= Y) break
      }
    } else {
      state <- "val"
      repeat {
        hsv.col <- rgb2hsv (col, maxColorValue = 255)
        if (state == "val") {
          hsv.col[3] <- hsv.col[3] + .01
          if (hsv.col[3] > 1) hsv.col[3] <- 1
          state = "col"
        } else {
          hsv.col[2] <- hsv.col[2] - .01
          if (hsv.col[2] < 0) hsv.col[2] <- 0
          state = "val"
        }
        col <- col2rgb (hsv (hsv.col[1], hsv.col[2], hsv.col[3]))
        Y.col <- get.Y (col)
        if (Y.col >= Y) break
      }
    }
    return (rgb (col[1], col[2], col[3], maxColorValue = 255))
  }
  
  
  HU.range <- seq (-1000, 1000, length.out=n)
  HU.region <- c(-1000, -600, -400, -100, -60, 40, 80, 400, 1000)
  Y <- c(0, .45, .45, .65, .65, .1, .15, 1, 1) * 254
  r <- c(0, 194, 194, 194, 194, 102, 153, 255, 255)
  g <- c(0, 105, 105, 166, 166, 0, 0, 255, 255)
  b <- c(0, 82, 82, 115, 115, 0, 0, 255, 255)
  
  HU.Y <- approx (HU.region, Y, HU.range)$y
  r.out <- approx (HU.region, r, HU.range)$y
  g.out <- approx (HU.region, g, HU.range)$y
  b.out <- approx (HU.region, b, HU.range)$y
  
  HU.col <- cbind (r.out, g.out, b.out)
  
  HU.col.pal <- c()
  for (i in 1:length (HU.Y)) HU.col.pal <- c(HU.col.pal, adjust.col(HU.Y[i], HU.col[i, ]))
  HU.col.pal[1] <- min.col
  HU.col.pal[length (HU.col.pal)] <- max.col
  attr(HU.col.pal,"label") <-"RVV"
  if (!is.null (alpha)) {
    alpha[alpha > 1] <- 1
    alpha[alpha < 0] <- 0
    alpha.channel <- toupper (sapply (alpha, function (a) sprintf ("%02x", as.integer (255 * a))))
    return (paste (HU.col.pal, alpha.channel, sep=""))
  } else return (HU.col.pal)
}