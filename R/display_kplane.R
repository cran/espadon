#' Display of a plane of a volume
#' @description The \code{display.kplane} function displays the requested plane
#' of a "volume" class object. This function is low-level, used for example
#' in the function \link[espadon]{display.plane} with more intuitive arguments.
#' @param vol "volume" class object to display. See \link[espadon]{espadon.class} 
#' for class definitions.
#' @param k Number of the plane to display. By default \code{k} is located at
#' mid-plane of the volume.
#' @param pt00 Origin point of the displayed plane. By default \code{pt00 = c (0, 0)},
#' corresponding to the bottom left of the displayed non-flipped image.
#' @param dxy width and height of a pixel in the plane. If \code{dxy = c (1, 1)} 
#' (default) abcissa and ordinate correspond to pixel number in the plane.
#' @param col Vector, representing the color palette of the image.
#' @param breaks One of :
#' \itemize{
#' \item \code{NULL} : the minimum and the maximum value of the \code{vol} define 
#' the range.
#' \item Vector giving the breakpoints of each color. Outside values are transparent,
#' leaving the background visible, depending on \code{sat.transp}.
#' }
#' @param sat.transp Boolean. If \code{TRUE}, outside values are transparent, else 
#' set to \code{breaks} limits colors.
#' @param add Boolean indicating whether to display the background image.
#' @param main Title of the background image. If \code{main = NULL},
#' the title just indicates the value of \code{k}.
#' @param abs.lab Label of the image abcissa. By default \code{abs.lab = 'i'}.
#' @param ord.lab Label of the image ordinate. By default \code{ord.lab = 'j'}.
#' @param abs.flip Boolean defaults to \code{FALSE} flipping the horizontal axis 
#' of the background image.
#' @param ord.flip Boolean defaults to \code{FALSE} flipping the vertical axis 
#' of the background image.
#' @param bg Background color of the image. By default, this color is black.
#' @param abs.rng Vector of 2 elements indicating the minimum and maximum 
#' background image abscissa to display.
#' @param ord.rng Vector of 2 elements indicating the minimum and maximum 
#' background image ordinate to display.
#' @param interpolate Boolean, indicating whether to apply linear interpolation 
#' to the image.
#' @return Returns a display of the  \mjeqn{k^{th}}{ascii} image plane of \code{vol}.
#' @seealso \link[espadon]{display.plane}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct","mr", "rtdose"),
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' MR <- patient$mr[[1]]
#' CT <- patient$ct[[1]]
#' D <- patient$rtdose[[1]]
#' 
#' # display
#' 
#' display.kplane (CT)
#' 
#' display.kplane (MR, k = floor (length(MR$k.idx)*5/8), 
#'                 col = grey.colors (256, start = 0, end = 1),
#'                 breaks = seq (0, 500, length.out = 257), bg = "darkblue")
#' 
#' 
#' display.kplane (D, k = floor (length(D$k.idx)*3/8),
#'                 col = rainbow (256, s = seq (1, 0, length.out = 256),
#'                                start = 0, end = 4/6,
#'                                alpha = seq (0.8, 0, length.out=256),
#'                                rev = TRUE),
#'                 bg = "darkblue", ord.flip = TRUE, sat.transp = FALSE,
#'                 interpolate = FALSE)
#' 
#' display.kplane (CT, k = floor (length(CT$k.idx)/3), col = pal.RVV (1000),
#'                 breaks = seq(-1000, 1000, length.out = 1001),
#'                 bg = "darkblue", ord.flip = TRUE, interpolate = FALSE)

#' @export
#' @importFrom stats ecdf
#' @importFrom grDevices rainbow grey.colors
#' @importFrom methods is


display.kplane <- function(vol, k = vol$k.idx [ceiling (length (vol$k.idx) / 2)],
                           pt00 = c (0, 0), dxy = c (1, 1),
                           col = grey.colors (255, start = 0, end = 1), 
                           breaks = NULL,
                           sat.transp = FALSE, add = FALSE, main = NULL,
                           abs.lab = "i", ord.lab = "j", abs.flip = FALSE,
                           ord.flip = FALSE, bg="#000000", abs.rng = NULL,
                           ord.rng = NULL, interpolate = FALSE) {

  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if(is.null(vol$vol3D.data)){
    stop ("empty vol$vol3D.data.")
  }
      
  k <- k[1]
  if (!(k %in% vol$k.idx)) stop ("k plane does not exist.")
  
  
  
  # i.idx <- which((c(1, 0, 0, 0) %*% t(vol$xyz.from.ijk))!=0)[1]
  # j.idx <- which((c(0, 1, 0, 0) %*% t(vol$xyz.from.ijk))!=0)[1]
  # k.idx <- which((c(0, 0, k, 1) %*% t(vol$xyz.from.ijk))!=0)[1]
  
  map <- aperm (as.matrix(vol$vol3D.data[,,which(k==vol$k.idx)], 
                          dim =vol$n.ijk[1:2] ), perm=c (2, 1))
  
  rg.min <- pt00 - 0.5 * dxy[1:2] 
  rg.max <- pt00 + (vol$n.ijk[1:2]-0.5) * dxy[1:2] 
  # rg.min <- c (min (rg.min_[1],rg.max_[1]), min(rg.min_[2],rg.max_[2]))
  # rg.max <- c (max (rg.min_[1],rg.max_[1]), max(rg.min_[2],rg.max_[2]))
  
  
  if (!add){
    
    if (is.null(abs.rng)){
      abs.rng_ <- range(c(rg.min[1],rg.max[1]))
    } else {
      abs.rng_ <- abs.rng
    }
    if (abs.flip) abs.rng_[1:2] <- abs.rng_[2:1]
    
    if (is.null(ord.rng)){
      ord.rng_ <- range(c(rg.min[2],rg.max[2]))
    } else {
      ord.rng_ <- ord.rng
    }
    if (ord.flip) ord.rng_[1:2] <- ord.rng_[2:1]
    
    # par(mar=c(4,4,4,4))
    if (is.null(main)) main = paste("@ k =",k)
    
    plot (abs.rng_, ord.rng_, type="n", xlim = abs.rng_, ylim = ord.rng_, asp=1,
          main=main ,xlab=abs.lab, ylab = ord.lab, xaxs="i", yaxs="i")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
    
  }
  if (!is.na(vol$min.pixel) & !is.na(vol$max.pixel)) {
    # abs.left <- rg.min[1] - dxy[1]/2.0
    # ord.bottom <- rg.min[2] - dxy[2]/2.0
    # abs.right <- rg.max[1] + dxy[1]/2.0
    # ord.top <- rg.max[2] + dxy[2]/2.0
    abs.left <- rg.min[1] 
    ord.bottom <- rg.min[2] 
    abs.right <- rg.max[1]
    ord.top <- rg.max[2]
    
    imin <- as.numeric(vol$min.pixel) #r[1]
    imax <- as.numeric(vol$max.pixel) #r[2]
    
    if (!is.null(breaks) & !any(is.na(breaks))) {
      if (length(breaks)!= length (col) + 1) stop("length(breaks) must be equal to length(col) + 1.")
      
    } else {
      breaks <- .pixel.scale (imin,imax,length(col))
    }
    
    
    if(!any(is.na(breaks))){
      if (!sat.transp) {
        # map[which(map <= imin)] <- imin
        # map[which(map >= imax)] <- imax
        breaks[1] <- imin - 1
        breaks[length (breaks)] <- imax + 1
      }
      map.layer <- matrix (col[cut (as.numeric (map),breaks, include.lowest=TRUE)], nrow=nrow(map))
      
      rasterImage (map.layer[nrow(map.layer):1,], xleft = abs.left, ybottom = ord.bottom,
                   xright = abs.right, ytop = ord.top, interpolate = interpolate)
    }
    
  }
  
}