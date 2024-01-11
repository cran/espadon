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
#' @param flip Boolean defaults to \code{FALSE} flipping the horizontal axis 
#' of the background image.
#' @param flop Boolean defaults to \code{FALSE} flipping the vertical axis 
#' of the background image.
#' @param bg Background color of the image. By default, this color is black.
#' @param abs.rng Vector of 2 elements indicating the minimum and maximum 
#' background image abscissa to display.
#' @param ord.rng Vector of 2 elements indicating the minimum and maximum 
#' background image ordinate to display.
#' @param interpolate Boolean, indicating whether to apply linear interpolation 
#' to the image.
#' @param ... others argument of plot function like xaxt, yaxt...
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
                           abs.lab = "i", ord.lab = "j", flip = FALSE,
                           flop = FALSE, bg="#000000", abs.rng = NULL,
                           ord.rng = NULL, interpolate = FALSE, ...) {
  args = list(...)
  k <- k[1]
  if (!is.null(args[["abs.flip"]])) { 
    flip <- args[["abs.flip"]]
    args[["abs.flip"]] <- NULL
  }
  if (!is.null(args[["ord.flip"]]))  {
    flop <- args[["ord.flip"]]
    args[["ord.flip"]] <- NULL
  }
  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if(is.null(vol$vol3D.data)){
    stop ("empty vol$vol3D.data.")
  }
  vol$dxyz <- c(dxy,1)
  vol$orientation <- c(1,0,0,0,1,0)
  vol$xyz.from.ijk[1,] <- c(dxy[1],0,0,pt00[1])
  vol$xyz.from.ijk[2,] <- c(0,dxy[2],0,pt00[2])
  vol$xyz.from.ijk[3,] <- c(0,0,1,0)
  vol$xyz0 <- (as.matrix(expand.grid(0,0,vol$k.idx,1)) %*% t(vol$xyz.from.ijk))[,1:3]
  
  expt.pt <-get.extreme.pt(vol)
  expt.pt_ <- expt.pt[1:2,] + vol$dxyz[c(1,2,1,2)]* c(-0.5,-0.5,0.5, 0.5)
  
  
  if (is.null(abs.rng)){
    abs.rng <- as.numeric(expt.pt_[1,])
    if (flip) abs.rng[1:2] <- abs.rng[2:1]
  } 
  if (is.null(ord.rng)){
    ord.rng <- as.numeric(expt.pt_[2,])
    if (flop) ord.rng[1:2] <- ord.rng[2:1]
  } 
  
  if (!add){
    if (is.null(main)) main = paste("@ k =",k)
    args[["x"]] <- abs.rng
    args[["y"]] <- ord.rng
    args[["type"]] <- "n"
    args[["xlim"]] <- abs.rng
    args[["ylim"]] <- ord.rng
    args[["main"]] <- main
    args[["xlab"]] <- abs.lab
    args[["ylab"]] <- ord.lab
    args[["xaxs"]] <- "i"
    args[["yaxs"]] <- "i"
    args[["asp"]] <- 1
    do.call(plot,args)
    
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
  }
  
  if (!(k %in% vol$k.idx)){
    text((par("usr")[1]+par("usr")[2])/2, (par("usr")[3]+par("usr")[4])/2, "missing plane",
         col= "red", cex = 2)
  } else {
    usr.x <- par("usr")[1:2]
    usr.y <- par("usr")[3:4]
    pt.min <- c(max(min(usr.x),min(abs.rng), expt.pt[1,1]),max(min(usr.y), min(ord.rng), expt.pt[2,1]), k )
    pt.max <- c(min(max(usr.x), max(abs.rng), expt.pt[1,2]), min(max(usr.y), max(ord.rng), expt.pt[2,2]), k )
    
    imin <- as.numeric(vol$min.pixel) #r[1]
    imax <- as.numeric(vol$max.pixel) #r[2]
    
    vol <- nesting.cube(vol, pt.min, pt.max)
    if (!is.na(vol$min.pixel) & !is.na(vol$max.pixel)) {
      expt.pt <- get.extreme.pt(vol)[1:2,] + vol$dxyz[c(1,2,1,2)]* c(-0.5,-0.5,0.5, 0.5)
      
      if (!is.null(breaks) & !any(is.na(breaks))) {
        if (length(breaks)!= length (col) + 1) stop("length(breaks) must be equal to length(col) + 1.")
      } else {
        breaks <- .pixel.scale (imin,imax,length(col))
      }
      
      if(!any(is.na(breaks))){
        if (!sat.transp) {
          breaks[1] <- imin - 1
          breaks[length (breaks)] <- imax + 1
        }
        map <- aperm (as.matrix(vol$vol3D.data[,,1], dim =vol$n.ijk[1:2] ), perm=c (2, 1))
        map.layer <- matrix (col[cut (as.numeric (map),breaks, include.lowest=TRUE)], nrow=nrow(map))
        
        if (usr.x[2]<usr.x[1]) expt.pt[1,] <- rev(expt.pt[1,])
        if (usr.y[2]<usr.y[1]) expt.pt[2,] <- rev(expt.pt[2,])
        if ((usr.x[2]-usr.x[1])* vol$dxyz[1] < 0) map.layer <- map.layer[,ncol(map.layer):1]
        if ((usr.y[2]-usr.y[1])* vol$dxyz[2] > 0) map.layer <- map.layer[nrow(map.layer):1,]
        rasterImage (map.layer, xleft = expt.pt[1,1], ybottom = expt.pt[2,1],
                     xright = expt.pt[1,2], ytop = expt.pt[2,2], interpolate = interpolate)
      }
      
    }
  }
  
}
