#' Display of a 2D histogram
#' @description The \code{display.2D.histo} function displays the density map of
#' a "histo2D" class object.
#' @param histo.2D "histo2D" class object.
#' @param add Boolean indicating whether to display the background image.
#' @param main Title of the background image. If \code{main = NULL},
#' the title indicates "2D histogram".
#' @param x.lab Label of the x-axis of the background image. If \code{x.lab = NULL},
#' this label is \code{histo.2D$x.file.src}
#' @param y.lab Label of the y-axis of the background image. If \code{y.lab = NULL},
#' this label is \code{histo.2D$y.file.src}.
#' @param x.lim Vector, representing the range of the x-axis.
#' @param y.lim Vector, representing the range of the y-axis.
#' @param bg Background color of the image. By default, this color is black.
#' @param i.rng Vector of 2 elements giving the minimum and the maximum intensity
#' of the image. If \code{i.rng = NULL}, then the minimum is 0 and the maximum the
#' maximum density.
#' @param display.mode function display mode. See Details.
#' @param col Color of the color gradient in \code{display.mode = "mono.color"}, 
#' or of the level lines in \code{display.mode = "line"}. By default, this color 
#' is white.
#' @param alpha A number from 0 to 1, indicating the opacity of the image in
#' \code{"rainbow.color"} mode. Default alpha = 1.
#' @param line.pc.levels Vector of level lines in percent of maximum density
#' in \code{display.mode = "line"}. By default lines 1% and 100% are displayed.
#' @param line.lwd Line thickness of the level lines in \code{display.mode = "line"}.
#' @param line.lty Type of lines for level lines in \code{display.mode = "line"}.
#' @details The \code{display.mode} argument can be set to three values: 
#' \code{"mono.color"}, \code{"rainbow.color"}, or \code{"line"}. The 2D 
#' histogram graph is displayed by default in \code{"mono.color"} mode.
#' \itemize{
#' \item The \code{"mono.color"} mode displays a gradient of the color defined
#' by the col argument, depending on the intensity of \code{$density.map} 2-dimensional
#' array.
#' \item The \code{"rainbow.color"} mode makes a display according to the
#' \code{"rainbow"} palette, while managing the opacity of the colors.
#' \item The \code{"line"} mode draws level lines defined in percent by the
#' \code{line.pc.levels} argument.
#' }
#' @return Returns a display of the density map of \code{histo.2D}. This one must 
#' be an object of class "histo2D". See \link[espadon]{espadon.class} for 
#' class definitions.
#' @seealso \link[espadon]{histo.2D}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "mr", "rtstruct"), 
#'                              roi.name =  "brain", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' MR <- patient$mr[[1]]
#' S <- patient$rtstruct[[1]]
#' T.MAT <- patient$T.MAT
#'
#' # restriction of the volume around the RoI
#' CT.on.roi <- nesting.roi (CT, S, roi.name = "brain", vol.restrict = TRUE,
#'                           xyz.margin = c (1, 1, 1), alias = CT$description)
#' MR.on.CT <- vol.regrid (vol = MR, back.vol = CT.on.roi, interpolate = TRUE,
#'                         T.MAT = T.MAT, alias = CT$description,
#'                         description = NULL)

#' # selection of voxels included in the RoI.
#' roi.bin <- bin.from.roi (vol = CT.on.roi, struct = S, roi.sname = "brain")
#' MR.select <- vol.from.bin (MR.on.CT, roi.bin, alias = MR$description)
#' CT.select <- vol.from.bin (CT.on.roi, roi.bin, alias = CT$description)

#' # 2D histogram
#' H2D <- histo.2D (MR.select, CT.select, x.breaks = seq (50, 400, 10),
#' 			  y.breaks = seq (50, 400, 10), alias = "H2D MR1 MR2")
#' display.2D.histo (H2D, display.mode = "mono.color", col = "#ffff00", 
#'                   main ="mono color mode")
#' display.2D.histo (H2D, display.mode = "rainbow.color", main ="rainbow mode")
#' display.2D.histo (H2D, display.mode = "line", main ="level lines mode",
#'                   line.pc.levels = c (0, 25, 50, 75, 100), col = "#ff0000")

#' @importFrom stats ecdf
#' @importFrom grDevices contourLines rainbow
#' @importFrom colorspace hex2RGB
#' @importFrom methods is
#' @export
display.2D.histo  <- function (histo.2D, add=TRUE, main= NULL, 
                               x.lab = NULL, y.lab = NULL, 
                               x.lim=NULL, y.lim = NULL,
                               bg="#000000", i.rng = NULL, 
                               display.mode=c("mono.color","rainbow.color","line"),
                               col= "#ffffff", alpha = 1,
                               line.pc.levels = c (1, 100),
                               line.lwd = 2, line.lty = 1){
  
  
  .apply.color.map <- function (M, alpha) {
    p <- rev(rainbow(256, s = seq(256, 1, by=-1)/256, v = 1, start = 0, end = 4/6, alpha = alpha * seq (256, 1, by=-1)/256))
    p[M * 255 + 1]
  }
  
  if (!is (histo.2D, "histo2D")) stop ("histo2D must be a histo2D class object.")

  

  scale1 <- histo.2D$x.mids[c(1,length(histo.2D$x.mids))]
  
  scale2 <- histo.2D$y.mids[c(1,length(histo.2D$y.mids))]
  
  if(is.null(main)) main <-"2D histogram"
  if(is.null(x.lab)) x.lab <- histo.2D$x.file.src
  if(is.null(y.lab)) y.lab <- histo.2D$y.file.src
  abs.left <- histo.2D$x.breaks[1]
  ord.bottom <- histo.2D$y.breaks[1]
  abs.right <- histo.2D$x.breaks[length(histo.2D$x.breaks)]
  ord.top <- histo.2D$y.breaks[length(histo.2D$y.breaks)]
 
  xlim <- suppressWarnings(range (x.lim, na.rm=T))
  ylim <- suppressWarnings(range (y.lim, na.rm=T))
  if (any(is.infinite(xlim)))  xlim <- c(abs.left,abs.right)
  if (any(is.infinite(ylim)))  ylim <- c(ord.bottom,ord.top)
  if (add){
    plot (scale1, scale2, type="n", xlim = xlim, ylim = ylim, #asp = 1,
          main=main ,xlab=x.lab, ylab = y.lab, xaxs="i", yaxs="i")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
  }
  
  
  if (!is.null(i.rng)) {
    imin <- i.rng[1]
    imax <- i.rng[2]
  } else {
    imin <- 0
    imax <- max(histo.2D$density.map, na.rm=TRUE)
  }
  
  
  map.layer <- (t(histo.2D$density.map) -imin)/(imax - imin)
  map.layer[map.layer>1] <- 1
  map.layer[map.layer<0] <- 0
  
  display.mode  <- display.mode[1]
  if (display.mode[1]=="rainbow.color"){
    
    map.layer <- matrix (.apply.color.map(map.layer, alpha), nrow=nrow(map.layer))
    rasterImage (map.layer[nrow(map.layer):1,], xleft = abs.left,
                 ybottom = ord.bottom, xright = abs.right, ytop = ord.top,
                 interpolate=FALSE)
  } else if (display.mode[1]=="line") {
    
    # z <- histo.2D
    imax <- max(histo.2D$density.map, na.rm=TRUE)
    z <-(histo.2D$density.map/imax) *100
    z[is.na(z)]<- 0
    ctr <- contourLines(x = histo.2D$x.mids, y= histo.2D$y.mids ,
                        z=z, levels = line.pc.levels)
    if (length(ctr)>0) for (i in 1:length(ctr)) {
      lines(c (ctr[[i]]$x,ctr[[i]]$x[1]), 
            c (ctr[[i]]$y,ctr[[i]]$y[1]), col=col, lwd=line.lwd, lty=line.lty)
    }
  } else {
    
    
    ramp<-seq(0,1,1/256)
    full.scale <- hex2RGB(col[1])@coords - hex2RGB(bg[1])@coords
    col.pal <- sapply(ramp, function(k) round((full.scale*k + hex2RGB(bg[1])@coords)*255))
    
    pal <-sapply(1:ncol(col.pal), function (c.idx) {
      vect <- as.character (as.hexmode(col.pal[,c.idx]))
      paste ("#",paste(sapply(vect, function(e)
        ifelse(nchar(e)==1,paste("0",e,sep=""),e)),collapse=""),sep="")
    })
    map.layer <- matrix (pal[round(map.layer*255 +1)], nrow=nrow(map.layer))
    rasterImage (map.layer[nrow(map.layer):1,], xleft = abs.left,
                 ybottom = ord.bottom, xright = abs.right, ytop = ord.top,
                 interpolate=FALSE)
    
  }
  
  
  
}