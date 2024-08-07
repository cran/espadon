#' Display 3D sections of a patient
#' @description The \code{display.3D.sections} function displays transverse,
#' sagittal and frontal views at a point in 3D.
#' @param obj "volume" class object to display. See \link[espadon]{espadon.class} for 
#' class definitions.
#' @param cross.pt Vector of x, y, z coordinates, representing the cross point of the 3 planes in display.plane.
#' @param display.ref Character string. Pseudonym of the frame of reference used for display.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} or
#' \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{obj} must be 
#' displayed in \code{display.ref = obj$ref.pseudo}.
#' @param col Vector, representing the color palette of the image. Transparent 
#' colors are not represented. 
#' @param breaks  One of :
#' \itemize{
#' \item \code{NULL} : the minimum and the maximum value of the \code{obj} define the range.
#' \item Vector giving the breakpoints of each color.
#' }
#' @param trans Boolean. If \code{TRUE} (default), the transverse view is displayed.
#' @param sagi Boolean. If \code{TRUE} (default), the sagittal view is displayed.
#' @param front Boolean. If \code{TRUE} (default), the frontal view is displayed.
#' @param border Boolean. If \code{TRUE} (default), the borders of the planes are displayed
#' @param border.col Color of planes borders
#' @param ... Argument for deprecated arguments
#' @return Returns a display of transverse, sagittal and frontal views of \code{obj} 
#' at \code{cross.pt}  in the current \pkg{RGL} window if it exists, in a new 
#' window otherwise. Palette colors are managed by \code{col} and \code{breaks}.  
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "ct", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' 
#' library (rgl)
#' open3d()
#' display.3D.sections(CT, cross.pt= c(0, 50, 80),
#'                     col= pal.RVV(200, alpha= c(rep(0,90), rep(1,110))))


#' @export
#' @importFrom methods is
display.3D.sections <- function(obj, cross.pt = c (0, 0, 0), display.ref = obj$ref.pseudo, 
                                T.MAT = NULL,  
                                col = grey.colors(10, start = 0, end = 1, alpha = c(rep(0,1),rep(1,9))),
                                breaks = NULL, 
                                trans = TRUE,
                                sagi = TRUE,
                                front = TRUE,
                                border =TRUE, border.col= "#379DA2",...) {
  passed <- names(as.list(match.call())[-1])
  args <- list(...)
  if (!("obj" %in% passed)){
    if (is.null(args[['vol']])) stop('argument "obj" is missing, with no default')
    obj <- args[['vol']]
  }
  
  if (!is (obj, "volume")) stop ("obj should be a volume class object.")
  
  ept <-get.extreme.pt(obj, ref.pseudo = display.ref, T.MAT= T.MAT)
  if (any(sign(c(cross.pt-ept[,1], ept[,2]-cross.pt))==-1)) 
    stop("cross.pt is not included in the obj")
  
  
  if (is (obj, "volume")) {
    obj_ <- vol.in.new.ref(obj, new.ref.pseudo = display.ref, T.MAT= T.MAT)
    
    if (is.null (obj_)) stop ("cannot display anything in selected frame of reference.")
    
    if (is.null(obj$vol3D.data)) stop  ("empty obj$vol3D.data.")
    
    imin <- as.numeric(obj$min.pixel) #r[1]
    imax <- as.numeric(obj$max.pixel) #r[2]
    if (!is.null(breaks) & !any(is.na(breaks))) {
      if (length(breaks)!= length (col) + 1) stop("length(breaks) must be equal to length(col) + 1.")
      
    } else {
      breaks <- .pixel.scale (imin,imax,length(col))
    }
    if (breaks[1]>imin) breaks[1] <- imin - 1
    if (breaks[length (breaks)]<imax) breaks[length (breaks)] <- imax + 1
    
    
    
    
    if (trans) {
      pl.XY <- get.plane (obj_, cross.pt, c(1, 0, 0, 0, 1, 0), interpolate = TRUE)
      if (!is.null(pl.XY))
        display.3D.stack(pl.XY, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                         breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                         line.col = border.col)
    }
    if (sagi) {
      pl.YZ <- get.plane (obj_, cross.pt, c(0, 0, 1, 0, 1, 0), interpolate = TRUE) 
      if (!is.null(pl.YZ))
        display.3D.stack(pl.YZ, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                         breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                         line.col = border.col)
    }
    if (front) {
      pl.ZX <- get.plane (obj_, cross.pt, c(1, 0, 0, 0, 0, 1), interpolate = TRUE) 
      if (!is.null(pl.ZX))
        display.3D.stack(pl.ZX, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                         breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                         line.col = border.col)
    } 
  }
 
}