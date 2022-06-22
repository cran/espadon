#' Display 3D sections of a patient
#' @description The \code{display.3D.sections} function displays transverse,
#' sagittal and frontal views at a point in 3D.
#' @param vol "volume" class object to display. See \link[espadon]{espadon.class} for 
#' class definitions.
#' @param cross.pt Vector of x, y, z coordinates, representing the cross point of the 3 planes.
#' @param display.ref Character string. Pseudonym of the frame of reference used for display.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} or
#' \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{vol} must be 
#' displayed in \code{display.ref = vol$ref.pseudo}.
#' @param col Vector, representing the color palette of the image. Transparent 
#' colors are not represented. 
#' @param breaks  One of :
#' \itemize{
#' \item \code{NULL} : the minimum and the maximum value of the \code{vol} define the range.
#' \item Vector giving the breakpoints of each color.
#' }
#' @param trans Boolean. If \code{TRUE} (default), the transverse view is displayed.
#' @param sagi Boolean. If \code{TRUE} (default), the sagittal view is displayed.
#' @param front Boolean. If \code{TRUE} (default), the frontal view is displayed.
#' @param border Boolean. If \code{TRUE} (default), the borders of the planes are displayed
#' @param border.col Color of planes borders
#' @return Returns a display of transverse, sagittal and frontal views of \code{vol} 
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
#'                     col= pal.RVV(200, alpha= c(rep(0,90), rep(1,110))),
#'                     breaks = seq(-1000, 1000, length.out = 201))


#' @export
#' @importFrom methods is
display.3D.sections <- function(vol, cross.pt = c (0, 0, 0), display.ref = vol$ref.pseudo, 
                                T.MAT = NULL,  
                                col = grey.colors(10, start = 0, end = 1, alpha = c(rep(0,1),rep(1,9))),
                                breaks = NULL, 
                                trans = TRUE,
                                sagi = TRUE,
                                front = TRUE,
                                border =TRUE, border.col= "#379DA2") {
  
  if (!is (vol, "volume")) 
    stop ("vol should be a volume class object.")
  
  if (is.null(vol$vol3D.data)) 
    stop  ("empty vol$vol3D.data.")
  
  
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo = display.ref, T.MAT= T.MAT)
  
  if (is.null (vol_)) stop ("cannot display anything in selected frame of reference.")
  
  if (trans) {
    pl.XY <- get.plane (vol_, cross.pt, c(1, 0, 0, 0, 1, 0), interpolate = TRUE)
    display.3D.stack(pl.XY, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                     breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                     line.col = border.col)
    }
  if (sagi) {
    pl.YZ <- get.plane (vol_, cross.pt, c(0, 0, 1, 0, 1, 0), interpolate = TRUE)  
    display.3D.stack(pl.YZ, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                     breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                     line.col = border.col)
    }
  if (front) {
    pl.ZX <- get.plane (vol_, cross.pt, c(1, 0, 0, 0, 0, 1), interpolate = TRUE) 
    display.3D.stack(pl.ZX, 0, display.ref = display.ref, T.MAT = T.MAT, col =col, 
                     breaks = breaks, cube=FALSE, border=border,ktext = FALSE,
                     line.col = border.col)
    } 

 
}