#' Display the 3D contours of the RoI
#' @description The \code{display.3D.contour} function performs a 3D display of the selected RoI in the chosen coordinate system.
#' @param struct "struct" class object. See \link[espadon]{espadon.class} for 
#' class definitions.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} object.
#' By default \code{roi.idx = NULL}. See Details.
#' @param roi.col Color of the RoI. If \code{roi.col = NULL} (default),the RoI 
#' colors are specified in the \code{struct$roi.info}.
#' @param roi.print Boolean vector indicating whether to display the pseudonym
#' of the RoI.
#' @param roi.cex Numeric character expansion factor of RoI name if 
#' \code{roi.print = TRUE}, defaults to 1.
#' @param roi.lwd Line width of the RoI, by default at 1. 
#' @param display.ref Pseudonym of frame of reference of the display.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} or
#' \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, \code{display.ref} must be equal to \code{NULL} or to
#' \code{struct$ref.pseudo}.
#' @param FoR.axis Boolean or numeric, by default set to \code{FALSE}. If \code{FoR.axis = TRUE},
#' the function displays 200 mm length director vectors of the frame of reference.
#' If \code{FoR.axis} is numeric, it represent the length in mm of the director vectors.
#' @param FoR.col Color of the frame of reference.
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all
#' \code{NULL}, then all of the RoI are selected.
#' @return If the concerned regions of interest (RoI) \code{struct} exist, 
#' it displays the 3D contours of these RoI in the current \pkg{RGL} window if it exists, 
#' in a new window otherwise.
#' @examples

#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtstruct", roi.name = "eye",
#'                              dxyz = rep (step, 3))
#' library (rgl)
#' open3d()
#' bg3d ("black")
#' display.3D.contour (struct = patient$rtstruct[[1]], roi.print = TRUE)

#' @importFrom rgl lines3d text3d cur3d points3d
#' @importFrom stats ecdf
#' @importFrom grDevices col2rgb rgb
#' @importFrom methods is
#' @export
display.3D.contour <- function (struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                                roi.col = NULL, roi.print = FALSE, roi.lwd=1, roi.cex=1,
                                display.ref = struct$ref.pseudo, T.MAT = NULL,
                                FoR.axis=FALSE, FoR.col="black") {
  
  if (!is (struct, "struct")){
    stop ("struct should be a struct class object.")
    # return (NULL)
  }
  M <- get.rigid.M (T.MAT, struct$ref.pseudo, display.ref)
  if (is.null (M)) stop ("cannot display anything in selected frame of reference.")
  
  M <- M %*% struct$ref.from.contour
  s.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  
 
  if (is.null(s.idx)){
    stop ("no RoI to display.")
    # return (NULL)
  }
  
  if (is.null(roi.col)) {
    col <- struct$roi.info$color[s.idx]
  } else {
    col <- rep(roi.col,length(s.idx))
  }
  for (sel in 1:length(s.idx)) { #struct$struct[[s.idx]]) {
    i <- s.idx[sel]
    if (!is.null(struct$roi.data[[i]])){
      for (j in 1:length (struct$roi.data[[i]])) {
        pt <- as.matrix(cbind(struct$roi.data[[i]][[j]]$pt,1)) %*% t(M)
        if (castup.str(struct$roi.data[[i]][[j]]$type) != "POINT") {lines3d (pt[ ,1], pt[ ,2], pt[ ,3], col=col[sel], lwd=roi.lwd)
        } else { points3d (pt[ ,1], pt[ ,2], pt[ ,3], col=col[sel], lwd=roi.lwd)}
      }
      
      pt <- as.numeric(c( struct$roi.info$Gx[i], struct$roi.info$Gy[i],struct$roi.info$Gz[i], 1) %*% t(M))
      
      if (roi.print) text3d (pt[1], pt[2], pt[3], struct$roi.info$roi.pseudo[i], 
                             col=col[sel], cex=roi.cex, offset = 0, adj=0.5)
    }
  }
  
  
  if (is.numeric(FoR.axis)) {
    
    lines3d (c(0, abs(FoR.axis)), c(0, 0), c(0, 0), col=FoR.col)
    text3d ((abs(FoR.axis) + 10), 0, 0, "X", col=FoR.col)
    
    lines3d (c(0, 0), c(0, abs(FoR.axis)), c(0, 0), col=FoR.col)
    text3d (0, (abs(FoR.axis) + 10), 0, "Y", col=FoR.col)
    
    lines3d (c(0, 0), c(0, 0), c(0, abs(FoR.axis)), col=FoR.col)
    text3d (0, 0, (abs(FoR.axis) + 10), "Z", col=FoR.col)
    
    
  } else if (is.logical(FoR.axis)){
    if (FoR.axis){
      lines3d (c(0, 200), c(0, 0), c(0, 0), col=FoR.col)
      text3d (210, 0, 0, "X", col=FoR.col)
      
      lines3d (c(0, 0), c(0, 200), c(0, 0), col=FoR.col)
      text3d (0, 210, 0, "Y", col=FoR.col)
      
      lines3d (c(0, 0), c(0, 0), c(0, 200), col=FoR.col)
      text3d (0, 0, 210, "Z", col=FoR.col)
    }
  }
}