#' Display in 3D the selected planes of an \pkg{espadon} class volume
#' @description The \code{display.3D.stack} function displays in 3D the requested 
#' cutting planes of a "volume" class object.
#' @param vol "volume" class object to display.
#' @param k.idx vector of cutting plane numbers to be displayed, to be chosen in 
#' \code{vol$k.idx}. By default \code{k.idx} is a vector of 10 uniformly 
#' distributed cutting planes in the volume.
#' @param display.ref Character string. Pseudonym of the frame of reference used 
#' for display.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{vol} must 
#' be displayed in \code{display.ref = vol$ref.pseudo}.
#' @param col Vector, representing the color palette of the image. Transparent 
#' colors are not represented. 
#' @param breaks  One of :
#' \itemize{
#' \item \code{NULL} : The minimum and the maximum value of the \code{vol} define 
#' the range.
#' \item Vector giving the breakpoints of each color.
#' }
#' @param cube Boolean. If \code{TRUE} (default), the "volume" edges are displayed.
#' @param border Boolean. If \code{TRUE} (default), the borders of the planes defined 
#' in \code{k.idx} are displayed.
#' @param ktext Boolean. If \code{TRUE} (default), the  selected cutting plane numbers 
#' are displayed.
#' @param line.col Color of cube, planes and texts displayed.
#' @param line.lwd Line width of the border and cube, by default at 1. 
#' @param cex Numeric character expansion factor of displayed plan numbers.

#' @return  Returns a display of the \code{k.idx} cutting planes of \code{vol}, 
#' in the current \pkg{RGL} window if it exists, in a new window otherwise. The 
#' colors of the palettes are managed by \code{col} and \code{breaks}.  

#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "ct", dxyz = rep (step, 3))
#' 
#' # display o 3 planes
#' library (rgl)
#' open3d()
#' display.3D.stack (patient$ct[[1]],
#'                   col = pal.RVV (200, alpha = c(rep(0,90), rep (1, 110))))



#' @export
#' @importFrom rgl mesh3d wire3d segments3d texts3d material3d
#' @importFrom methods is
display.3D.stack <- function (vol,  
                             k.idx = unique(vol$k.idx[seq(1,vol$n.ijk[3], length.out=10)]), 
                             display.ref = vol$ref.pseudo, 
                             T.MAT = NULL,  
                             col = grey.colors(10, start = 0, end = 1, alpha = c(rep(0,1),rep(1,9))),
                             breaks = NULL,
                             cube = TRUE, 
                             border = TRUE,
                             ktext = TRUE,
                             line.col = "#379DA2",
                             line.lwd = 1,
                             cex = 1){
  
  if (!is (vol, "volume")) 
    stop ("vol should be a volume class object.")
  
  if (is.null(vol$vol3D.data)) 
    stop  ("empty vol$vol3D.data.")
  
  M <- get.rigid.M (T.MAT, vol$ref.pseudo, display.ref)
  if (is.null (M)) stop ("cannot display anything in selected frame of reference.")
  
  imin <- as.numeric(vol$min.pixel)
  imax <- as.numeric(vol$max.pixel) 
  
  if (!is.null(breaks) & !any(is.na(breaks))) {
    if (length(breaks)!= length (col) + 1) stop ("length(breaks) must be equal to length(col) + 1.")
    
  } else {
    breaks <- .pixel.scale (imin,imax,length(col))
  }
  if (breaks[1]>imin) breaks[1] <- imin - 1
  if (breaks[length (breaks)]<imax) breaks[length (breaks)] <- imax + 1
  
  color  <- substr(col,1,7)
  alpha <- substr(col,8,9)
  if (nchar(alpha[1])==0) alpha= rep ("ff",length(col))
  cube_disp <- vol$cube.idx
  cube_disp[c(1:3,6:7,11,13,15,17:18,22,29)] <- cube_disp[c(1:3,6:7,11,13,15,17:18,22,29)] - 0.5
  cube_disp[c(5,9:10,14,19,21,23,25:27,30:31)] <- cube_disp[c(5,9:10,14,19,21,23,25:27,30:31)] + 0.5
  
  if (cube & !any (apply(vol$cube.idx == 0 ,1, all))){

    line <-unique (do.call(rbind.data.frame,lapply(1:ncol(vol$cube.idx), function(i) 
      data.frame(a=rep(i,3),b=which(apply(sweep(vol$cube.idx,1,vol$cube.idx[,i])==0,2,sum)==3)))))

    pt <- (M %*% vol$xyz.from.ijk %*% cube_disp)[1:3,]
    for (i in 1:nrow(line)) rgl::segments3d (t (pt[ ,c (line[i,1], line[i,2])]), color = line.col, lwd=line.lwd)
  }

  for (k in k.idx){
    p <- cube_disp[,1:4]
    p[3,] <- k
    pt <- (M %*% vol$xyz.from.ijk %*% p)[1:3,]
    if (border){    
      segments3d (t (pt[,c(1,2)]), color=line.col, lwd=line.lwd)
      segments3d (t (pt[,c(2,3)]), color=line.col, lwd=line.lwd)
      segments3d (t (pt[,c(3,4)]), color=line.col, lwd=line.lwd)
      segments3d (t (pt[,c(4,1)]), color=line.col, lwd=line.lwd)
      
      if (ktext) texts3d(t (pt[,c(1)]), texts=paste0("k=",k," "), 
                         color=line.col,cex=cex, adj=1)  
    }
      
    idx <- which (k == vol$k.idx)
    map <-  array(vol$vol3D.data[ , ,idx],dim=vol$n.ijk[1:2])
    
    if(!all(is.na(map))){
      cut.map <- cut (as.numeric (map),breaks, include.lowest=TRUE)
      alpha.layer <-  matrix (alpha[cut.map], nrow=nrow(map))
      map.layer <- matrix (color[cut.map], nrow=nrow(map))
      ind <- which(alpha.layer!="00", arr.ind = T) 
      
      if(nrow(ind)>0){
        pt.ijk <- as.matrix(cbind(as.data.frame(ind)-1,k=vol$k.idx[idx],t=1))
        
        pt <- (as.matrix(cbind(as.data.frame(ind)-1,k=vol$k.idx[idx],t=1)) %*% t(M %*% vol$xyz.from.ijk))
        A <- (as.matrix (sweep(pt.ijk,2, c ( 0.5,  0.5, 0, 0))) %*% t(M %*% vol$xyz.from.ijk))
        B <- (as.matrix (sweep(pt.ijk,2, c (-0.5,  0.5, 0, 0))) %*% t(M %*% vol$xyz.from.ijk))
        C <- (as.matrix (sweep(pt.ijk,2, c (-0.5, -0.5, 0, 0))) %*% t(M %*% vol$xyz.from.ijk))
        D <- (as.matrix (sweep(pt.ijk,2, c ( 0.5, -0.5, 0, 0))) %*% t(M %*% vol$xyz.from.ijk))
        
        mesh <- list()
        mesh$vb <- t((rbind(A,B,C,D)))
        mesh$qu <- rbind(1:nrow(A),(1:nrow(A))+nrow(A),(1:nrow(A))+2*nrow(A),(1:nrow(A))+3*nrow(A))
        material <-material3d()
        material$specular = "black"
        material$color <- as.character(map.layer[ind])
        material$lit <- FALSE
        m <- mesh3d(vertices=mesh$vb, quads=mesh$qu , meshColor="faces", material=material)
        
        shade3d(m)
      }
    }
  }

}
