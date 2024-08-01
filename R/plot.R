#' @name plot
#' @aliases plot.volume
#' @aliases plot.struct
#' @aliases plot.mesh
#' @title plot a 2D cut of a 3D object
#' @description The \code{plot} function displays the requested map
#' of espadon objects of class "volume", "struct", "mesh".

#' @param x object of class "volume", "struct" or "mesh". See \link[espadon]{espadon.class} 
#' @param ... others parameters of plot functions. See details 
#' @param view.type character string among the values 'ij', 'ji', 'ik', 'ki', 
#' 'jk', 'kj', 'yx', 'xz', 'zx', 'yz', 'zy', 'trans', 'front' or 'sagi" 
#' representing the map to be displayed.
#' @param view.coord value representing the coordinate where the map is displayed.
#' This parameter can be a 3D-vector, representing the coordinate of the point on 
#' the displayed map. If \code{NULL}, the display is located in the center of the object.

#' @param flip Boolean defaults to \code{FALSE} flipping the horizontal axis 
#' of the background image.
#' @param flop Boolean defaults to \code{FALSE} flipping the vertical axis 
#' of the background image.


#' @return Returns a display of the  \mjeqn{k^{th}}{ascii} image plane of \code{x}.
#' @seealso \link[espadon]{display.plane}, \link[espadon]{display.kplane}, 
#' \link[espadon]{display.palette}, \link[espadon]{pal.RVV}, \link[espadon]{pal.rainbow}.
#' @importFrom graphics layout
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct","mr", "rtdose", "rtstruct"),
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' MR <- patient$mr[[1]]
#' CT <- patient$ct[[1]]
#' D <- patient$rtdose[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' # display 1
#' layout (matrix(c(1,1,2,3), ncol=2), widths=c(1,0.2))
#' plot (CT, view.coord = 0, col = pal.RVV(255))
#' S_plot <- plot (S, view.coord = 0, add = TRUE, lwd = 2)
#' display.palette (col = pal.RVV(255), main="HU")
#' display.legend (S_plot, bg="white", text.col="black", lwd = 2, cex = 1.1)

#' # display 2
#' layout (matrix(c (1, 1, 2, 3), ncol = 2), widths = c (1, 0.2))
#' # Coordinates of the PTV barycenter in CT frame of reference
#' G <- as.numeric (S$roi.info[S$roi.info$roi.pseudo == "ptv", 
#'                             c ("Gx", "Gy", "Gz")])
#' # Coordinates of the PTV baricenter in MR frame of reference
#' G_MR <- as.numeric (c(G, 1) %*% 
#'                    t(get.rigid.M(CT$ref.pseudo, MR$ref.pseudo, 
#'                                T.MAT = patient$T.MAT)))[1:3]
#' plot (MR, view.type = "sagi", view.coord =  G_MR,
#'       col = grey.colors (255, start = 0, end = 1),
#'       breaks = seq (0, 500, length.out = 256) , bg = "darkblue")
#' plot (vol.in.new.ref(D, MR$ref.pseudo, T.MAT = patient$T.MAT),
#'       view.type = "sagi", view.coord = G_MR,
#'       col = pal.rainbow(255), add = TRUE)
#' display.palette (col = grey.colors (255, start = 0, end = 1),
#'                  breaks = seq (0, 500, length.out = 256), main="MR")
#' display.palette (col =  pal.rainbow(255),
#'                  breaks = seq (D$min.pixel, D$max.pixel, length.out = 256),
#'                  main="Gy")
#'layout(1)

#################################################################################
#' @rdname plot
#' @param cut.interpolate Boolean, indicating whether to calculate the volume cut 
#' using linear interpolation.
#' @param display.interpolate Boolean, indicating whether to apply linear 
#' interpolation when displaying the cut.
#' @param col  Vector, representing the color palette of the image, if \code{x} 
#' is of class 'volume'. Color of the mesh outline if object \code{x} is of class 'mesh'.
#' @param breaks One of :
#' \itemize{
#' \item \code{NULL} : the minimum and the maximum value of the object \code{x} define 
#' the range.
#' \item Vector giving the breakpoints of each color. Outside values are transparent,
#' leaving the background visible, depending on \code{sat.transp}.
#' }
#' @param sat.transp Boolean. If \code{TRUE}, outside values are transparent, else 
#' set to \code{breaks} limits colors.
#' @details ... can be xlim, ylim, add, bg etc. If \code{view.type} is egal to 
#' 'trans' or 'front' or 'sagi', the direction of xlim and ylim is ignored.

#' @importFrom grDevices rainbow grey.colors
#' @export
plot.volume <- function (x, ..., view.type = "trans", 
                         view.coord = NULL, 
                         flip = FALSE, flop = FALSE, 
                         cut.interpolate=TRUE,
                         display.interpolate=FALSE,
                         col = grey.colors (255, start = 0, end = 1),
                         breaks = NULL, sat.transp = FALSE) {
  
  args <- tryCatch(list(...), error = function(e)list())
  add <- args[["add"]]; if (is.null(add)) add <- FALSE; args[["add"]]<-NULL
  bg <- args[["bg"]]; if (is.null(bg)) bg <- "#000000";args[["bg"]] <- NULL
  
  if(is.null(x$vol3D.data)){
    stop ("empty x$vol3D.data.")
  }
  
  coord.map <-c('ij', 'ji', 'ik', 'ki', 'jk','kj' ,'xy', 'yx', 'xz', 'zx', 'yz', 'zy', "trans","front","sagi")
  m <- match(view.type[1],coord.map)
  if (is.na(m)) stop("Undefined abcisse or ordinate. 
  view.type must be set to 'ij' 'ji' 'ik' 'ki' 'jk' 'kj'  or 'xy' 'yx' 'xz' 'zx' 'yz' 'zy' 'trans' 'front' or 'sagi")
  
  
  
  abs.rng <- args[['xlim']]
  ord.rng <- args[['ylim']]
  
  if (m<7 ){
    if (is.null(view.coord)) view.coord <-  c(floor((x$n.ijk[1:2]-1)/2), 
                                              x$k.idx[floor(x$n.ijk[3]/2)])
    coord.label <- c("i","j","k")
    coord <- match(unique(unlist(strsplit(view.type,""))), coord.label)
    x$dxyz <- c(1,1,1)
    x$xyz.from.ijk <- diag(4)
    x$orientation <-  c(1, 0, 0, 0, 1, 0)
    x$xyz0 <- matrix(c(rep(0, 2*x$n.ijk[3]), x$k.idx), ncol=3)
    postfix <- ""
    if (!is.null(abs.rng)) {flip <- abs.rng[1]> abs.rng[2]}
    if (!is.null(ord.rng)) {flop <- ord.rng[1]> ord.rng[2]}
    
  }  else {
    if (is.null(view.coord))
      view.coord <- (x$xyz.from.ijk %*% 
                       c (floor (get.ijk.from.xyz (apply (get.extreme.pt (x), 1, mean),x)), 1))[1:3, 1]
    if (view.type=="trans") {
      view.type <- "xy"
      flop <- TRUE
      flip <- FALSE
    } else if (view.type=="front") {
      flip <- FALSE
      flop <- FALSE
      view.type <- "xz"
    } else if (view.type=="sagi")  {
      view.type <- "zy" 
      flip <- FALSE
      flop <- TRUE
    }else{
      if (!is.null(abs.rng)) {flip <- abs.rng[1]> abs.rng[2]}
      if (!is.null(ord.rng)) {flop <- ord.rng[1]> ord.rng[2]}
    }
    
    coord.label <- c("x","y","z")
    coord <- match(unique(unlist(strsplit(view.type,""))), coord.label)
    postfix <- " mm"
  }
  if (!is.null(ord.rng)) ord.rng <- range(ord.rng)
  if (!is.null(abs.rng)) abs.rng <- range(abs.rng)
  
  coord_ <- coord[!is.na(coord)]
  cut.lab <- coord.label[-coord_]
  cut.idx <- (1:3)[-coord_]
  
  # view.coord <- args[["view.coord"]]; 
  # if (is.null(view.coord)) view.coord <- mean(as.numeric(get.extreme.pt(x)[cut.idx, ])); 
  if (is.na(match(length(view.coord),c(1,3)))) stop("view.coord must be of length 1 or 3")
  # args[["view.coord"]]<-NULL
  
  if (!(length(view.coord) ==1 | length(view.coord)==3)) stop ("view.coord must be of length 1 or 3")
  if (length(view.coord) == 1) {
    origin = x$xyz0[1,]
    if (view.type=="ij" | view.type=="ji" | view.type=="xy" | view.type=="yx" | view.type == "trans") origin[3] <- view.coord
    if (view.type=="ik" | view.type=="ki" | view.type=="xz" | view.type=="zx" | view.type == "front") origin[2] <- view.coord
    if (view.type=="jk" | view.type=="kj" | view.type=="yz" | view.type=="zy" | view.type == "sagi") origin[1] <- view.coord
  } else {origin=view.coord; origin[coord] = x$xyz0[1,coord]}
  
  cut.pt <- origin[-coord_]
  coord_ <- coord.label[coord_]
  orientation <- switch (view.type,
                         "xy" = c(1,0,0,0,1,0),
                         "ij" = c(1,0,0,0,1,0),
                         "yx" = c(0,1,0,1,0,0),
                         "ji" = c(0,1,0,1,0,0),
                         "xz" = c(1,0,0,0,0,1),
                         "ik" = c(1,0,0,0,0,1),
                         "zx" = c(0,0,1,1,0,0),
                         "ki" = c(0,0,1,1,0,0),
                         "yz" = c(0,1,0,0,0,1),
                         "jk" = c(0,1,0,0,0,1),
                         "zy" = c(0,0,1,0,1,0),
                         "kj" = c(0,0,1,0,1,0)
  )
  
  
  vol_ <- get.plane(x, origin = origin, plane.orientation=orientation, interpolate = cut.interpolate)
  vol_$local.gridx <- NULL
  vol_$local.gridy <- NULL
  
  if (is.na(vol_$max.pixel)){
    message("No view @", cut.lab," = ", round(cut.pt,2),postfix)
    return(invisible(vol_))
    }
  t.mat <- ref.cutplane.add(vol_,ref.cutplane = "display.ref")
  vol_ <- vol.in.new.ref(vol_,"display.ref",t.mat)
  ext.pt <- get.extreme.pt(vol_)
  
  
  if (!is.null(abs.rng)) ext.pt[1,] <- range(abs.rng) #+ vol_$dxyz[1]*c(0.5,-0.5)
  if (!is.null(ord.rng)) ext.pt[2,] <- range(ord.rng) #+ vol_$dxyz[2]*c(0.5,-0.5)
  if (!is.null(abs.rng) |!is.null(ord.rng)) {
    vol_ <- nesting.cube(vol_, pt.min = as.numeric(ext.pt[,1]), pt.max = as.numeric(ext.pt[,2]))
    ext.pt <- get.extreme.pt(vol_)
  }
  
  if (!add){  
    if (is.null(args[['main']]))
      args[['main']] <- paste0(x$modality," (",x$description, ") @", cut.lab," = ", round(cut.pt,2),postfix)
    args[['type']] <- "n"
    args[['x']] <- as.numeric(ext.pt[1,1:2])
    args[['y']] <- as.numeric(ext.pt[2,1:2])
    if (is.null(args[['xaxs']])) args[['xaxs']] <- "i"
    if (is.null(args[['yaxs']])) args[['yaxs']] <- "i"
    if (is.null(args[['xlab']])) args[['xlab']] <- coord_[1]
    if (is.null(args[['ylab']])) args[['ylab']] <- coord_[2]
    args[['xlim']] <- as.numeric(ext.pt[1,1:2])+ abs(vol_$dxyz[1])*c(-0.5,+0.5)
    args[['ylim']] <- as.numeric(ext.pt[2,1:2])+ abs(vol_$dxyz[2])*c(-0.5,+0.5)
    
    if (flip) args[['xlim']] <- rev( args[['xlim']])
    if (flop) args[['ylim']] <- rev( args[['ylim']])
    
    args[['asp']] <- 1
    do.call(plot,args)
    abs.usr <- par("usr")[1:2]
    ord.usr <- par("usr")[3:4]
    rect(abs.usr[1], ord.usr[1], abs.usr[2], ord.usr[2], col = bg)
    ext.pt <- ext.pt[1:2,] + c(-0.5,-0.5,0.5,0.5)* vol_$dxyz[c(1:2, 1:2)]
  }else {
    abs.usr <- par("usr")[1:2]
    ord.usr <- par("usr")[3:4]
    pt.min <- c (max(ext.pt[1,1], min(abs.usr)), max(ext.pt[2,1], min(ord.usr)), ext.pt[3,1])
    pt.max <- c (min(ext.pt[1,2], max(abs.usr)), min(ext.pt[2,2], max(ord.usr)), ext.pt[3,2])
    
    vol_ <- nesting.cube(vol_, pt.min = pt.min, pt.max = pt.max)
    ext.pt <- get.extreme.pt(vol_)[1:2,] + c(-0.5,-0.5,0.5,0.5)* vol_$dxyz[c(1:2, 1:2)]
    
  }
  
  
  if (!is.na(vol_$max.pixel)){
    
    if (abs.usr[1]<abs.usr[2]){
      abs.left <- ext.pt[1,1] 
      abs.right <- ext.pt[1,2] 
    } else {
      abs.left <- ext.pt[1,2] 
      abs.right <- ext.pt[1,1] 
    }
    if (ord.usr[1]<ord.usr[2]){
      ord.bottom <- ext.pt[2,1] 
      ord.top <- ext.pt[2,2] 
    } else {
      ord.bottom <- ext.pt[2,2] 
      ord.top <- ext.pt[2,1] 
    }
    
    imin <- as.numeric(x$min.pixel) #r[1]
    imax <- as.numeric(x$max.pixel) #r[2]
    
    if (!is.null(breaks) & !any(is.na(breaks))) {
      if (length(breaks)!= length (col) + 1) stop("length(breaks) must be equal to length(col) + 1.")
      
    } else {
      pal <- attributes(col)$label
      if (!is.null(pal)) if (pal=="RVV"){imin=-1000; imax = 1000}
      breaks <- .pixel.scale (imin,imax,length(col))
    }
    map <- aperm (vol_$vol3D.data[,,1], perm=c (2, 1))
    
    if(!any(is.na(breaks))){
      if (!sat.transp) {
        # map[which(map <= imin)] <- imin
        # map[which(map >= imax)] <- imax
        breaks[1] <- imin - 1
        breaks[length (breaks)] <- imax + 1
      }
    }
    map.layer <- matrix (col[cut (as.numeric (map),breaks, include.lowest=TRUE)], nrow=nrow(map))
    
    if ((ord.usr[1] - ord.usr[2]) * vol_$dxyz[2] < 0) {
      map.layer <- map.layer[nrow(map.layer):1,]
    }
    
    if ((abs.usr[2] - abs.usr[1]) * vol_$dxyz[1] < 0) {
      map.layer <- map.layer[,ncol(map.layer):1]
    }
    
    rasterImage (map.layer, xleft = abs.left, ybottom = ord.bottom,
                 xright = abs.right, ytop = ord.top, interpolate = display.interpolate)
    
  } else{
    message("No view @", cut.lab," = ", round(cut.pt,2),postfix)
  }
  

  return(invisible(vol.in.new.ref(vol_,x$ref.pseudo,t.mat)))
}

#################################################################################

.suppress.roi <- function(x, kept.flag){
  x$nb.of.roi <- sum(kept.flag)
  x$roi.info <- x$roi.info[kept.flag, ,drop = FALSE]
  rownames( x$roi.info) <- NULL
  x$roi.obs <- x$roi.obs[kept.flag, ,drop = FALSE]
  rownames( x$roi.obs) <- NULL
  x$roi.data <- x$roi.data[kept.flag]
  return(x)
}
.range.intersection <- function(v1,v2,egal=c(TRUE,FALSE)){
  v1 <- range(v1)
  v2 <- range(v2)
  v <- c(v1[((v1>v2[1])| (egal[1] & v1==v2[1])) & (v1<v2[2] | (egal[2] & v1==v2[2]))],
         v2[((v2>v1[1])| (egal[1] & v2==v1[1])) & (v2<v1[2] | (egal[2] & v2==v1[2]))])
  if (length(v)==0) return(c(NA,NA))
  range(v)
}


#' @rdname plot
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} 
#' object. By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @param back.dxyz 3D vector of voxel size, used to calculate contours in frontal 
#' or sagittal view.

#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are
#' all set to \code{NULL}, all closed planar or point RoI, present at \code{view.coord}
#' are selected.
#' @export
plot.struct <- function(x, ..., view.type = "trans", 
                        view.coord =  NULL, 
                        flip = FALSE, flop = FALSE, 
                        roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                        back.dxyz =  c(0.5,0.5,x$thickness)) {
  
  args <- tryCatch(list(...), error = function(e)list())

  add <- args[["add"]]; if (is.null(add)) add <- FALSE; args[["add"]]<-NULL
  bg <- args[["bg"]]; if (is.null(bg)) bg <- "#000000";args[["bg"]] <- NULL
    
  coord.map <-c('xy', 'yx', 'xz', 'zx', 'yz', 'zy', "trans","front","sagi")
  m <- match(view.type[1],coord.map)
  if (is.na(m)) stop("Undefined abcisse or ordinate. view.type must be set to ", paste(coord.map, collapse=", "))
  
  x <- .suppress.roi(x, !is.na(x$roi.info[,"min.x"]))
  list.roi.idx <- select.names (x$roi.info$roi.pseudo, 
                                roi.name =roi.name, roi.sname = roi.sname, 
                                roi.idx = roi.idx)
    # args[["roi.name"]] <- args[["roi.sname"]] <- args[["roi.idx"]] <- NULL
  
  if (is.null(list.roi.idx)){
    message("no ROI to display")
    return(x)
  } 
  if (is.null(x$roi.data)) stop("empty x$roi.data")
  if (is.null(view.coord)) 
    view.coord <-  apply(get.extreme.pt(x, roi.idx = list.roi.idx),1,mean)
  
  dxyz <- c(back.dxyz, 0)
  cut.dxyz <- dxyz %*% t(x$ref.from.contour)
  
  abs.rng <- args[['xlim']]
  ord.rng <- args[['ylim']]
  
  if (view.type=="trans") {
    view.type <- "xy"
    flop <- TRUE
    flip <- FALSE
  } else if (view.type=="front") {
    flip <- FALSE
    flop <- FALSE
    view.type <- "xz"
  } else if (view.type=="sagi")  {
    view.type <- "zy" 
    flip <- FALSE
    flop <- TRUE
  }else{
    if (!is.null(abs.rng)) {flip <- abs.rng[1]> abs.rng[2]}
    if (!is.null(ord.rng)) {flop <- ord.rng[1]> ord.rng[2]}
  }
  if (!is.null(ord.rng)) ord.rng <- range(ord.rng)
  if (!is.null(abs.rng)) abs.rng <- range(abs.rng)
  

  coord.label <- c("x","y","z")
  coord <- match(unique(unlist(strsplit(view.type,""))), coord.label)
  postfix <- " mm"
  
  
  coord_ <- coord[!is.na(coord)]
  cut.lab <- coord.label[-coord_]
  cut.idx <- (1:3)[-coord_]
  
  # view.coord <- args[["view.coord"]]; 
  # if (is.null(view.coord)) view.coord <- mean(as.numeric(get.extreme.pt(x)[cut.idx, ])); 
  if (is.na(match(length(view.coord),c(1,3)))) stop("view.coord must be of length 1 or 3")
  # args[["view.coord"]]<-NULL
  
  if (length(view.coord) == 1) {
    if (view.type=="xy" | view.type=="yx" | view.type == "trans") origin <- c(0,0,view.coord)
    if (view.type=="xz" | view.type=="zx" | view.type == "front") origin <- c(0,view.coord,0)
    if (view.type=="yz" | view.type=="zy" | view.type == "sagi") origin <-c(view.coord,0,0)
  } else {origin=view.coord; origin[coord] = 0}
  
  cut.pt <- round(origin[-coord_],6)
  coord_ <- coord.label[coord_]
  orientation <- switch (view.type,
                         "xy" = c(1,0,0,0,1,0),
                         "yx" = c(0,1,0,1,0,0),
                         "xz" = c(1,0,0,0,0,1),
                         "zx" = c(0,0,1,1,0,0),
                         "yz" = c(0,1,0,0,0,1),
                         "zy" = c(0,0,1,0,1,0))
  

  col.minx <- which(colnames(x$roi.info)=="min.x")-1
  x$roi.info [,col.minx+1:10] <-round( x$roi.info [,col.minx+1:10],6)
  f.roi <- ((x$roi.info[list.roi.idx,col.minx+2*cut.idx-1] - cut.dxyz[cut.idx]/2 <=cut.pt) & 
              (x$roi.info[list.roi.idx,col.minx+2*cut.idx] + cut.dxyz[cut.idx]/2 > cut.pt) ) | 
    x$roi.info[list.roi.idx,col.minx+2*cut.idx - 1] == cut.pt 
  f <- rep(FALSE, x$nb.of.roi)
  f[list.roi.idx[f.roi]] <- TRUE
  x <- .suppress.roi(x,f)
  list.roi.idx <- (0:sum(f))[-1]
  
  if (length(list.roi.idx)==0) {
    message("no ROI to display @", cut.lab," = ", round(cut.pt,2),postfix)
    return(invisible(x))
    }
  t.mat <- ref.add (x$ref.pseudo,orientation = orientation, origin = origin,
                    new.ref.pseudo = "display.ref")
  obj_ <- struct.in.new.ref(x, "display.ref", t.mat)
  ext.pt <- round(get.extreme.pt(obj_),6)
  
  if (!is.null(abs.rng))  ext.pt[1, ] <- .range.intersection (abs.rng,ext.pt[1, ], c(TRUE,TRUE))
  if (!is.null(ord.rng))  ext.pt[2, ] <- .range.intersection (ord.rng,ext.pt[2, ], c(TRUE,TRUE))
  message1 <- F
  if(all(!is.na(ext.pt))){
    if (add) {
      ext.pt[1, ] <- .range.intersection (par("usr")[1:2],ext.pt[1, ], c(TRUE,TRUE))
      ext.pt[2, ] <- .range.intersection (par("usr")[3:4],ext.pt[2, ], c(TRUE,TRUE))
    }
  } else{ message1 <- T}
  if(any(is.na(ext.pt))){
    f <- rep(FALSE, x$nb.of.roi)
    x <- .suppress.roi(x,f)
    list.roi.idx <- (0:sum(f))[-1]
    if (message1){
      chr1 <- chr2 <-""
      if (!is.null(abs.rng)) chr1 <- paste0(", ",coord_[1], " in [", round(abs.rng[1],2), " ", round(abs.rng[2],2),"]")
      if (!is.null(ord.rng)) chr2 <- paste0(", ",coord_[2], " in [", round(ord.rng[1],2), " ", round(ord.rng[2],2),"]")
    } else {
      chr1 <- paste0(", ",coord_[1], " in [", round(par("usr")[1],2), " ", round(par("usr")[2],2),"]")
      chr2 <- paste0(", ",coord_[2], " in [", round(par("usr")[3],2), " ", round(par("usr")[4],2),"]")
    }
    message("no ROI to display @", cut.lab," = ", round(cut.pt,2),postfix, chr1, chr2)
    return(invisible(x))
  }
  M <- get.rigid.M(t.mat, x$ref.pseudo, "display.ref") %*% x$ref.from.contour
  dxyz <-round(as.vector(abs(dxyz %*% t(M)))[1:3],6)
  obj_$roi.data <- lapply(obj_$roi.data ,function(roi.L) {roi.L[sapply(roi.L,function(l){
    rg.xyz <- round(as.matrix(cbind(l$pt,1)) %*% t(obj_$ref.from.contour),6)
    (!is.na(.range.intersection(range(rg.xyz[,1]),  ext.pt[1, ], egal = c(TRUE,TRUE))[1])) &
      (!is.na(.range.intersection(range(rg.xyz[,2]), ext.pt[2, ], egal = c(TRUE,TRUE))[1])) &
      (!is.na(.range.intersection(range(rg.xyz[,3]),  0 + dxyz[3]*c(-0.5,0.5))[1]) |
         all(range(rg.xyz[,3]) == 0))})]})
  le.f <- sapply(obj_$roi.data, length) != 0
  obj_ <- .suppress.roi(obj_,le.f)
  list.roi.idx <- (0:sum(le.f))[-1]

  if (length(list.roi.idx)==0) {
    message("no ROI to display @", cut.lab," = ", round(cut.pt,2),postfix, ", with ",
            coord_[1], " in [",paste(round(ext.pt[1, ],2), collapse=","),"] and ",coord_[2], 
            " in [", paste(round(ext.pt[2, ],2), collapse=","),"]")
    return(invisible(x))
  }

  
  if (!all(obj_$ref.from.contour[1:2,3]== c(0,0))){
    
    for(j in list.roi.idx){
      pt <- round(as.matrix(cbind(do.call(rbind.data.frame, lapply(obj_$roi.data[[j]], function(l) l$pt)),1)) %*% 
        t(obj_$ref.from.contour),6)
      rgx <-.range.intersection(pt[,1], ext.pt[1,], egal = c(TRUE,TRUE)) 
      rgy <-.range.intersection(pt[,2], ext.pt[2,], egal = c(TRUE,TRUE))
      if (is.na(rgx[1]) | is.na(rgy[1])){
        obj_$roi.data[[j]] <- list()
      } else if(nrow(pt)>1){
        pt.max <- c(rgx[length(rgx)]+ dxyz[1],rgy[length(rgy)]+dxyz[2],0)
        pt000 =c(rgx[1]-dxyz[1],rgy[1]-dxyz[2],0)
        back.vol <- vol.create(n.ijk=ceiling((pt.max-pt000)/dxyz)+1, dxyz = dxyz,
                               pt000 =pt000,
                               ref.pseudo = "display.ref",default.value = FALSE)
        b <- bin.from.roi(back.vol,obj_,roi.idx =j, verbose = FALSE, 
                          description = paste0(obj_$roi.info$roi.pseudo[j],"@", cut.lab," = ", round(cut.pt,2)))
        if (any(b$vol3D.data)) {
          obj_$roi.data[[j]] <- .display.roi.data.from.bin(b)
          
        } else{obj_$roi.data[[j]] <- list()}
      }
      
    }
    le.f <- sapply(obj_$roi.data, length) != 0
    obj_ <- .suppress.roi(obj_,le.f)
    list.roi.idx <- (0:sum(le.f))[-1]
    
    if (length(list.roi.idx)==0) {
      message("no ROI to display @", cut.lab," = ", round(cut.pt,2),postfix)
      return(invisible(x))
    }
    
    
    obj_$ref.from.contour <- diag(4)
    obj_$thickness <- dxyz[3]
    
    for (i in 1:obj_$nb.of.roi) {
      if (length(obj_$roi.data[[i]])>1) {
        same.k.roi <- 1:length(obj_$roi.data[[i]])
        for (j in same.k.roi){
          ptj<- obj_$roi.data[[i]][[j]]$pt
          roi.index.k <-same.k.roi[same.k.roi!=j]
          r <- unique (sapply (roi.index.k, function (k) {
            if ( castlow.str (obj_$roi.data[[i]][[k]]$type) == "point") return(k)
            if ( castlow.str (obj_$roi.data[[i]][[k]]$type) != "closedplanar") return(NA)
            ptk <- obj_$roi.data[[i]][[k]]$pt
            keep <- .pt.in.polygon (ptj[ 1,1], ptj[ 1,2],ptk[ ,1], ptk[ ,2]) > 0.5
            return (ifelse (any(keep), k,NA))}))
          r <- r[!is.na (r)]
          obj_$roi.data[[i]][[j]]$level <- ifelse (length(r)!=0, length(r), 0)
        } 
      } else {obj_$roi.data[[i]][[1]]$level <- 0}
      
    }
  }
  
  ext.pt <- do.call(rbind,lapply(obj_$roi.data, function(L) 
    do.call(rbind,lapply(L, function(l) l$pt))))
  ext.pt <- t(apply(ext.pt,2,range))
  ext.pt[1, ] <- ext.pt[1, ]+c(-0.5,0.5)
  ext.pt[2, ] <- ext.pt[2, ]+c(-0.5,0.5)

  if (!is.null(abs.rng)) ext.pt[1,] <- range(abs.rng) #+ vol_$dxyz[1]*c(0.5,-0.5)
  if (!is.null(ord.rng)) ext.pt[2,] <- range(ord.rng) #+ vol_$dxyz[2]*c(0.5,-0.5)
  
  if (!add){  
    if (is.null(args[['main']]))
      args[['main']] <- paste0(x$modality," (",x$description, ") @", cut.lab," = ", round(cut.pt,2),postfix)
    args[['type']] <- "n"
    args[['x']] <- as.numeric(ext.pt[1,1:2])
    args[['y']] <- as.numeric(ext.pt[2,1:2])
    if (is.null(args[['xaxs']])) args[['xaxs']] <- "i"
    if (is.null(args[['yaxs']])) args[['yaxs']] <- "i"
    if (is.null(args[['xlab']])) args[['xlab']] <- coord_[1]
    if (is.null(args[['ylab']])) args[['ylab']] <- coord_[2]
    args[['xlim']] <- as.numeric(ext.pt[1,1:2])
    args[['ylim']] <- as.numeric(ext.pt[2,1:2])
    
    if (flip) args[['xlim']] <- rev( args[['xlim']])
    if (flop) args[['ylim']] <- rev( args[['ylim']])
    
    args[['asp']] <- 1
    do.call(plot,args)
    abs.usr <- par("usr")[1:2]
    ord.usr <- par("usr")[3:4]
    rect(abs.usr[1], ord.usr[1], abs.usr[2], ord.usr[2], col = bg)
    
  }
  
  # col <- args[["col"]];
  # if (!is.null(col)) obj_$roi.info$color <- c(col,rep(col[length(col)], nrow(obj_$roi.info)-length(col)))
  if (is.null(args[['lwd']])) args[['lwd']] <- 1
  for (j in list.roi.idx) {
    
    for (nb in 1:length(obj_$roi.data[[j]])){
      type <- castlow.str (obj_$roi.data[[j]][[nb]]$type)
      pt <- as.matrix(cbind(obj_$roi.data[[j]][[nb]]$pt,1)) %*% t(obj_$ref.from.contour)
      if (type=="closedplanar" | type=="openplanar"){
        lines (pt[,1], pt[,2], 
               col= obj_$roi.info$color[j], lwd=args[['lwd']])
      } else if ((type=="point") &
                 (round(obj_$roi.data[[j]][[nb]]$pt$z[1],6) == 0)) {
        points(pt[,1], pt[,2], 
               col= obj_$roi.info$color[j], pch="+",cex=1)
      }
    }
  }
  
  return(invisible(struct.in.new.ref(obj_, x$ref.pseudo,t.mat)))
}

#################################################################################
#' @rdname plot
#' @export
plot.mesh <- function (x,..., view.type = "trans", 
                       view.coord = NULL, 
                       flip = FALSE, flop = FALSE, col = "#ff0000") {
  
  args <- tryCatch(list(...), error=function(e)list())
  add <- args[["add"]]; if (is.null(add)) add <- FALSE; args[["add"]]<-NULL
  bg <- args[["bg"]]; if (is.null(bg)) bg <- "#000000"; args[["bg"]] <- NULL
  
  coord.map <-c('xy', 'yx', 'xz', 'zx', 'yz', 'zy', "trans","front","sagi")
  m <- match(view.type[1],coord.map)
  if (is.na(m)) stop("Undefined abcisse or ordinate. view.type must be set to ", paste(coord.map, collapse=", "))
  
  if (is.null(x$mesh)) stop("empty x$mesh")
  
  if (is.null(view.coord)) view.coord <-  apply(get.extreme.pt(x),1,mean)
  
  abs.rng <- args[['xlim']]
  ord.rng <- args[['ylim']]
  
  if (view.type=="trans") {
    view.type <- "xy"
    flop <- TRUE
    flip <- FALSE
    if (!is.null(ord.rng)) ord.rng <- rev(range(ord.rng))
  } else if (view.type=="front") {
    flip <- FALSE
    flop <- FALSE
    view.type <- "xz"
  } else if (view.type=="sagi")  {
    view.type <- "zy" 
    flip <- FALSE
    flop <- TRUE
    if (!is.null(ord.rng)) ord.rng <- rev(range(ord.rng))
  } else{
    if (!is.null(abs.rng)) {flip <- abs.rng[1]> abs.rng[2]}
    if (!is.null(ord.rng)) {flop <- ord.rng[1]> ord.rng[2]}
  }
  coord.label <- c("x","y","z")
  coord <- match(unique(unlist(strsplit(view.type,""))), coord.label)
  postfix <- " mm"
  
  if (!is.null(ord.rng)) ord.rng <- range(ord.rng)
  if (!is.null(abs.rng)) abs.rng <- range(abs.rng)
  
  coord_ <- coord[!is.na(coord)]
  cut.lab <- coord.label[-coord_]
  cut.idx <- (1:3)[-coord_]
  
  # view.coord <- args[["view.coord"]]; 
  # if (is.null(view.coord)) view.coord <- mean(as.numeric(get.extreme.pt(x)[cut.idx, ])); 
  if (is.na(match(length(view.coord),c(1,3)))) stop("view.coord must be of length 1 or 3")
  # args[["view.coord"]]<-NULL
  
  if (length(view.coord) == 1) {
    if (view.type=="xy" | view.type=="yx" | view.type == "trans") origin <- c(0,0,view.coord)
    if (view.type=="xz" | view.type=="zx" | view.type == "front") origin <- c(0,view.coord,0)
    if (view.type=="yz" | view.type=="zy" | view.type == "sagi") origin <-c(view.coord,0,0)
  } else {origin=view.coord; origin[coord] = 0}
  
  cut.pt <- origin[-coord_]
  coord_ <- coord.label[coord_]
  orientation <- switch (view.type,
                         "xy" = c(1,0,0,0,1,0),
                         "yx" = c(0,1,0,1,0,0),
                         "xz" = c(1,0,0,0,0,1),
                         "zx" = c(0,0,1,1,0,0),
                         "yz" = c(0,1,0,0,0,1),
                         "zy" = c(0,0,1,0,1,0))
  
  t.mat <- ref.add (x$ref.pseudo,orientation = orientation, origin = c(0,0,0),#origin,
                    new.ref.pseudo="display.ref")
  obj_ <- mesh.in.new.ref(x, "display.ref", T.MAT= t.mat,  alias = x$object.alias)
  ext.pt <- get.extreme.pt(obj_)
  obj_ <- struct.from.mesh(obj_, z = cut.pt, thickness = NULL, roi.color = col, alias = paste(x$object.alias, "cut"),
                           verbose = FALSE)
  
  list.roi.idx <- (0:length(obj_$roi.data[[1]]))[-1]
  
  if (length(list.roi.idx)==0) stop("no ROI to display")
  ext.pt <- get.extreme.pt(obj_)
  ext.pt[1, ] <- ext.pt[1, ]+c(-0.5,0.5)
  ext.pt[2, ] <- ext.pt[2, ]+c(-0.5,0.5)
  
  if (!is.null(abs.rng))  ext.pt[1, ] <- range(abs.rng)
  if (!is.null(ord.rng))  ext.pt[2, ] <- range(ord.rng) 
  

  if (!add){  
    if (is.null(args[['main']]))
      args[['main']] <- paste0(x$modality," (",x$description, ") @", cut.lab," = ", round(cut.pt,2),postfix)
    args[['type']] <- "n"
    args[['x']] <- as.numeric(ext.pt[1,1:2])
    args[['y']] <- as.numeric(ext.pt[2,1:2])
    if (is.null(args[['xaxs']])) args[['xaxs']] <- "i"
    if (is.null(args[['yaxs']])) args[['yaxs']] <- "i"
    if (is.null(args[['xlab']])) args[['xlab']] <- coord_[1]
    if (is.null(args[['ylab']])) args[['ylab']] <- coord_[2]
    args[['xlim']] <- as.numeric(ext.pt[1,1:2])
    args[['ylim']] <- as.numeric(ext.pt[2,1:2])
    
    if (flip) args[['xlim']] <- rev( args[['xlim']])
    if (flop) args[['ylim']] <- rev( args[['ylim']])
    
    args[['asp']] <- 1
    do.call(plot,args)
    abs.usr <- par("usr")[1:2]
    ord.usr <- par("usr")[3:4]
    rect(abs.usr[1], ord.usr[1], abs.usr[2], ord.usr[2], col = bg)
    
  }
  
  if (is.null(args[['lwd']])) args[['lwd']] <- 1

  for (nb in list.roi.idx){
    type <- castlow.str (obj_$roi.data[[1]][[nb]]$type)
    if (type=="closedplanar" | type=="openplanar"){
      lines (obj_$roi.data[[1]][[nb]]$pt$x, 
             obj_$roi.data[[1]][[nb]]$pt$y, 
             col= obj_$roi.info$color[1], lwd=args[['lwd']])
    } else if ((type=="point") &
               (round(obj_$roi.data[[1]][[nb]]$pt[1,3],6) == 0)) {
      points(obj_$roi.data[[1]][[nb]]$pt$x, 
             obj_$roi.data[[1]][[nb]]$pt$y, 
             col= obj_$roi.data$roi.info$color[1], pch="+",cex=1)
    }
  }
  
  return(invisible(struct.in.new.ref(obj_, x$ref.pseudo,t.mat)))
}
