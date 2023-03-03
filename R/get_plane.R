#' Extracting a plane from a volume
#' @description The \code{get.plane} function extracts a plane from a "volume" 
#' class object.
#' @param vol "volume" class object.
#' @param origin Vector of x, y, z coordinates, representing the origin of the 
#' plane to extract. If \code{interpolate = FALSE}, these coordinates are replaced 
#' by the coordinates of the voxel closest to \code{origin}.
#' @param plane.orientation Vector orientation of the plane in the \code{vol} 
#' frame of reference, composed by the 2 vectors coordinates of the orthonormal
#' basis of the plane. First vector is x-axis, and second one is y-axis.
#' @param alias \code{$object.alias} of the created object.
#' @param xgrid Vector, representing the grid of the plane abscissa. See note.
#' @param ygrid Vector, representing the grid of the plane ordinates. See note.
#' If \code{ygrid = NULL}, the ordinate is the line intercepting the volume and the
#' step is set to the projection of \code{vol$dxyz} onto the ordinate orientation.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a
#' trilinear interpolation of the value of the voxels, relative to the values of
#' adjacent voxels, is performed.
#' @note \emph{Determination of axes} : 
#' \itemize{
#' \item the x-axis has \code{plane.orientation[1:3]} as unit vector. 
#' \item the y-axis has \code{plane.orientation[4:6]} as unit vector. 
#' \item If \code{xgrid} is not {NULL}, \code{origin + x.grid * plane.orientation[1:3]}  
#' are the coordinates of the points on the x axis.
#' \item If \code{ygrid} is not {NULL}, \code{origin + y.grid * plane.orientation[4:6]}  
#' are the coordinates of the points on the y axis.
#' \item If \code{xgrid} or \code{ygrid} are NULL, they are determined to represent 
#' as closely as possible the initial volume in the required cut. 
#' }
#' @return Returns a "volume" class object, containing only a single plane, 
#' at \code{k = 0}, in the same frame of reference as \code{vol}.
#' This returned object has 2 new fields \code{local.xgrid}, and \code{local.ygrid}, 
#' representing the local grids of the abscissa (columns) and ordinate (rows) 
#' of the plane.
#' @return Returns \code{NULL} if plane doesn't exist.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "mr", dxyz = rep (step, 3))
#' MR <- patient$mr[[1]] 
#'     
#' # mid-volume point 
#' mid.point <- apply (get.extreme.pt (MR),1,mean)
#' 
#' plane <- get.plane (MR, origin = mid.point, interpolate = TRUE)
#' display.kplane (plane, interpolate = FALSE)
#' 
#' plane <- get.plane (MR, origin = mid.point, xgrid = seq (-50, 50, 1), 
#'                     ygrid = seq (-50, 50, 1), interpolate = TRUE)
#' display.kplane (plane, interpolate = FALSE)
#'
#' # 3 points on the inclined plane
#' pts <- t ((MR$xyz.from.ijk %*% MR$cube.idx) [1:3 , c (1, 2, 7)])
#' orientation <- orientation.create (A = pts[1,], B = pts[2,], C = pts[3,])
#' origin <- apply (pts, 2, mean)
#' plane <- get.plane (MR, origin = origin, 
#'                     plane.orientation = orientation, interpolate = TRUE)
#' display.kplane (plane, interpolate = FALSE)
#' 
#' orientation <- orientation.create (A = c (0, 0, 0) , B = c (1, 1, 0), 
#'                                    C = c (-1, 1, 0))
#' plane <- get.plane (MR, origin = origin, 
#'                     plane.orientation = orientation, interpolate = TRUE) 
#' display.kplane (plane, interpolate = FALSE)

#' @export
#' @importFrom methods is
get.plane <- function(vol, origin = c (0, 0, 0),
                      plane.orientation = c (1, 0, 0, 0, 1, 0), 
                      alias = "plane.n",
                      xgrid = NULL, ygrid = NULL, interpolate = TRUE){
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  plane.ref="plane.ref"

  
  if (is.null (vol)) return (NULL)
  
  if (length (plane.orientation)!=6 & length (plane.orientation)!=9) {
    warning ("bad plane.orientation length.")
    return (NULL)
  }
  #si pas d'interpolation, on prend le point le plus proche qui existe dans le volume
  real.origin <- origin
  if (!interpolate) {
    
    idx.c <- which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
    idx.r <-  which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),1,sum)==0)
    u <- vol$xyz.from.ijk 
    if (length(idx.c)>0) {
      u[idx.r,idx.c]<- 1
    }
    real.origin <- (floor ((c(origin,1) %*% t (solve (u))) + 0.5) %*% t (vol$xyz.from.ijk ))[1:3]
  }
  orientation <- plane.orientation
  #if (rev.k) orientation <- c(plane.orientation,-vector.product(plane.orientation[1:3],plane.orientation[4:6]))
  
  t.mat <-ref.add (vol$ref.pseudo, origin = real.origin, orientation=orientation,
                   new.ref.pseudo = plane.ref)
  t.mat$ref.info[t.mat$ref.info$ref.pseudo==vol$ref.pseudo, ]$ref <- vol$frame.of.reference
  # TM <- t.mat$matrix.list[[paste (plane.ref,vol$ref.pseudo,sep="<-")]]

  v <- vol.in.new.ref(vol, plane.ref, t.mat)
  new.dxyz <- v$dxyz
  
  #2D
  idx.c <- which(apply(abs(v$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
  idx.r <-  which(apply(abs(v$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    if (abs(v$xyz0[1,idx.r])>1e-6) return(NULL)
    u <- v$xyz.from.ijk 
    u[idx.r,idx.c]<- 1
    ijk.from.xyz <- solve(u)
    ijk.from.xyz[idx.r,idx.c] <- 0
  } else {ijk.from.xyz <- solve(v$xyz.from.ijk)}
  
  # ijk.from.xyz <- solve(v$xyz.from.ijk)
  ext.pt <- as.matrix(rbind(get.extreme.pt(v),c(1,1)))
  dif.ext.pt <- ext.pt[,2]-ext.pt[,1]
  

  mat <- ijk.from.xyz %*% (ext.pt*matrix(c(1,0,0,1,1,0,0,1), ncol=2))
  ds <- max(abs(mat[,2]-mat[,1]))
  if (ds>0) new.dxyz[1] <- dif.ext.pt[1]/max(abs(mat[,2]-mat[,1]))
  
  mat <- ijk.from.xyz %*% (ext.pt*matrix(c(0,1,0,1,0,1,0,1), ncol=2))
  ds <- max(abs(mat[,2]-mat[,1]))
  if (ds>0) new.dxyz[2] <- dif.ext.pt[2]/max(abs(mat[,2]-mat[,1]))
  
  mat <- ijk.from.xyz %*% (ext.pt*matrix(c(0,0,1,1,0,0,1,1), ncol=2))
  ds <- max(abs(mat[,2]-mat[,1]))
  if (ds>0) new.dxyz[3] <- dif.ext.pt[3]/max(abs(mat[,2]-mat[,1]))
  
  ext.pt <- ext.pt[1:2,]/ matrix(new.dxyz[c(1,2,1,2)],ncol=2) 
  # round.ext.pt <- round(ext.pt + sign(ext.pt) * 0.5*(abs(ext.pt - round(ext.pt)) > 1e-6))
  round.ext.pt <- round(ext.pt - sign(ext.pt) * 0.5*(abs(ext.pt - round(ext.pt)) > 1e-6))
  if (is.null(xgrid)){
    # creation de la grille incluant le point origin
    i.idx <- round.ext.pt[1,1]:round.ext.pt[1,2]
    new.grid.i <- i.idx * abs(new.dxyz [1])
    if (sign(new.dxyz [1])==-1) new.grid.i <- rev(new.grid.i)
  } else{
    new.grid.i <- xgrid
    if (length(new.grid.i)>1) new.dxyz [1] <- new.grid.i[2]-new.grid.i[1]
  }
  if (is.null(ygrid)){
    j.idx <- round.ext.pt[2,1]:round.ext.pt[2,2]
    new.grid.j <- j.idx * abs(new.dxyz [2])
    if (sign(new.dxyz [2])==-1) new.grid.j <- rev(new.grid.j)
  } else{
    new.grid.j <- ygrid
    if (length(new.grid.j)>1) new.dxyz [2] <- new.grid.j[2]-new.grid.j[1]
  }
  TM <- get.rigid.M(t.mat,plane.ref, vol$ref.pseudo)
  if (all (vector.product (orientation[1:3], orientation[4:6])==- as.numeric (TM[1:3,3]))) 
    new.dxyz[3] <- -new.dxyz[3]

  
  back.vol <- vol.create (pt000=  c (new.grid.i[1],new.grid.j[1], 0),
                          dxyz = new.dxyz,
                          n.ijk= c (length(new.grid.i),length(new.grid.j),1),
                          ref.pseudo = plane.ref, frame.of.reference = t.mat$ref.info[t.mat$ref.info$ref.pseudo == plane.ref, ]$ref,
                          alias = alias, number = vol$number, modality = vol$modality,  description = vol$description)
  
  new.vol <- vol.regrid (v, back.vol, T.MAT=t.mat, interpolate = interpolate, alias = alias, verbose = FALSE)
  
  new.vol$object.alias <- vol$object.alias
  new.vol$object.info <- vol$object.info

  
  new.vol <- vol.in.new.ref (new.vol, vol$ref.pseudo, T.MAT=t.mat, alias = alias)
  new.vol$local.gridx <- new.grid.i
  new.vol$local.gridy <- new.grid.j
  
  
  
  return (new.vol)
  
}

