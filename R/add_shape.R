#' Adding shape to a 3D volume.
#' @description The \code{add.shape} function adds the shape defined by espadon volume 
#' object of the modality "binary" or "weight" to a 3D volume.
#' @param obj Object of the "volume" class and a modality other than "binary" or "weight".
#' @param shape.bin  Object of the "volume" class and "binary" or "weight" modality, 
#' sharing the same voxels as \code{obj}..
#' @param shape.mean,shape.sd Positive numbers, representing the mean value  and 
#' the standard deviation of voxels identified by \code{shape.bin}.
#' @return Returns the ‘volume’ class object \code{obj}, in which the voxels 
#' identified by \code{shape.bin} have been replaced by a linear interpolation 
#' between the initial value and new values. These new values follow a normal 
#' distribution with mean \code{shape.mean} and standard deviation \code{shape.sd}.
#' The weights of the linear interpolation are defined by the voxels values of \code{shape.bin}.
#' @seealso \link[espadon]{bin.cuboid}, \link[espadon]{bin.cylinder}, \link[espadon]{bin.ellipsoid},
#' \link[espadon]{bin.from.roi}.
#' @examples
#' # Creation of a volume
#' CT <- vol.create (c(80, 80,40), c(1.2, 1.2, 2), 
#'                   pt000 = c(-50.4,-50.4,-75), modality = "ct", 
#'                   default.value = as.integer(-997), value.sd = 1)
#' # Creation of a new shape shape                
#' ellipsoid <- bin.ellipsoid(CT, center = c(-20.1, 0.1, -15), 
#'                            radius =  c(19.3, 20.2, 15.3))  
#'                            
#' # Incorporating form into the CT
#' CT <- add.shape (CT, shape.bin = ellipsoid, shape.mean = 100, shape.sd = 10)
#' 
#' plot(CT, view.coord =  c(-20.1, 0.1, -15))      
#' @export 
add.shape <- function(obj, shape.bin, shape.mean, shape.sd){
  
  if (!is (obj, "volume")) stop ("obj should be a volume class object.")
  if (!is (shape.bin, "volume")) stop ("shape.bin should be a volume class object.")
  if ((obj$modality=="binary") | (obj$modality=="weight")) 
    stop ("obj must be of modality other than binary or weight.")
  if ((shape.bin$modality!="binary") & (shape.bin$modality!="weight")) 
    stop ("shape.bin must be of binary or weight modality.")
  if (obj$ref.pseudo != shape.bin$ref.pseudo)
    stop("obj and shape.bin must have same ref.pseudo.")
  if (any (abs(obj$xyz.from.ijk[1:3, 1:3] - shape.bin$xyz.from.ijk[1:3, 1:3]) > 1e-6))   
    stop("obj and shape.bin must share the same grid.")
  
  vol3D <- shape.bin$vol3D.data
  idx <- which(vol3D>0, arr.ind = T)
  if (nrow(idx)==0) return( obj)
  vol3D[idx] <- vol3D[idx] *  rnorm(nrow(idx), abs(shape.mean), abs(shape.sd))
  #2D
  idx.c <- which(apply(abs(obj$xyz.from.ijk[1:3,1:3]),2,sum)==0)
  idx.r <-  which(apply(abs(obj$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    if (abs(obj$xyz0[1,idx.r])>1e-6) return(NULL)
    u <- obj$xyz.from.ijk
    u[idx.r,idx.c]<- 1
    volijk.from.xyz <- solve(u)
    volijk.from.xyz[idx.r,idx.c] <- 0
    u <- shape.bin$xyz.from.ijk
    u[idx.r,idx.c]<- 1
    shapeijk.from.xyz <- solve(u)
    shapeijk.from.xyz[idx.r,idx.c] <- 0
  } else {#3D
    volijk.from.xyz <- solve(obj$xyz.from.ijk)
    shapeijk.from.xyz <-  solve(shape.bin$xyz.from.ijk)
  }
  
  rg.volindex <- (volijk.from.xyz %*% shape.bin$xyz.from.ijk %*% shape.bin$cube.idx[,c(1,7)])[1:3,]
  
  if (any(abs(rg.volindex-round(rg.volindex))> 1e-6)) stop("obj and shape.bin must share voxels")
  rg.volindex <- round(rg.volindex)
  rg.volindex[,1] <- apply(cbind(rg.volindex[,1],obj$cube.idx[1:3,1]), 1 , max)
  rg.volindex[,2] <- apply(cbind(rg.volindex[,2],obj$cube.idx[1:3,7]), 1 , min)
  rg.shapeindex <- (shapeijk.from.xyz %*% obj$xyz.from.ijk %*% rbind(rg.volindex, c(1,1)))[1:3,]
  rg.shapeindex <- round(rg.shapeindex)
  obj$vol3D.data[(rg.volindex[1,1]:rg.volindex[1,2]) + 1,
                 (rg.volindex[2,1]:rg.volindex[2,2]) + 1,
                 (rg.volindex[3,1]:rg.volindex[3,2]) + 1] <- 
    (1 - shape.bin$vol3D.data[(rg.shapeindex[1,1]:rg.shapeindex[1,2]) + 1,
                              (rg.shapeindex[2,1]:rg.shapeindex[2,2]) + 1,
                              (rg.shapeindex[3,1]:rg.shapeindex[3,2]) + 1]) *  
    obj$vol3D.data[(rg.volindex[1,1]:rg.volindex[1,2]) + 1,
                   (rg.volindex[2,1]:rg.volindex[2,2]) + 1,
                   (rg.volindex[3,1]:rg.volindex[3,2]) + 1] +
    vol3D[(rg.shapeindex[1,1]:rg.shapeindex[1,2]) + 1,
          (rg.shapeindex[2,1]:rg.shapeindex[2,2]) + 1,
          (rg.shapeindex[3,1]:rg.shapeindex[3,2]) + 1] 
  
  obj$max.pixel <- max(obj$vol3D.data, na.rm=TRUE)
  obj$min.pixel <- min(obj$vol3D.data, na.rm=TRUE)
  
  return(obj)
}