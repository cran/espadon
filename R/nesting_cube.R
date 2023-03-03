#' Restriction of a volume to a rectangular parallelepiped
#' @description The \code{nesting.cube} function restricts or increases
#' a volume to the rectangular parallelepiped defined by its 2 extreme vertices.

#' @param obj object of class volume or mesh.
#' @param pt.min minimum x, y, z coordinates of the rectangular parallelepiped vertex.
#' @param pt.max maximum x, y, z coordinates of the rectangular parallelepiped vertex.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. If
#' the \code{description = NULL} (default value), it will be set to \code{obj$description}.
#' @param ... Additional arguments \code{vol} (depracated), replaced by \code{obj}.
#' @return Returns a "volume" class object, in which 3D volume is restricted 
#' or increased to be circumscribed to the requested rectangular parallelepiped. 
#' If the created volume exceeds the initial volume, new voxels are set to \code{NA}.
#' @seealso \link[espadon]{add.margin}, \link[espadon]{nesting.roi} and  
#' \link[espadon]{nesting.bin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = "ct", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' # Calculation of the new CT restricted to the parallelepiped reduced by 10 mm.
#' pt.CT <- get.extreme.pt (CT) # extreme points of CT
#' new.pt.CT <- pt.CT + matrix (rep (c (+ 12, -12), 3), ncol = 2, byrow = TRUE)
#' new.CT <- nesting.cube (CT, new.pt.CT$min, new.pt.CT$max, alias = "new CT")
#' \dontrun{  
#' # check for change
#' display.3D.stack (CT)
#' display.3D.stack (new.CT, line.col="red")
#' }
#' @export
#' @importFrom methods is

nesting.cube <- function (obj, pt.min, pt.max, alias = "", description = NULL,...){
  
  passed <- names(as.list(match.call())[-1])
  args <- list(...)
  if (!("obj" %in% passed)){
    if (is.null(args[['vol']])) stop('argument "obj" is missing, with no default')
    obj <- args[['vol']]
  }
  
  if (!(is (obj, "volume") | is (obj, "mesh"))) {
    warning ("obj must be an object of class volume or mesh.")
    return (NULL)
  }
  
  if (is (obj, "volume")){
    if(is.null(obj$vol3D.data)){
      warning ("empty obj$vol3D.data.")
      return (NULL)
    }
    
    struct.cube <- obj$cube.idx
    struct.cube[1, c(1,4,5,8)] <- pt.min[1]
    struct.cube[1, c(2,3,6,7)] <- pt.max[1]
    struct.cube[2, c(1,2,5,6)] <- pt.min[2]
    struct.cube[2, c(3,4,7,8)] <- pt.max[2]
    struct.cube[3, 1:4] <- pt.min[3]
    struct.cube[3, 5:8] <- pt.max[3]
    
    struct.in.obj.cube.idx <- round(solve (obj$xyz.from.ijk) %*% struct.cube,6)
    rg.i <- c(floor(min(struct.in.obj.cube.idx[1,])):ceiling(max(struct.in.obj.cube.idx[1,])))
    rg.j <- c(floor(min(struct.in.obj.cube.idx[2,])):ceiling(max(struct.in.obj.cube.idx[2,])))
    rg.k <- c(floor(min(struct.in.obj.cube.idx[3,])):ceiling(max(struct.in.obj.cube.idx[3,])))
    pt000 <- c(rg.i[1],rg.j[1],rg.k[1], 1) %*% t(obj$xyz.from.ijk)
    
    flag.i <- rg.i>=0 & rg.i<obj$n.ijk[1]
    flag.j <- rg.j>=0 & rg.j<obj$n.ijk[2]
    flag.k <- !is.na (match(rg.k, obj$k.idx))
    
    if (!is.null(obj$local.gridx)) {
      obj$local.gridx <- NULL
      obj$local.gridy <- NULL
    }
    
    new.obj <- vol.copy (vol = obj, alias = alias, description = description)
    
    new.obj$n.ijk <- c (length( rg.i), length (rg.j), length (rg.k))
    new.obj$xyz.from.ijk[ ,4] <- pt000
    new.obj$k.idx <- 0:(new.obj$n.ijk[3]-1)
    new.obj$xyz0  <- matrix ((as.matrix (expand.grid (0, 0, new.obj$k.idx,1))%*% t(new.obj$xyz.from.ijk))[ ,1:3],ncol=3)
    
    new.obj$cube.idx <- matrix ( c(0, 0, 0, 1, new.obj$n.ijk[1]-1, 0, 0, 1,
                                   new.obj$n.ijk[1]-1, new.obj$n.ijk[2]-1, 0, 1, 0, new.obj$n.ijk[2]-1, 0, 1,
                                   0, 0, new.obj$n.ijk[3]-1, 1, new.obj$n.ijk[1]-1, 0, new.obj$n.ijk[3]-1, 1,
                                   new.obj$n.ijk[1]-1, new.obj$n.ijk[2]-1, new.obj$n.ijk[3]-1, 1, 0, new.obj$n.ijk[2]-1, new.obj$n.ijk[3]-1, 1),
                                 nrow=4, byrow= FALSE)
    
    new.obj$vol3D.data <- array(NA, dim=new.obj$n.ijk)
    
    rg.i.to.complete <- (1:new.obj$n.ijk[1])[flag.i]
    rg.j.to.complete <- (1:new.obj$n.ijk[2])[flag.j]
    rg.k.to.complete <- (1:new.obj$n.ijk[3])[flag.k]
    
    new.obj$vol3D.data[rg.i.to.complete, rg.j.to.complete, rg.k.to.complete] <- obj$vol3D.data[rg.i[flag.i] + 1,
                                                                                               rg.j[flag.j] + 1,
                                                                                               rg.k[flag.k] + 1]
    #new.obj$vol3D.data[new.obj$vol3D.data<=1e-6] <- 0
    if (any(!is.na(new.obj$vol3D.data))){
      new.obj$min.pixel <- min(new.obj$vol3D.data, na.rm = TRUE)
      new.obj$max.pixel <- max(new.obj$vol3D.data, na.rm = TRUE)
    } else {
      new.obj$min.pixel <- NA
      new.obj$max.pixel <- NA
    }
    
    return(new.obj)
  } else {
    vb.idx<- which(as.numeric( obj$mesh$vb[1,])>=as.numeric(pt.min[1]) & 
                     as.numeric( obj$mesh$vb[2,])>=as.numeric(pt.min[2]) & 
                     as.numeric( obj$mesh$vb[3,])>=as.numeric(pt.min[3]) &
                     as.numeric( obj$mesh$vb[1,])<=as.numeric(pt.max[1]) & 
                     as.numeric( obj$mesh$vb[2,])<=as.numeric(pt.max[2]) & 
                     as.numeric( obj$mesh$vb[3,])<=as.numeric(pt.max[3]))
    
 
    f1 <- !is.na(match(obj$mesh$it[1,],vb.idx)) | 
      !is.na(match( obj$mesh$it[2,],vb.idx)) | 
      !is.na(match( obj$mesh$it[3,],vb.idx))
    pt.idx1 <- sort(unique(as.numeric(obj$mesh$it[,f1])))
    
    
    m1 <-obj
    
    
    m1$mesh$vb <- obj$mesh$vb[,pt.idx1]
    m1$mesh$it <- obj$mesh$it[,f1]
    m1$mesh$it[1,] <- match(obj$mesh$it[1,f1],pt.idx1)
    m1$mesh$it[2,] <- match(obj$mesh$it[2,f1],pt.idx1)
    m1$mesh$it[3,] <- match(obj$mesh$it[3,f1],pt.idx1)
    m1$mesh$normals <- obj$mesh$normals[,pt.idx1]
    m1$mesh$remface <- obj$mesh$remface[f1]
    m1$nb.faces<- ncol(m1$mesh$it)
    if (!is.null(description)) m1$description <- description
    
    m1$file.basename <- ""
    m1$file.dirname <- ""
    m1$object.name <- alias
    m1$object.alias <- alias
    m1$object.info <- NULL
    m1$ref.object.alias <- NULL
    m1$ref.object.info <- NULL
    
    if (alias=="") return(m1)
    return(.set.ref.obj(m1,list(obj)))
  }
  
  
}