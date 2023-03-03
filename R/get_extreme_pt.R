#' Coordinates of the extreme points
#' @description The \code{get.extreme.pt} function returns the x, y, z coordinates
#' of the 2 extreme voxels of the rectangular parallelepiped, containing 
#' the objet \code{obj} of class volume, struct or mesh. These coordinates are given in 
#' the \code{ref.pseudo} frame of reference.
#' @param obj  object of class volume or struct or mesh.
#' @param ref.pseudo Pseudonym of the frame of reference in which you want the 
#' coordinates.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{ref.pseudo} must be equal to \code{obj$ref.pseudo}.
#' @param ... Additional arguments \code{min}, \code{max} if \code{obj} is of class 'volume'. 
#' Arguments \code{roi.name}, \code{roi.sname}, \code{roi.idx} if \code{obj} is 
#' of class 'struct'. Arguments \code{vol} (depracated), replaced by \code{obj}.
#' @return Returns a dataframe of min and max columns, and x, y and z rows. 
#' 
#' \itemize{
#' \item If \code{obj} is a member of the class volume: the returned dataframe 
#' represents the coordinates of the 2 extreme points of the rectangle parallelepiped 
#' including all the voxels such as \code{min <= obj$vol3D.data <= max}, 
#' if the arguments \code{min} or \code{max} exist, or including all the voxels otherwise.
#' \item If \code{obj} is a member of the class struct: the returned dataframe 
#' represents the coordinates of the 2 extreme points of the rectangular parallelepiped 
#' including all the selected RoI.
#' \item if \code{obj} is a member of the class mesh: the returned dataframe 
#' represents the coordinates of the 2 extreme points of the rectangular parallelepiped 
#' including all the mesh.
#' }
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = "ct", roi.name = "", dxyz = c (5, 5, 5))
#' CT <- patient$ct[[1]]
#'
#' # xyz extreme coordinate
#' get.extreme.pt (CT)
#' get.extreme.pt (CT, min = 0)
#' @export
#' @importFrom methods is
get.extreme.pt <- function (obj,ref.pseudo = obj$ref.pseudo, T.MAT = NULL, ...) {
  

  passed <- names(as.list(match.call())[-1])
  args <- list(...)
  if (!("obj" %in% passed)){
    if (is.null(args[['vol']])) stop('argument "obj" is missing, with no default')
    obj <- args[['vol']]
  }
  
  if (!(is (obj, "volume") | is (obj, "mesh") | is (obj, "struct"))) {
    warning ("obj must be an object of class volume or struct or mesh.")
    return (NULL)
  }
  M_ <- diag(4)
  if (ref.pseudo != obj$ref.pseudo) {
    M_ <- get.rigid.M(T.MAT,obj$ref.pseudo,ref.pseudo)
    if (is.null(M_)){
      warning ("different ref.pseudo. Enter T.MAT")
      return (NULL)
    }
  }
    
  if (is (obj, "volume")){
    if (is.null(args[['min']]) & is.null(args[['max']])){
      M <- M_ %*% obj$xyz.from.ijk
      ext_ <- (M %*% obj$cube.idx)[1:3, ]
      ext <- data.frame (min= apply (ext_,1,min),max= apply (ext_,1,max))
    } else {
      min=-Inf
      max=Inf
      if(!is.null(args[['min']])) min <- args[['min']]
      if(!is.null(args[['max']])) min <- args[['max']]  
      obj_ <- vol.in.new.ref(obj, new.ref.pseudo = ref.pseudo, T.MAT= T.MAT)
      pt <- get.xyz.from.index(which(obj_$vol3D.data >= min & obj_$vol3D.data <= max), obj_)
      if (is.null(pt))return (NULL)
      ext <- data.frame(min=apply(pt,2,min),max = apply(pt,2,max))
    }
    
    
  } else if (is (obj, "struct")){
    roi.name <- roi.sname <- roi.idx <- NULL
    if(!is.null(args[['roi.name']])) roi.name <- args[['roi.name']]
    if(!is.null(args[['roi.sname']])) roi.sname <- args[['roi.sname']]
    if(!is.null(args[['roi.idx']])) roi.idx <- args[['roi.idx']]
    
    list.roi.idx <- select.names (obj$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
    if (length(list.roi.idx)==0)   return (NULL)
    obj_ <- obj
    obj_$nb.of.roi <- length(list.roi.idx)
    obj_$roi.info <- obj$roi.info[list.roi.idx,]
    obj_$roi.obs <-  obj$roi.obs[list.roi.idx,]
    row.names(obj_$roi.info) <- NULL
    row.names(obj_$roi.obs) <- NULL
    obj_$roi.data <- obj$roi.data[list.roi.idx]
    obj_ <- struct.in.new.ref(obj_, new.ref.pseudo = ref.pseudo, T.MAT= T.MAT)
    ext <- data.frame(min = as.numeric(apply(obj_$roi.info[,c("min.x","min.y","min.z")],2,min)),
                      max = as.numeric(apply(obj_$roi.info[,c("max.x","max.y","max.z")],2,max)))
    
    
  } else {
    pt <-  M_ %*% obj$mesh$vb
    ext <- data.frame(min=apply(pt[1:3,],1,min),max = apply(pt[1:3,],1,max))
  }
  row.names (ext) <- c("x", "y", "z")
  return(ext)
}