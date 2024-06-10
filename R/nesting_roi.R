#' Restrict volume to RoI
#' @description The \code{nesting.roi} function restricts a "volume" class 
#' object to the rectangular parallelepiped circumscribed to the chosen RoI.
#' @param obj object of class volume or mesh.
#' @param struct "struct" class object.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Names or parts of names of the RoI in the \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Index of the RoI that belong to the \code{struct} object.
#' By default \code{roi.idx = NULL}. See Details.
#' @param xyz.margin Vector of length 3, representing the distances in mm to be added
#' to the x, y and z directions of the rectangular parallelepiped circumscribed
#' to the chosen RoI, in the cutting planes frame of reference. By default 
#' xyz.margin = c (0, 0, 0).
#' @param obj.restrict Boolean. Used if \code{obj} is of class"volume". If 
#' \code{obj.restrict = TRUE}, the rectangular parallelepiped circumscribed to 
#' the selected voxels, enlarged by xyz.margin cannot exceed the initial volume.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, \code{struct$ref.pseudo}
#' must be equal to \code{obj$ref.pseudo}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be that of the \code{obj}, plus 
#' "restricted to" the selected RoI.
#' @param ... Additional arguments such as \code{vol} (depracated), replaced by \code{obj}.
#' @return Returns a "volume" class object, in which 3D volume is limited to the
#' rectangular parallelepiped circumscribed to the chosen RoI, increased by the
#' requested margins.
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all set
#' to NULL, all RoI are selected.
#' @seealso \link[espadon]{add.margin}, \link[espadon]{nesting.cube} and  
#' \link[espadon]{nesting.bin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name = "brain", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' CT.brain <- nesting.roi (CT, S, roi.sname = "brain")
#' CT.brain.with.margin <- nesting.roi (CT, S, roi.sname = "brain",
#'                                          xyz.margin = c (10,10,10))
#' 
#' # display at the center of gravity of the cerebellum Gz
#' Gz <- S$roi.info [grep("^brain",S$roi.info$roi.pseudo),]$Gz
#' display.plane (bottom = CT.brain, view.coord = Gz,
#'                struct = S, bottom.col = pal.RVV (1000),
#'                bottom.breaks = seq (-1000, 1000, length.out = 1001),
#'                bg = "#00ff00",  interpolate  = FALSE, legend.shift = -20)
#' display.plane (bottom = CT.brain.with.margin,view.coord = Gz,
#'                struct = S,  bottom.col = pal.RVV (1000),
#'                bottom.breaks = seq(-1000, 1000, length.out = 1001),
#'                bg = "#00ff00", interpolate  = FALSE, legend.shift = -20)

#' @export
#' @importFrom methods is
nesting.roi <- function (obj, struct, roi.name = NULL, roi.sname = NULL, 
                         roi.idx = NULL, xyz.margin = c (0, 0, 0),
                         obj.restrict = FALSE, 
                         T.MAT = NULL, alias = "", description = NULL,...) {
  passed <- names(as.list(match.call())[-1])
  args <- list(...)
  if (!("obj" %in% passed)){
    if (is.null(args[['vol']])) stop('argument "obj" is missing, with no default')
    obj <- args[['vol']]
  }
  if (!is.null(args[['vol.restrict']])) obj.restrict <- args[['vol.restrict']]
  
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  # if (length (roi.idx) != 1) {
  #   warning ("multiple names or no names forbidden.")
  #   return (NULL)
  # }
  
  if (!(is (obj, "volume") | is (obj, "mesh"))) stop ("obj must be an object of class volume or mesh.")
  if (!is (struct, "struct")) stop ("struct should be a struct class object.")
  if (is.null(struct$roi.data)) stop ("empty roi.data.")
  if (all (sapply (roi.idx,function(i) {length(struct$roi.data[[i]])==0}))) stop ("empty roi.data.")
  if (length (xyz.margin)==1) xyz.margin <- rep(xyz.margin,3)[1:3]
  if (length (xyz.margin)!=3) stop ("xyz.margin must have a length of 1 or 3.")
 
  
  
  
  if (is.null(description)) description<- paste (obj$description,"restricted to", 
                                                 paste (struct$roi.info$roi.pseudo[roi.idx], collapse="|"))
  
  if(is (obj, "volume")){
  #--------------------- 
    if(is.null(obj$vol3D.data))  stop ("empty vol3D.data.")

    if (!is.null(obj$local.gridx)) {
      obj$local.gridx <- NULL
      obj$local.gridy <- NULL
    }
    
    t.mat <- ref.cutplane.add(obj, T.MAT=T.MAT)
    S_ <- struct.in.new.ref(struct,new.ref.pseudo=paste0(obj$ref.pseudo,"m"), t.mat)
    vol_ <- vol.in.new.ref(obj, new.ref.pseudo=paste0(obj$ref.pseudo,"m"), t.mat)
    
    if (is.null (S_) | is.null (vol_)){
      warning ("different ref.pseudo.")
      return (NULL)
    }
    pt.min <- c (min (S_$roi.info$min.x[roi.idx], na.rm=TRUE) - xyz.margin[1],
                 min (S_$roi.info$min.y[roi.idx], na.rm=TRUE) - xyz.margin[2],
                 min (S_$roi.info$min.z[roi.idx], na.rm=TRUE) - xyz.margin[3])
    pt.max <- c (max (S_$roi.info$max.x[roi.idx], na.rm=TRUE) + xyz.margin[1],
                 max (S_$roi.info$max.y[roi.idx], na.rm=TRUE) + xyz.margin[2],
                 max (S_$roi.info$max.z[roi.idx], na.rm=TRUE) + xyz.margin[3])
    
    if (obj.restrict){
      vol.pt <- vol_$xyz.from.ijk %*% vol_$cube.idx[ , c (1,7)]
      pt.min <- c (max(pt.min[1], min(vol.pt [1,])), max(pt.min[2], min(vol.pt [2,])), max(pt.min[3], min(vol.pt [3,])))
      pt.max <- c (min(pt.max[1], max(vol.pt [1,])), min(pt.max[2], max(vol.pt [2,])), min(pt.max[3], max(vol.pt [3,])))
    }

    if (pt.min[1]<=pt.max[1] & pt.min[2]<=pt.max[2] & pt.min[3]<=pt.max[3]){
      vol_ <- nesting.cube(obj=vol_, pt.min=pt.min, pt.max=pt.max, alias = alias, description=description )
      vol <- vol.in.new.ref(vol_,new.ref.pseudo=obj$ref.pseudo, t.mat,alias = alias, description=description)
      
      if (alias == "") return(vol)
      return(.set.ref.obj(vol,list(obj,struct),add=FALSE))
    }
    
  } else {#mesh
  ##############  
    ext <- get.extreme.pt (struct,roi.idx = roi.idx, ref.pseudo = obj$ref.pseudo, 
                           T.MAT = T.MAT)
    if (is.null(ext))  return (NULL)
    pt.min <- as.numeric(ext[,1]) - xyz.margin 
    pt.max <- as.numeric(ext[,2]) + xyz.margin
    
    if (pt.min[1]<=pt.max[1] & pt.min[2]<=pt.max[2] & pt.min[3]<=pt.max[3]){
      new.obj <- nesting.cube(obj=obj, pt.min=pt.min, pt.max=pt.max, alias = alias, description=description)
      if (alias =="") return(new.obj)
      return(.set.ref.obj(new.obj, list(struct)))
    }
  }
  return (NULL)
}