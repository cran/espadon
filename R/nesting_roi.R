#' Restrict volume to RoI
#' @description The \code{nesting.roi} function restricts a "volume" class 
#' object to the rectangular parallelepiped circumscribed to the chosen RoI.
#' @param vol "volume" class object
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
#' @param vol.restrict Boolean. If \code{vol.restrict = TRUE}, the rectangular
#' parallelepiped circumscribed to the chosen RoI, enlarged by xyz.margin cannot
#' exceed the initial volume.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, \code{struct$ref.pseudo}
#' must be equal to \code{vol$ref.pseudo}.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be that of the \code{vol}.
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
nesting.roi <- function (vol, struct, roi.name = NULL, roi.sname = NULL, 
                         roi.idx = NULL, xyz.margin = c (0, 0, 0),
                         vol.restrict = FALSE, 
                         T.MAT = NULL, alias = "", description = NULL) {
  
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  # if (length (roi.idx) != 1) {
  #   warning ("multiple names or no names forbidden.")
  #   return (NULL)
  # }
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if (!is (struct, "struct")){
    warning ("struct should be a struct class object.")
    return (NULL)
  }
  
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol3D.data.")
    return (NULL)
  }
  if(is.null(struct$roi.data)){
    warning ("empty roi.data.")
    return (NULL)
  }
  if (length (xyz.margin)<3) {
    warning ("xyz.margin must have a length of 3.")
    return (NULL)
  }
  
  
  if (all (sapply (roi.idx,function(i) {length(struct$roi.data[[i]])==0}))) return (NULL)
  
  if (!is.null(vol$local.gridx)) {
    vol$local.gridx <- NULL
    vol$local.gridy <- NULL
  }
  
  t.mat <- ref.cutplane.add(vol, T.MAT=T.MAT)
  S_ <- struct.in.new.ref(struct,new.ref.pseudo=paste0(vol$ref.pseudo,"m"), t.mat)
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo=paste0(vol$ref.pseudo,"m"), t.mat)

  if (is.null (S_) | is.null (vol_)){
    warning ("different ref.pseudo.")
    return (NULL)
  }
  pt.min <- c (min (S_$roi.info$min.x[roi.idx]) - xyz.margin[1],
               min (S_$roi.info$min.y[roi.idx]) - xyz.margin[2],
               min (S_$roi.info$min.z[roi.idx]) - xyz.margin[3])
  pt.max <- c (max (S_$roi.info$max.x[roi.idx]) + xyz.margin[1],
               max (S_$roi.info$max.y[roi.idx]) + xyz.margin[2],
               max (S_$roi.info$max.z[roi.idx]) + xyz.margin[3])
  
  if (vol.restrict){
    vol.pt <- vol_$xyz.from.ijk %*% vol_$cube.idx[ , c (1,7)]
    pt.min <- c (max(pt.min[1], min(vol.pt [1,])), max(pt.min[2], min(vol.pt [2,])), max(pt.min[3], min(vol.pt [3,])))
    pt.max <- c (min(pt.max[1], max(vol.pt [1,])), min(pt.max[2], max(vol.pt [2,])), min(pt.max[3], max(vol.pt [3,])))
  }
  

  if (is.null(description)) description<- paste (vol$description,"restricted to", paste (struct$roi.info$roi.pseudo[roi.idx], collapse="|"))
  if (pt.min[1]<=pt.max[1] & pt.min[2]<=pt.max[2] & pt.min[3]<=pt.max[3]){
    vol_ <- nesting.cube(vol=vol_, pt.min=pt.min, pt.max=pt.max, alias = alias, description=description )
    return(vol.in.new.ref(vol_,new.ref.pseudo=vol$ref.pseudo, t.mat,alias = alias, description=description))
  }
  return (NULL)
}