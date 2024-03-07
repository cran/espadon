################################################################################
#' Adding or removing a margin to a volume
#' @description The \code{add.margin} function adds or subtracts a margin of the 
#' rectangular parallelepiped circumscribed by a volume.
#' @param vol "volume" class object.
#' @param xyz.margin Vector of the 3 positive or negative x, y and z margins in 
#' mm, in the frame of reference of volume cut planes.
#' @param alias Character string, \code{$alias} of the created object
#' @param description Character string, describing the created object.
#' If \code{description = NULL} (default value), it will be set to \code{vol$description}
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), in which 3D volume is restricted  or increased by the 
#' requested margins.
#' If the created volume exceeds the initial volume, new voxels are set to \code{NA}.
#' @seealso \link[espadon]{nesting.cube}, \link[espadon]{nesting.roi} and  
#' \link[espadon]{nesting.bin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "ct", roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#'
#' # Calculation of new volumes decreased by 10 mm in all directions.
#' new.CT <- add.margin (CT, xyz.margin = c (-10, -10, 10), alias = "new CT")

#' # display of the CT before and after, in the middle plane
#' z.mid <- apply (get.extreme.pt (CT), 1, mean)[3]
#' display.plane (bottom = CT, view.coord = z.mid, bottom.col = pal.RVV(1000),
#'                bg = "#00ffff", interpolate = FALSE)
#' display.plane (bottom = new.CT, view.coord = z.mid, bottom.col = pal.RVV(1000),
#'                bg = "#00ffff", interpolate = FALSE)
#'

#' @export
#' @importFrom methods is
add.margin <- function (vol, xyz.margin, alias="", description=NULL) {
  
  if (!is(vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  
  if (length (xyz.margin)<3) {
    warning ("xyz.margin must have a length of 3.")
    return (NULL)
  }
  if (!is.null(vol$local.gridx)) {
    vol$local.gridx <- NULL
    vol$local.gridy <- NULL
  }
  t.mat <- ref.cutplane.add(vol)
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo=paste0(vol$ref.pseudo,"m"), t.mat)
  
  
  vol.pt <- vol_$xyz.from.ijk %*% vol_$cube.idx[ , c (1,7)]
  pt.min <- c (min(vol.pt [1,]) - xyz.margin [1] , min(vol.pt [2,]) - xyz.margin [2], min(vol.pt [3,]) - xyz.margin [3])
  pt.max <- c (max(vol.pt [1,]) + xyz.margin [1] , max(vol.pt [2,]) + xyz.margin [2], max(vol.pt [3,]) + xyz.margin [3])
  if (pt.min[1]<=pt.max[1] & pt.min[2]<=pt.max[2] & pt.min[3]<=pt.max[3]){
    if (is.null(description)) 
      description<- paste (vol$description,"restricted to xyz.margin", paste (xyz.margin, collapse="|"))
    new.vol <- nesting.cube(vol_, pt.min, pt.max, alias, description )
    new.vol$object.alias <- vol$object.alias
    new.vol$object.info <- vol$object.info
    return(vol.in.new.ref(new.vol, new.ref.pseudo=vol$ref.pseudo, t.mat, 
                              alias = alias, description=description))

  }
  return (NULL)
}