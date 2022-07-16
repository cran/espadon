#' Volume creating
#' @description The \code{vol.create} function creates a volume object from a 
#' user-defined grid.
#' @param n.ijk Vector of length 3, representing the number of elements on the i, 
#' j and k axes.
#' @param dxyz Vector of length 3, representing the x, y, z steps in mm, between 
#' voxels. See details.
#' @param mid.pt Vector of length 3, representing the x, y, z coordinates of the
#' midpoint  of the volume. See details.
#' @param pt000 Vector of length 3, representing the x, y, z coordinates of the
#' first voxel of the first plane.
#' @param default.value Numerical or boolean value, representing the default 
#' value of the voxels.
#' @param ref.pseudo Character string, frame of reference pseudonym of the 
#' created object.
#' @param frame.of.reference Character string, frame of reference of the 
#' created object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param modality Character string, \code{$modality} of the created object.
#' @param description Character string, describing the the created object.
#' @param number Integer, by default set to 0, number of the created object.
#' @details If \code{mid.pt} and \code{pt000} are both equal to \code{NULL}, 
#' then \code{mid.pt = c (0, 0, 0)} by default. 
#' If \code{mid.pt} and \code{pt000} are both different from \code{NULL}, then 
#' only \code{mid.pt} is taken into account. 


#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions), in which the grid is defined by  \code{pt000} or 
#' \code{mid.pt}, \code{dxyz} \code{n.ijk}. If \code{default.value} are 
#' initialized to \code{FALSE}, then \code{modality = "binary"}. 
#' The orientation of the patient is orthonormal to the grid.

#' @examples
#' new.vol <- vol.create (pt000 = c(1,10,10), dxyz = c (1 , 1, 1),
#'                        n.ijk = c(100, 100, 100),
#'                        ref.pseudo = "ref1",
#'                        frame.of.reference = "toyref1",
#'                        alias = "new ct", modality = "ct",
#'                        description = "")
#' str (new.vol)

#' @export
vol.create <- function (n.ijk, dxyz, mid.pt = NULL, pt000 = NULL, 
                        default.value = NA, ref.pseudo = "", 
                        frame.of.reference = "",
                        alias = "", modality = "",  
                        description = "", number = 0) {
  
  if (is.null(pt000) & is.null(mid.pt)) {
    mid.pt <- c (0, 0, 0)
    pt000 <- mid.pt - n.ijk*dxyz/2.0
  } else if(is.null(pt000)){
    pt000 <- mid.pt - n.ijk*dxyz/2.0
  }
  
  back.vol <- list()
  back.vol$patient <- ""
  back.vol$patient.bd <- ""
  back.vol$patient.sex <- ""
  
  back.vol$file.basename <- ""
  back.vol$file.dirname <- ""
  back.vol$object.name <- alias
  back.vol$object.alias <- alias
  

  back.vol$frame.of.reference <- frame.of.reference
  back.vol$ref.pseudo <- ref.pseudo
  back.vol$modality <- modality
  back.vol$description <- description
  if (!is.na(default.value) & modality == "" & is.logical(default.value)) back.vol$modality = "binary"
  
  
  back.vol$acq.date <- ""
  back.vol$study.date <- ""
  back.vol$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  # back.vol$error <- ""
  
  back.vol$number <- number

 
  back.vol$unit <- ""
  back.vol$n.ijk <- n.ijk
  back.vol$slice.thickness <- abs (dxyz[3])

  back.vol$min.pixel <- default.value
  back.vol$max.pixel <- default.value

  back.vol$dxyz<- dxyz
  back.vol$patient.orientation <- c(1,0,0,0,1,0)
  back.vol$patient.xyz0 <- matrix(pt000,ncol=3)
 

  
  
  back.vol$xyz.from.ijk <- as.matrix( cbind (c (1,0,0,0) * back.vol$dxyz[1],
                                             c (0,1,0,0) * back.vol$dxyz[2],
                                             c (0,0,1,0) * back.vol$dxyz[3],
                                             c (pt000,1)))
  
  back.vol$k.idx <- 0:(n.ijk[3]-1)
  back.vol$missing.k.idx <- FALSE
  
  back.vol$patient.xyz0  <- matrix((as.matrix (expand.grid (0, 0, back.vol$k.idx,1))%*% t(back.vol$xyz.from.ijk))[ ,1:3],ncol=3)
  
  back.vol$cube.idx <- matrix ( c(0,0,0,1,
                                  n.ijk[1]-1, 0, 0, 1,
                                  n.ijk[1]-1, n.ijk[2]-1, 0, 1,
                                  0, n.ijk[2]-1, 0, 1,
                                  0, 0, n.ijk[3]-1, 1,
                                  n.ijk[1]-1, 0, n.ijk[3]-1, 1,
                                  n.ijk[1]-1, n.ijk[2]-1, n.ijk[3]-1, 1,
                                  0, n.ijk[2]-1, n.ijk[3]-1, 1), nrow=4, byrow= FALSE)
  

  back.vol$vol3D.data <- array(default.value, dim=back.vol$n.ijk)

  class (back.vol) <- "volume"
  return(back.vol)
}