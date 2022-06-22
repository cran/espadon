#' Transform the grid of a volume class object into the grid of another
#' @description The \code{vol.regrid} function transforms the grid of a volume 
#' according to the grid of another.
#' @param vol "volume" class object to regrid.
#' @param back.vol "volume" class object whose grid will be used for regriding. 
#' Its \code{$ref.pseudo} must exist in the \code{T.MAT} list.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}. If \code{T.MAT = NULL}, 
#' \code{back.vol$ref.pseudo} must be equal to \code{vol$ref.pseudo}.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be that of \code{vol}.
#' @param verbose Boolean. if \code{TRUE} (default) a progress bar is displayed.
#' @return Returns a copy of \code{vol}, in which grid is that of \code{back.vol}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c ("mr", "rtdose"), 
#'                              dxyz = rep (step, 3), beam.nb = 4)
#' MR <- patient$mr[[1]]
#' D <- patient$rtdose[[1]]
#'
#' # change grid
#' D.on.MR <- vol.regrid (vol = D, back.vol = MR, interpolate = TRUE,
#'                        T.MAT = patient$T.MAT, alias = "",
#'                        description = NULL, verbose = FALSE)
#'
#' # maximum dose location
#' max.dose.in.MR <- get.xyz.from.index (which.max (D.on.MR$vol3D.data), D.on.MR)
#' display.plane (bottom = MR, view.coord = max.dose.in.MR[3],
#'                top= D.on.MR, bottom.col = grey.colors(255, start = 0, end = 1),
#'                bottom.breaks = seq (0, 500, length.out = 256),
#'                bg = "#00ff00", interpolate = FALSE)

#' @import progress
#' @importFrom methods is
#' @export
vol.regrid <- function (vol, back.vol, T.MAT = NULL, interpolate = TRUE, alias = "", 
                        description=NULL, verbose = TRUE){
  
  if (!is (vol, "volume")) {
    warning ( "vol should be a volume class object.")
    return (NULL)
  }
  
  if (!is (back.vol, "volume")){
    warning ("back.vol should be a volume class object.")
    return (NULL)
  }
  
  if (grid.equal (vol, back.vol)) return (vol.copy (vol, alias = alias, modality = vol$modality, description =description))
  M.Tref <- get.rigid.M (T.MAT, back.vol$ref.pseudo, vol$ref.pseudo)
  
  if (is.null (M.Tref)) {
    warning ("vol$ref.pseudo and back.vol$ref.pseudo are different.")
    return (NULL)
  }
  if (is.null(description)) description<- vol$description
  new.vol <- vol.copy (back.vol, alias = alias, modality = vol$modality, description =description)
  new.vol$patient <- vol$patient
  new.vol$patient.bd <- vol$patient.bd
  new.vol$patient.sex <- vol$patient.sex
  new.vol$acq.date <- vol$acq.date
  new.vol$study.date <- vol$study.date
  new.vol$creation.date <- vol$creation.date
  new.vol$study.time <- vol$study.time
  new.vol$number <- vol$number
  new.vol$unit <- vol$unit
  
  new.vol$vol3D.data[] <- NA
  
  if (is.na(vol$min.pixel) | is.na(vol$max.pixel)) return (new.vol)
  if (abs (det(vol$xyz.from.ijk))<1e-6) {
    warning ("non invertible matrix, check arguments.")
    return (NULL)
  }
  M.old.from.new <- solve (vol$xyz.from.ijk) %*% M.Tref %*% new.vol$xyz.from.ijk
  if (abs (det(M.old.from.new))<1e-6) {
    warning ("non invertible matrix, check arguments.")
    return (NULL)
  }
  new.cube.idx <- solve (M.old.from.new) %*% vol$cube.idx
  rg.i <- floor (min(new.cube.idx[1,])):ceiling(max(new.cube.idx[1,]))
  rg.j <- floor (min(new.cube.idx[2,])):ceiling(max(new.cube.idx[2,]))
  rg.k <- floor (min(new.cube.idx[3,])):ceiling(max(new.cube.idx[3,]))
  rg.i <- rg.i[!is.na (match (rg.i,(1:new.vol$n.ijk[1])-1))]
  rg.j <- rg.j[!is.na (match (rg.j,(1:new.vol$n.ijk[2])-1))]
  rg.k <- rg.k[!is.na (match (rg.k, new.vol$k.idx))]
  n <- ceiling ((2^18) / (length(rg.i) * length(rg.j)))
  i <- 0
  N <- ceiling(length(rg.k)/n)
# print(N)
  if (verbose) pb <- progress_bar$new(format = " processing [:bar] :percent",
                                      total = N, clear = FALSE, width= 60)
  repeat {
    k.idx <- rg.k[i + (1:n)]
    k.idx <- k.idx[!is.na(k.idx)]
 
    if (length (k.idx)==0) break
    ijk.selection <- as.matrix (expand.grid (rg.i, rg.j, k.idx, 1))
    ijk.Rselection <- ijk.selection[ , 1:3] +1
    ijk <-(ijk.selection %*% t(M.old.from.new))[ ,1:3]
    new.vol$vol3D.data[ijk.Rselection] <- get.value.from.ijk (ijk, vol, interpolate)
    i <- i+n
    if (verbose) pb$tick()
  }
  #new.vol$vol3D.data[new.vol$vol3D.data<=1e-6] <- 0
  if (any(!is.na(new.vol$vol3D.data))){
    new.vol$min.pixel <- min ( new.vol$vol3D.data, na.rm=TRUE)
    new.vol$max.pixel <- max ( new.vol$vol3D.data, na.rm=TRUE)
  } else {
    new.vol$min.pixel <- NA
    new.vol$max.pixel <- NA
  }
  return (new.vol)
  
}
