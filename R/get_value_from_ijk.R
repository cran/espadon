#' Value of the volume at a selection of DICOM indices
#' @description The \code{get.value.from.ijk} function calculates the value of 
#' a "volume" class object at DICOM indices i, j, k, whether they are 
#' integers or not.
#' @param ijk Vector or 3-column matrix of DICOM indices.
#' @param vol "volume" class object.
#' @param interpolate Boolean, default to \code{TRUE}. If \code{interpolate = TRUE}, a 
#' trilinear interpolation of the value of the voxels, relative to the values of 
#' adjacent voxels, is performed.
#' @return Returns a vector of the values of the volume at the requested DICOM 
#' indices.
#' @seealso \link[espadon]{get.ijk.from.index}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#' # isodose
#' Dmax <- max (D$vol3D.data, na.rm = TRUE)
#' Dmax
#' idx <- which (D$vol3D.data >= (Dmax -1) & D$vol3D.data <= (Dmax - 0.2))
#' ijk <- get.ijk.from.index (idx, D)
#' get.value.from.ijk (ijk, vol = D, interpolate = FALSE)

#' @export
#' @importFrom methods is
get.value.from.ijk <- function (ijk, vol, interpolate = TRUE)  {
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  toC3M <- function (vect) {return(matrix(vect,ncol=3))}
  ijk <- toC3M (ijk)

  k.idx<- match(0:max(vol$k.idx),vol$k.idx)
  k.loc <- k.idx-1
  fna <- is.na(k.idx)
  k.idx[!fna] <- vol$k.idx
  k.loc[fna] <- max(vol$k.idx)+1
  k.idx[fna] <- max(vol$k.idx)+1
  
  value <- .getvaluefromijkC (vol3D = as.numeric(vol$vol3D.data),
                              interpolate = interpolate,
                              i = as.numeric(ijk[ ,1]),
                              j = as.numeric(ijk[ ,2]),
                              k = as.numeric(ijk[ ,3]),
                              k_idx = k.idx,
                              k_loc = k.loc, n_ijk=vol$n.ijk)
  
  value [is.nan(value)] <- NA
  return(value)


  # value <- rep(NA,nrow(ijk))
  # if (!interpolate) {
  #   round.idx <- floor (ijk + 0.5) + 1
  # 
  #   keep <- !is.na (match (round.idx[,1],(1:vol$n.ijk[1]))) &
  #     !is.na (match (round.idx[,2],(1:vol$n.ijk[2]))) &
  #     !is.na (match (round.idx[,3],vol$k.idx + 1))
  #   keep.map.ind <- toC3M(round.idx[keep,])
  #   if (vol$missing.k.idx) keep.map.ind[,3] <- (1:vol$n.ijk[3])[match(keep.map.ind[,3],vol$k.idx + 1)]
  # 
  #   value[keep] <- vol$vol3D.data[keep.map.ind]
  # 
  # } else {
  #   idx000 <- floor (ijk) + 1
  #   idx111 <- floor (ijk) + 2
  #   keep <- rep(TRUE, nrow(ijk))
  #   if (vol$n.ijk[1]>1) {
  #     keep <- keep & !is.na (match (idx000[,1],(1:vol$n.ijk[1]))) &
  #       !is.na (match (idx111[,1],(1:vol$n.ijk[1])))
  #   } else {
  #     keep <- keep & (idx000[,1]==1)
  #   }
  #   if (vol$n.ijk[2]>1) {
  #     keep <- keep & !is.na (match (idx000[,2],(1:vol$n.ijk[2]))) &
  #       !is.na (match (idx111[,2],(1:vol$n.ijk[2])))
  #   } else {
  #     keep <- keep & (idx000[,2]==1)
  #   }
  #   if (vol$n.ijk[3]>1) {
  #     keep <- keep & !is.na (match (idx000[,3],vol$k.idx + 1)) &
  #       !is.na (match (idx111[,3],vol$k.idx + 1))
  #   } else {
  #     keep <- keep & (idx000[,3]== vol$k.idx + 1)
  #   }
  # 
  #   idx000.keep <- toC3M(idx000[keep,])
  #   idx111.keep <- toC3M(idx111[keep,])
  #   ijk.keep <- toC3M(ijk[keep,]) + 1
  # 
  #   muvw <- 1 + idx000.keep - ijk.keep
  #   puvw <- 1 - idx111.keep + ijk.keep
  # 
  #   if (vol$missing.k.idx) {
  #     idx000.keep[,3] <- (1:vol$n.ijk[3])[match(idx000.keep[,3],vol$k.idx + 1)]
  #     idx111.keep[,3] <- (1:vol$n.ijk[3])[match(idx111.keep[,3],vol$k.idx + 1)]
  #     # ijk.keep[,3] <- (1:vol$n.ijk[3])[match(ijk.keep[,3],vol$k.idx + 1)]
  #   }
  # 
  # 
  #   value[keep] <- vol$vol3D.data [idx000.keep] * muvw[ , 1] * muvw[ , 2] * muvw[ , 3]
  #   if (vol$n.ijk[3]>1) {
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx000.keep[, 1], idx000.keep[, 2], idx111.keep[, 3])] * muvw[ , 1] * muvw[ , 2] * puvw[ , 3]
  #   }
  #   
  #   if (vol$n.ijk[2]>1) {
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx000.keep[, 1], idx111.keep[, 2], idx000.keep[, 3])] * muvw[ , 1] * puvw[ , 2] * muvw[ , 3]
  #   }
  #   
  #   if ((vol$n.ijk[2]>1) & (vol$n.ijk[3]>1) ){
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx000.keep[, 1], idx111.keep[, 2], idx111.keep[, 3])] * muvw[ , 1] * puvw[ , 2] * puvw[ , 3]
  #   } 
  #   
  #   if (vol$n.ijk[1]>1) {
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx111.keep[, 1], idx000.keep[, 2], idx000.keep[, 3])] * puvw[ , 1] * muvw[ , 2] * muvw[ , 3]
  #   } 
  #   
  #   if ((vol$n.ijk[1]>1) & (vol$n.ijk[3]>1) ){
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx111.keep[, 1], idx000.keep[, 2], idx111.keep[, 3])] * puvw[ , 1] * muvw[ , 2] * puvw[ , 3]
  #   } 
  #   
  #   if ((vol$n.ijk[1]>1) & (vol$n.ijk[2]>1) ){
  #     value[keep] <- value[keep] + vol$vol3D.data [cbind (idx111.keep[, 1], idx111.keep[, 2], idx000.keep[, 3])] * puvw[ , 1] * puvw[ , 2] * muvw[ , 3]
  #   } 
  #   
  #   if ((vol$n.ijk[1]>1) & (vol$n.ijk[2]>1) & (vol$n.ijk[3]>1) ){
  #     value[keep] <- value[keep] + vol$vol3D.data [idx111.keep] * puvw[ , 1] * puvw[ , 2] * puvw[ , 3]
  #   } 
  # }
  # return(value)
}