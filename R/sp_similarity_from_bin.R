#' Volume-based spatial similarity metrics calculated from binary modality 3D volumes.
#' @description The \code{sp.similarity.from.bin} function computes volumetric Dice 
#' similarity coefficient, Dice-Jaccard coefficient and Dice surface similarity coefficient.
#' @param vol.A,vol.B "volume" class objects, of \code{"binary"} modality. \code{vol.B} is the reference for MDC calculation.
#' @param coeff Vector indicating the requested metrics from among 
#' 'DSC' (Dice similarity coefficient),'DJC' (Dice-Jaccard coefficient), 
#' and 'MDC' (mean distance to conformity). Equal to \code{NULL} if not requested.

#' @return returns a dataframe containing (if requested): 
#' \itemize{
#' \item volumetric Dice similarity coefficient \code{DSC} defined by :
#' \mjdeqn{DSC = 2 \frac{V_{A} ~\cap~ V_{B}}{V_{A} + V_{B}}}{ascii}
#' \item Dice-Jaccard coefficient \code{DJC} defined by :
#' \mjdeqn{DJC = \frac{V_{A} ~\cap~ V_{B}}{V_{A} ~\cup~ V_{B}}}{ascii}
#' \item mean distance to conformity \code{MDC}, over-contouring mean distance 
#' \code{over.MDC} and under-contouring mean distance \code{under.MDC}, defined by 
#' \emph{Jena et al} \strong{\[1\]}
#' }
#' 
#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{JENA201044}{espadon}

#' @seealso \link[espadon]{sp.similarity.from.mesh}
#' @examples
#' # creation of to volume" class objects, of "binary" modality
#' vol.A <- vol.create (pt000 = c(-25,-25,0), dxyz = c (1 , 1, 1),
#'                     n.ijk = c(50, 50, 1), default.value = FALSE,
#'                     ref.pseudo = "ref1",
#'                     alias = "vol.A", modality = "binary",
#'                     description = "") 
#' vol.B <- vol.copy (vol.A,alias = "vol.B")   
#' vol.A$vol3D.data [as.matrix(expand.grid(15:35,20:35,1))] <- TRUE
#' vol.A$max.pixel <- TRUE
#' vol.B$vol3D.data [as.matrix(expand.grid(16:36,18:37,1))] <- TRUE
#' vol.B$max.pixel <- TRUE
#' display.plane (vol.A, vol.B, interpolate = FALSE, 
#'                main = "vol.A & vol.B @ z = 0 mm") 
#'                
#' sp.similarity.from.bin (vol.A, vol.B)

#' @export
sp.similarity.from.bin <- function(vol.A, vol.B, 
                                   coeff = c('DSC', 'DJC', 'MDC', 'under.MDC', 'over.MDC')){
  if (!is (vol.A, "volume") | !is (vol.B, "volume")){
    warning ("vol.A or vol.B should be volume class objects.")
    return (NULL)
  }
  if ((vol.A$modality!="binary") | (vol.B$modality!="binary")) {
    warning ("both volumes must be of binary modality.")
    return (NULL)
  }
  if(is.null(vol.A$vol3D.data)){
    warning ("empty vol.A$vol3D.data.")
    return (NULL)
  }
  if(is.null(vol.B$vol3D.data)){
    warning ("empty vol.B$vol3D.data.")
    return (NULL)
  }
  #verifier que les volumes ont le même support
  if (!grid.equal (vol.A, vol.B)) {
    warning ("both volumes must share the same grid.")
    return (NULL)
  }
  
  f.idx <- match(coeff, c('DSC','DJC','MDC','under.MDC','over.MDC'))  
  f.idx <- f.idx[!is.na(f.idx)]
  label <- c()
  metrics <- NULL
  if (length(f.idx)!=0) { 
    label <- c('DSC', 'DJC', 'MDC', 'under.MDC', 'over.MDC')[f.idx]
    metrics <- matrix(0, nrow = 1, ncol =  length(label),  dimnames = list(NULL,label))
  }
  
  inter <- bin.intersection (vol.A,vol.B)
  if (is.na(inter$max.pixel)){
    metrics[] <- NA
    return(as.data.frame(metrics))
  }
  
  if (!(inter$max.pixel)) return(as.data.frame(metrics))
  
  metrics <- as.data.frame(metrics)
  inter.vol <- get.volume.from.bin (inter)
  vol.A_inter <- bin.subtraction (vol.A, inter)
  vol.B_inter <- bin.subtraction (vol.B, inter)
  
  vol.A.vol <- get.volume.from.bin (vol.A)
  vol.B.vol <- get.volume.from.bin (vol.B)
  if ('DSC' %in% label) metrics$DSC <- 2 * inter.vol / (vol.A.vol + vol.B.vol)
  if ('DJC' %in% label) metrics$DJC <- inter.vol/(vol.A.vol + vol.B.vol - inter.vol)
  
  
  if (any(!is.na(match(label,c('MDC', 'under.MDC', 'over.MDC'))))){
    L <- dist.to.conformity(vol.A_inter,vol.B_inter,inter)
    under.MDC <- ifelse(is.null(L$under.contouring),0, mean(L$under.contouring, na.rm = TRUE))
    over.MDC <- ifelse(is.null(L$over.contouring),0, mean(L$over.contouring, na.rm = TRUE))
    if ('MDC' %in% label) metrics$MDC <- under.MDC + over.MDC
    if ('under.MDC' %in% label) metrics$under.MDC <- under.MDC 
    if ('over.MDC' %in% label) metrics$over.MDC <- over.MDC
  }
  
 
  return(metrics)
}

