####################################################################
#' Error evaluation metrics for 2-volume comparisons 
#' @description The \code{err.metrics.from.roi} calculates various metrics 
#' (ME, MAE, MSE, RMSE) to compare 2 “volume” class objects in the zones delimited 
#' by the requested RoI.
#' @param obj "volume" class object to be compared.
#' @param obj.ref "volume" class reference object.
#' @param struct "struct" class object or NULL.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} object.
#' By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object.
#' By default \code{roi.idx = NULL}. See Details.
#' @param T.MAT "t.mat" class object to link the reference frames of \code{obj} 
#' \code{obj.ref} and \code{struct}. \code{T.MAT} can be created by 
#' \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If 
#' \code{T.MAT = NULL}, \code{struct$ref.pseudo} must be equal to 
#' \code{obj$ref.pseudo} and \code{obj.ref$ref.pseudo}.

#' @return Returns, in the zones delimited by the requested RoI, the following metrics:
#' \itemize{
#' \item ME: Mean Error
#' \item MAE: Mean Absolute Error
#' \item MSE: Mean Squared Error
#' \item RMSE: Root Mean Squared Error
#' }
#' 
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are
#' all set to \code{NULL}, all RoI whose volume is greater than 0 are selected.

#' @examples
#' # loading of toy-patient objects (decrease dxyz)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "sct","rtstruct"), 
#'                              roi.name = c("eye", "brain","gizzard"),
#'                              dxyz = rep (step, 3))
#' 
#' patient$ct[[1]]$description
#' patient$ct[[2]]$description

#' # Calculation of eye zone and ptv metrics
#' err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
#'                      struct = patient$rtstruct[[1]], roi.sname = c("eye","ptv"),
#'                      T.MAT= patient$T.MAT)
#' # Calculation of metrics on delimited zones on all RoIs
#' err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
#'                      struct = patient$rtstruct[[1]],
#'                      T.MAT= patient$T.MAT)
#'                      
#' # Calculation on all volume
#' err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
#'                      T.MAT= patient$T.MAT)
#'                      
#' @importFrom methods is
#' @export
err.metrics.from.roi <- function(obj, obj.ref, struct = NULL, roi.name = NULL, roi.sname = NULL, 
                      roi.idx = NULL, T.MAT = NULL){
  if (!(is (obj, "volume") & is (obj.ref, "volume"))) 
    stop ("obj & obj.ref should be volume class objects.")
  if (is.null(struct)){
    roi.idx <- integer(0)
  } else {
    if (!is (struct, "struct")) stop ("struct should be a struct class object.")
    roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
    roi.idx <- roi.idx[struct$roi.info$vol[roi.idx] >0]
  }

  
  if (length(roi.idx)>0) {
    new.obj1 <- nesting.roi(obj, struct, roi.idx = roi.idx, 
                            xyz.margin =c(5,5,5), T.MAT = T.MAT)
  } else {
    warning("No RoI available")
    new.obj1 <- obj 
    }
  new.obj2 <- vol.regrid (obj.ref, new.obj1, T.MAT = T.MAT) 
  
  error.vol <- vol.copy(new.obj1)
  error.vol$vol3D.data <- new.obj1$vol3D.data - new.obj2$vol3D.data
  error.vol$max.pixel <- max(error.vol$vol3D.data, na.rm = TRUE)
  error.vol$min.pixel <- min(error.vol$vol3D.data, na.rm = TRUE)
  
  if (length(roi.idx)>0){
  bin.ROI <- lapply(roi.idx, function(idx){
    bin.from.roi (new.obj1, struct = struct, roi.idx = idx, 
                  alias = struct$roi.info$roi.pseudo[idx], T.MAT = T.MAT, modality = "weight")})
  names(bin.ROI) <- roi.idx
  } else {
  bin.ROI  <- list(weight =vol.copy(error.vol, modality ="weight"))
  bin.ROI[[1]]$vol3D.data[] <- bin.ROI[[1]]$min.pixel <-bin.ROI[[1]]$max.pixel <- 1
  }
  
  error.list <- lapply(bin.ROI, function(bin){
    vol.from.bin(error.vol, bin, alias =  bin$object.alias, description =bin$object.alias )})
  names(error.list) <- roi.idx
  metrics.list <- lapply((0:length(error.list))[-1], function(i){
    v <- as.numeric(error.list[[i]]$vol3D.data)
    w <- as.numeric(bin.ROI[[i]]$vol3D.data)
    f <- !(is.na(v) | is.na(w))
    ntot <- sum(w[f])
    MSE <- sum( w[f] * (v[f]^2)) / ntot
    data.frame (ME = sum (w[f] * v[f]) / ntot,
                MAE = sum (w[f] * abs(v[f])) / ntot,
                MSE = MSE,
                RMSE = sqrt(MSE))
    })
  
   if (length(roi.idx)>0) {
     metrics.tab <- cbind(data.frame(ROI=struct$roi.info$name[roi.idx]),  
                       do.call(rbind, metrics.list))
   } else {
     metrics.tab <- cbind(data.frame(ROI="all volume"),  
                          do.call(rbind, metrics.list))
     }
  metrics.tab
  
}
