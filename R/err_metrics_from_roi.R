####################################################################
#' Error evaluation metrics for 2-volume comparisons 
#' @name err.metrics.from.roi
#' @aliases err.metrics.from.bin
#' @description The \code{err.metrics.from.roi} and \code{err.metrics.from.bin} 
#' calculate various metrics (ME, MAE, MSE, RMSE) to compare 2 “volume” class 
#' objects in the zones delimited by the requested RoI or binary selections
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
#' \item MIN: Min Error
#' \item MAX: Max Error
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
#' # Calculation using err.metrics.from.bin. The binary selection must first be 
#' # created.  
#' roi.idx <- select.names(patient$rtstruct[[1]]$roi.info$roi.pseudo,
#'                         roi.sname = c("eye","ptv"))  
#' bin.ROI <- lapply (roi.idx, function(idx){
#'    bin.from.roi (patient$ct[[1]], struct = patient$rtstruct[[1]], 
#'                  roi.idx = idx, T.MAT = patient$T.MAT, 
#'                  alias = patient$rtstruct[[1]]$roi.info$roi.pseudo[idx], 
#'                  description = patient$rtstruct[[1]]$roi.info$name[idx],
#'                  modality = "weight")})   
#'  names (bin.ROI) <-patient$rtstruct[[1]]$roi.info$name[roi.idx]   
#'  
#'  err.metrics.from.bin (obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
#'                        bin.list = bin.ROI, T.MAT= patient$T.MAT)                             
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
    roi.idx <- roi.idx[!is.na(struct$roi.info$vol[roi.idx])]
    roi.idx <- roi.idx[struct$roi.info$vol[roi.idx] > 0]
  }

  
  if (length(roi.idx)>0) {
    new.obj.ref <- nesting.roi(obj.ref, struct, roi.idx = roi.idx, 
                            xyz.margin =c(5,5,5), T.MAT = T.MAT)
  } else {
    # warning("No RoI available")
    new.obj.ref <- obj.ref 
    }

  
  if (length(roi.idx)>0){
  bin.ROI <- lapply(roi.idx, function(idx){
    bin.from.roi (new.obj.ref, struct = struct, roi.idx = idx, 
                  alias = struct$roi.info$roi.pseudo[idx], T.MAT = T.MAT, modality = "weight",
                  description = struct$roi.info$name[idx])})
  names(bin.ROI) <- struct$roi.info$name[roi.idx]
  } else {
  bin.ROI  <- list("all volume" =vol.copy(new.obj.ref, modality ="weight"))
  bin.ROI[["all volume"]]$vol3D.data[] <- bin.ROI[["all volume"]]$min.pixel <-bin.ROI[["all volume"]]$max.pixel <- 1
  bin.ROI[["all volume"]]$description <- "all volume"
  }
  .error.metrics (bin.ROI, obj, new.obj.ref, T.MAT)
}
.error.metrics <- function(sel.bin, obj, obj.ref, T.MAT){
  
  obj.to.test <- vol.regrid (obj, obj.ref, T.MAT = T.MAT) 
  
  error.vol <- vol.copy(obj.ref)
  error.vol$vol3D.data <- obj.to.test$vol3D.data - obj.ref$vol3D.data
  error.vol$max.pixel <- max(error.vol$vol3D.data, na.rm = TRUE)
  error.vol$min.pixel <- min(error.vol$vol3D.data, na.rm = TRUE)
  
  error.list <- lapply(sel.bin, function(bin){
    vol.from.bin(error.vol, bin, alias =  bin$object.alias, description =bin$object.alias )})
  names(error.list) <- names(sel.bin)
  metrics.list <- lapply((0:length(error.list))[-1], function(i){
    v <- as.numeric(error.list[[i]]$vol3D.data)
    w <- as.numeric(sel.bin[[i]]$vol3D.data)
    f <- !(is.na(v) | is.na(w))
    ntot <- sum(w[f])
    MSE <- sum( w[f] * (v[f]^2)) / ntot
    data.frame (ME = sum (w[f] * v[f]) / ntot,
                MAE = sum (w[f] * abs(v[f])) / ntot,
                MSE = MSE,
                RMSE = sqrt(MSE),
                MIN = min(v[f]),
                MAX = max(v[f]))
    })
  

  metrics.tab <- cbind(data.frame(ROI=names(sel.bin)),  
                       do.call(rbind, metrics.list))
   
  metrics.tab
}

#' @rdname err.metrics.from.roi
#' @param bin.list list of objects of class 'volume' and modality 'binary' or 'weight',
#' giving the selection of voxels in which metrics will be calculated.
#' @importFrom methods is
#' @export
err.metrics.from.bin <- function(obj, obj.ref, bin.list = NULL, T.MAT){
  
  if (!(is (obj, "volume") & is (obj.ref, "volume"))) 
    stop ("obj & obj.ref should be volume class objects.")
  if (is.null(bin.list)){
    bin.list  <- list("all volume" =vol.copy(obj.ref, modality ="weight"))
    bin.list[["all volume"]]$vol3D.data[] <- bin.list[["all volume"]]$min.pixel <-bin.list[["all volume"]]$max.pixel <- 1
    bin.list[["all volume"]]$description <- "all volume" 
  } else{
    
  if (!all(sapply (bin.list, function(l) is (l, "volume"))))
    stop ("bin.list must be a list of objects of class 'volume' and modality 'binary' or 'weight'.")
  if (!all(sapply (bin.list, function(l) l$modality %in% c("binary", "weight"))))
    stop ("bin.list must be a list of objects of class 'volume' and modality 'binary' or 'weight'.")
  
  }
  vgrid <- sapply (bin.list, function(l) !grid.equal(obj.ref,l))
  bin.list [vgrid]  <- lapply(bin.list [vgrid],function(l) vol.regrid(l,obj.ref,T.MAT = T.MAT))
  
  .error.metrics(bin.list, obj, obj.ref, T.MAT)

}