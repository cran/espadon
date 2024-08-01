#' Dosimetry, volume, conformity, homogeneity indices from binary selection
#' @description The \code{rt.indices.from.bin} function calculates, from a 
#' "volume" class object of modality "rtdose", all the standard 
#' indicators of radiotherapy, as long as their options are transmitted, for the 
#' target and healthy "volume" object of modality "binary".

#' @param vol "volume" class object, of "rtdose" modality.
#' @param target.bin.list list of "volume" class objects, of "binary" 
#' modality. The \code{$object.alias} field of each \code{target.bin.list} object 
#' represents the name of the selected region of interest of the target volume.
#' @param healthy.bin.list list of "volume" class objects, of "binary" 
#' modality. The \code{$object.alias} field of each \code{healthy.bin.list} object 
#' represents the name of the selected region of interest of the healthy tissues.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, all \code{$ref.pseudo} 
#' of \code{bin.list} elements must be equal to \code{vol$ref.pseudo}.
#' 
#' @param presc.dose vector of prescription doses that serve as reference doses 
#' for the target RoI.
#' @param healthy.tol.dose vector of tolerance dose of each healthy RoI. 
#' @param healthy.weight Vector of weight, indicating the importance of the healthy 
#' RoI.
  
#' @param dosimetry Vector indicating the requested dose indicators from among 
#' 'D.min', 'D.max', 'D.mean' and 'STD'. If \code{D.xpc} is different from 
#' \code{NULL}, it will be added.
#' @param volume.indices Vector indicating the requested volume indices from among 
#' 'V.tot',  'V.prescdose' (i.e. volume over \code{presc.dose}) and 'area'. If 
#' \code{V.xGy} is different from \code{NULL}, it will be added.
#' @param conformity.indices Vector. Requested conformity indices from among 'PITV',
#' 'PDS', 'CI.lomax2003', 'CN', 'NCI', 'DSC', 'CI.distance', 'CI.abs_distance', 
#' 'CDI', 'CS3', 'ULF', 'OHTF', 'gCI', 'COIN', 'COSI' and 'G_COSI'.
#' @param homogeneity.indices Vector. Requested homogeneity indices from among 
#' 'HI.RTOG.max_ref', 'HI.RTOG.5_95', 'HI.ICRU.max_min', 'HI.ICRU.2.98_ref', 
#' 'HI.ICRU.2.98_50', 'HI.ICRU.5.95_ref', 'HI.mayo2010' and 'HI.heufelder.'
#' @param gradient.indices Vector. Requested gradient indices from among 
#' 'GI.ratio.50', 'mGI'.

#' @param D.xpc Vector of the percentage of the volume, for which the dose coverage 
#' is requested.
#' @param D.xcc Vector of the volume in \mjeqn{cm^3}{ascii}, for which the dose 
#' coverage is requested.
#' @param V.xpc Vector of the percentage of the reference dose, received by the volume to be 
#' calculated.
#' @param V.xGy Vector of the minimum dose in Gy, received by the volume to be 
#' calculated.
#' @param DVH boolean. If \code{TRUE} (default), the function returns the DVHs of 
#' the target and healthy ROIs.
#' @param verbose Boolean. if \code{TRUE} (default) a progress bar is displayed.
#' @param ... others parameters such as \code{DVH.step}.
#' @return Return  a list of  indices dataframe. For details, see 
#' \link[espadon]{rt.indices.from.roi}.
#' @seealso \link[espadon]{rt.indices.from.roi}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for better
#' # result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), roi.name = "eye",
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]
#' struct <- patient$rtstruct[[1]]
#' T.MAT <- patient$T.MAT
#' 
#' # creation of the list of target binary volumes
#' target.roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.sname = "ptv")
#' healthy.roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.sname = "eye")
#' target.bin.list <- lapply (target.roi.idx , function (idx) {
#'   vr <- nesting.roi (D, struct, roi.idx = idx, xyz.margin = c (5, 5, 5),
#'                      T.MAT = T.MAT, alias = struct$roi.info$name[idx])
#'   b <- bin.from.roi(vr, struct, roi.idx = idx, T.MAT = T.MAT,
#'                     alias = struct$roi.info$name[idx], verbose = FALSE)
#'   })
#' names (target.bin.list) <- struct$roi.info$name[target.roi.idx]
#' 
#' healthy.bin.list <- lapply (healthy.roi.idx , function (idx) {
#'   vr <- nesting.roi (D, struct, roi.idx = idx, xyz.margin = c (5, 5, 5),
#'                      T.MAT = T.MAT, alias = struct$roi.info$name[idx])
#'   b <- bin.from.roi(vr, struct, roi.idx = idx, T.MAT = T.MAT,
#'                     alias = struct$roi.info$name[idx], verbose = FALSE)
#' })
#' names (healthy.bin.list) <- struct$roi.info$name[healthy.roi.idx]
#' 
#' indices <- rt.indices.from.bin (D, target.bin.list, healthy.bin.list,
#'                                 presc.dose = 50,
#'                                 conformity.indices = c("PITV", "PDS", "CI.lomax2003", 
#'                                                        "CN", "NCI", "DSC","COIN"),
#'                                 verbose = FALSE)
#' indices[c("dosimetry","volume", "conformity","homogeneity","gradient")]
#' head(indices$DVH)

#' @import progress
#' @importFrom methods is
#' @export

rt.indices.from.bin <- function (vol, 
                                 target.bin.list = NULL, 
                                 healthy.bin.list = NULL, 
                                 T.MAT = NULL, 
                                 presc.dose = NA,
                                 healthy.tol.dose = NA,
                                 healthy.weight = 1,
                                 dosimetry = c("D.min", "D.max", "D.mean", "STD"),
                                 volume.indices = c("V.tot", "area", "V.prescdose"),
                                 conformity.indices =  c("PITV", "CI.lomax2003", "CN",
                                                         "NCI", "DSC", "CI.distance", 
                                                         "CI.abs_distance", "CDI", "CS3", 
                                                         "ULF","OHTF", "gCI", 
                                                         "COIN", "G_COSI", "COSI"),
                                 homogeneity.indices = c("HI.RTOG.max_ref", "HI.RTOG.5_95", 
                                                         "HI.ICRU.max_min", 
                                                         "HI.ICRU.2.98_ref", "HI.ICRU.2.98_50", 
                                                         "HI.ICRU.5.95_ref", "HI.mayo2010", 
                                                         "HI.heufelder"),
                                 gradient.indices = c("GI.ratio.50"),
                                 D.xpc = NULL, D.xcc = NULL, V.xpc = NULL, V.xGy = NULL, 
                                 DVH = TRUE,
                                 verbose = TRUE,...){
  args <- list(...)
  DVH.step <- args[['DVH.step']]
  if (is.null(DVH.step)) DVH.step <-0.001
  if (DVH.step>0.01) DVH.step <- 0.01
  
  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  if (vol$modality !="rtdose") warning ("vol should be of rtdose modality.")
  
  
  
  if (!is.null(target.bin.list)){
    class.target <- unique(sapply(target.bin.list, class))
    if (length(class.target)!=1 | class.target!="volume") 
      stop ("target.bin.list should be a list of volume class objects.")
    
    modality.target <-  unique(sapply(target.bin.list, function(l) l$modality))
    if (length(modality.target)!=1 | modality.target!="binary") 
    stop ("target.bin.list should be a list of objects of binary modality.")
  }
  if (!is.null(healthy.bin.list)){
    class.target <- unique(sapply(healthy.bin.list, class))
    if (length(class.target)!=1 | class.target!="volume") 
      stop ("healthy.bin.list should be a list of volume class objects.")
    
    modality.target <-  unique(sapply(healthy.bin.list, function(l) l$modality))
    if (length(modality.target)!=1 | modality.target!="binary") 
      stop ("healthy.bin.list should be a list of objects of binary modality.")
  }
  vol.ref <-  vol$ref.pseudo
  
  bin.L_ <- c(target.bin.list, healthy.bin.list)
  bin.L <- lapply(bin.L_, function(l) vol.in.new.ref (l, vol.ref, T.MAT, alias=l$object.alias))
  bin.L <- bin.L[!sapply(bin.L, is.null)]
  if (length(bin.L)==0) bin.L <-NULL
  
  vol.L <- NULL
  vol.names <- NULL
  L.target.idx <- NULL
  L.healthy.idx <- NULL

  if (!is.null(bin.L)){ 
    vol.names<- as.character((trimws(sapply(bin.L, function(l) l$object.alias))))
    names(bin.L) <- vol.names
    
    target.names <- as.character((trimws(sapply(target.bin.list, function(l) l$object.alias))))
    healthy.names <- as.character((trimws(sapply(healthy.bin.list, function(l) l$object.alias))))
    L.target.idx <- match (target.names,vol.names)
    L.target.idx <- L.target.idx[!is.na(L.target.idx)]
    if (length(L.target.idx)==0) L.target.idx <- NULL
    L.healthy.idx <- match (healthy.names,vol.names)
    L.healthy.idx <- L.healthy.idx[!is.na(L.healthy.idx)]
    if (length(L.healthy.idx)==0) L.healthy.idx <- NULL
    
    if (verbose) pb <- progress_bar$new(format = " processing [:bar] :percent",
                                        total = length(bin.L), width= 60)
    vol.L <- list()
    for (l.idx in 1:length(bin.L)) {
      l <- bin.L[[l.idx]]
      vol_ <- vol
      if (!grid.equal(vol_,l)) vol_ <- vol.regrid(vol,l,T.MAT=T.MAT, verbose = FALSE)
      vol.L[[l.idx]] <- vol.from.bin(vol_,l, alias = l$object.alias, description=l$description)
      if (verbose) pb$tick()
    }
    names(vol.L) <- vol.names
  }
  if (is.null(presc.dose)) presc.dose <- NA
  if (is.null(healthy.tol.dose)) healthy.tol.dose <- NA
  
  .indices.from.bin(vol,presc.dose, bin.L, vol.L, vol.names,
                    L.target.idx, L.healthy.idx,
                    dosimetry,
                    volume.indices,
                    conformity.indices,
                    homogeneity.indices,
                    gradient.indices,
                    D.xpc, D.xcc, V.xpc, V.xGy, 
                    healthy.weight, healthy.tol.dose,DVH,DVH.step)
  
}