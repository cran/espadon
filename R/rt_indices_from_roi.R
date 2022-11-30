#' Dosimetry, volume, conformity, homogeneity indices from RoI
#' \loadmathjax
#' @description The \code{rt.indices.from.roi} function calculates, from a "volume"
#' class object of modality "rtdose", standard indicators of radiotherapy
#' in relation to the target and healthy RoI, as long as their options are transmitted.

#' @param vol "volume" class object, of "rtdose" modality.
#' @param struct "struct" class object.
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} 
#' or  \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, \code{struct$ref.pseudo}
#' must be equal to \code{vol$ref.pseudo}.
#' 
#' @param target.roi.name Exact name of target RoI in \code{struct} object. 
#' By default \code{target.roi.name = NULL}. See Details.
#' @param target.roi.sname Name or part of name of target RoI in \code{struct} 
#' object. By default \code{target.roi.sname = NULL}. See Details.
#' @param target.roi.idx Value of the index of target RoI that belong to the 
#' \code{struct} object. By default \code{target.roi.idx = NULL}. See Details.
#' @param healthy.roi.name Exact name of healthy RoI in \code{struct} object. 
#' By default \code{healthy.roi.name = NULL}. 
#' @param healthy.roi.sname Name or part of name of healthy RoI in \code{struct} 
#' object. By default \code{healthy.roi.sname = NULL}.
#' @param healthy.roi.idx Value of the index of healthy RoI that belong to the 
#' \code{struct} object. By default \code{healthy.roi.idx = NULL}.
#' 
#' @param presc.dose Vector of prescription doses that serve as reference doses 
#' for the target RoI.
#' @param healthy.tol.dose Vector of tolerance doses of each healthy RoI. 
#' @param healthy.weight Vector of weights, indicating the importance of the 
#' healthy RoI.

#' @param dosimetry Vector indicating the requested dose indicators from among 
#' 'D.min', 'D.max', 'D.mean' and 'STD.' If \code{D.xpc} is different from 
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
#' @param V.xGy Vector of the minimum dose in Gy, received by the volume to be 
#' calculated.
#' 
#' @param verbose Boolean. if \code{TRUE} (default) a progress bar is displayed.

#' @details If \code{target.roi.name}, \code{target.roi.sname}, and 
#' \code{target.roi.idx} are all set to \code{NULL}, all RoI containing 'ptv' 
#' (if they exist) are selected.

#' @return  Return  a list containing (if requested)

#' @return \mjeqn{-~dosimetry}{ascii} : dataframe containing, for all target and 
#' healthy structures:
#' \itemize{
#' \item the requested \code{dosimetry} : \code{D.min} (Gy), \code{D.max} (Gy), 
#' \code{D.mean} (Gy) and \code{STD} (Gy), respectively the minimum, maximum, 
#' mean and standard deviation of the dose in the regions of interest.
#' \item the requested \code{$D.x\%} : (Gy) Dose covering x percent of structure volume.
#' \item the requested \code{$D.xcc} : (Gy) Dose covering x (\mjeqn{cm^3}{ascii})
#'  of structure volume.
#' } 

#' @return \mjeqn{-~volume}{ascii} : dataframe containing, for all target and 
#' healthy structures, and isodoses:
#' \itemize{
#' \item the requested \code{volume.indices} : \code{V_tot} (\mjeqn{cm^3}{ascii}) 
#' (except for isodose) the total volume of the regions of interest, \code{area} 
#' (\mjeqn{cm^2}{ascii}) (except for isodose) their surface areas, 
#' \code{V.prescdose} (\mjeqn{cm^3}{ascii}) the volumes  receiving at least 
#' \code{presc.dose} Gy,
#' \item the requested \code{V.xGy} (\mjeqn{cm^3}{ascii}): 
#' volumes receiving at least x Gy.
#' }
#' 
#' @return \mjeqn{-~conformity}{ascii} : dataframe containing, if requested,
#' \itemize{
#' \item \code{PITV} : (1) Prescription Isodose Target Volume, or conformity index
#'  defined by \emph{E.Shaw} \strong{\[1\]}
#' \mjdeqn{PITV = \frac{V_{presc.dose}}{V_{target}}}{ascii}

#' \item \code{PDS} : (1) Prescription Dose Spillage
#'  defined by \emph{SABR UK Consortium 2019} \strong{\[2\]}
#' \mjdeqn{PDS = \frac{V_{presc.dose}}{V_{target ~\ge~ presc.dose}} = 
#' \frac{V_{presc.dose}}{V_{target} ~\cap~ V_{presc.dose}}}{ascii}


#' \item \code{CI.lomax2003} : (1) Conformity Index defined by \emph{Lomax and al} 
#' \strong{\[3\]}
#' \mjdeqn{CI_{lomax2003} = \frac{V_{target ~\ge~ presc.dose}}{V_{presc.dose}} = 
#' \frac{V_{target} ~\cap~ V_{presc.dose}}{V_{presc.dose}}}{ascii}

#' \item \code{CN} : (1) Conformation Number defined by \emph{Van't Riet and al} 
#' \strong{\[4\]}. It corresponds to conformity index defined by \emph{Paddick} 
#' \strong{\[5\]}
#' \mjdeqn{CN = CI_{paddick2000} =\frac{V_{target ~\ge~ presc.dose}^2}{V_{target}~\cdot~V_{presc.dose}} = 
#' \frac{(V_{target} ~\cap~ V_{presc.dose})^2}{V_{target}~\cdot~V_{presc.dose}}}{ascii}
 
#' \item \code{NCI} : (1) New conformity index, inverse of CN, defined by 
#' \emph{Nakamura and al} \strong{\[6\]}
#' \mjdeqn{NCI  =\frac{1}{CN}}{ascii}

#' \item \code{DSC} : (1) Dice Similarity Coefficient \strong{\[7\]}
#' \mjdeqn{DSC =  2 ~\cdot~\frac{V_{target ~\ge~ presc.dose}}{V_{target} + V_{presc.dose}} =
#' 2 ~\cdot~\frac{V_{target} ~\cap~ V_{presc.dose}}{V_{target} + V_{presc.dose}}}{ascii}
                       

#' \item \code{CI.distance} : (1) Conformity Index based on distance defined by 
#' \emph{Park and al} \strong{\[8\]}
#' \mjdeqn{CI.distance = \frac{100}{N} \sum^N \frac{dist_{S_{presc.dose}~\to~G_{target}} - 
#' dist_{S_{target}~\to~G_{target}}}{dist_{S_{target}~\to~G_{target}}}}{ascii}
#' where \mjeqn{dist_{S_{presc.dose}~\to~G_{target}}}{ascii} is the distance between
#' the surface of the prescription dose volume and the centroid of the target, 
#' and \mjeqn{dist_{S_{target}~\to~G_{target}}}{ascii} the surface of the target 
#' volume and the centroid of the target. 
#' \mjeqn{N}{ascii} is the number of directions where the distances are calculated. 
#' These directions are computed every 1°. If the centroid is not within the target 
#' volume, then \code{CI.distance = NA}.

#' \item \code{CI.abs_distance} : (1) Conformity Index based on distance defined 
#' by \emph{Park and al} \strong{\[8\]}
#' \mjdeqn{CI.abs_distance = \frac{100}{N} \sum^N \frac{|dist_{S_{presc.dose}~\to~G_{target}} - 
#' dist_{S_{target}~\to~G_{target}}|}{dist_{S_{target}~\to~G_{target}}}}{ascii}
# 
#' \item \code{CDI} : (1) Conformity Distance Index defined by \emph{Wu and al} 
#' \strong{\[9\]}
#' \mjdeqn{CDI = 2 \frac{V_{presc.dose} + V_{target} - 2~V_{target ~\ge~ presc.dose}}
#' {S_{target} + S_{presc.dose}} = \frac{V_{presc.dose} + V_{target} - 2~\cdot~V_{target} ~\cap~ V_{presc.dose}}
#' {S_{target} + S_{presc.dose}}}{ascii}
#' where \mjeqn{S_{target}}{ascii} is the surface of the target volume and 
#' \mjeqn{S_{presc.dose}}{ascii} is the surface of the prescription dose volume.
# 
#' \item \code{CS3} : (1) Triple Point Conformity Scale defined by \emph{Ansari 
#' and al} \strong{\[10\]}
#' \mjdeqn{CS3 = \frac{V_{0.95~\cdot~presc.dose} +  V_{presc.dose} +
#'  V_{1.05~\cdot~presc.dose}}{3~\cdot~V_{target}}}{ascii}

#' \item \code{ULF} : (1) Underdosed lesion factor defined by \emph{Lefkopoulos 
#' and al} \strong{\[11\]}
#' \mjdeqn{ULF = \frac{V_{target ~<~ presc.dose}}{V_{target}}=
#'  \frac{V_{target} ~\cap~ \overline{V_{presc.dose}}}{V_{target}}}{ascii}

#' \item \code{OHTF} :(1) Overdosed healthy tissues factor defined by \emph{Lefkopoulos 
#' and al} \strong{\[11\]}
#' \mjdeqn{OHTF = \frac{\sum V_{healthy ~\ge~ presc.dose}}{V_{target}} = 
#' \frac{\sum V_{healthy} ~\cap~ V_{presc.dose}}{V_{target}} }{ascii}

#' \item \code{gCI} : (1) Geometric Conformity Index  defined by 
#' \emph{Lefkopoulos and al} \strong{\[11\]}
#' \mjdeqn{gCI = ULF + OHTF}{ascii}

#' \item \code{COIN} : Conformity Index defined by \emph{Baltas and al} \strong{\[12\]}
#' \mjdeqn{COIN  =\frac{V_{target ~\ge~ presc.dose}^2}{V_{target}~\cdot~V_{presc.dose}}~\cdot~
#'   \prod^{N_{healthy}} \left( 1 -\frac{V_{healthy ~\ge~ presc.dose}}{V_{healthy}}\right)}{ascii}

#' \item \code{gCOSI} : generalized COSI, defined by \emph{Menhel and al} \strong{\[13\]}.
#' \mjdeqn{gCOSI  = 1- \sum^{N_{healthy}}  healthy.weight~\cdot~
#' \frac{\frac{V_{healthy ~\ge~ healthy.tol.dose}}{V_{healthy}}}{\frac{V_{target ~\ge~ presc.dose}}{V_{target}}}}{ascii} 
#' }

#' @return \mjeqn{-~COSI}{ascii} : if "COSI" is requested in \code{conformity.indices}, 
#' it returns a dataframe of Critical Organ Scoring Index for each healthy organ, 
#' at each \code{presc.dose}, and for each target. COSI is defined by 
#' \emph{Menhel and al} \strong{\[13\]}
#' \mjdeqn{COSI  = 1- 
#' \frac{\frac{V_{healthy ~\ge~ healthy.tol.dose}}{V_{healthy}}}{\frac{V_{target ~\ge~ presc.dose}}{V_{target}}}}{ascii} 

#' @return \mjeqn{-~homogeneity}{ascii} : dataframe containing
#' \itemize{
#' \item \code{HI.RTOG.max_ref} : (1) Homogeneity Index from RTOG defined by 
#' \emph{E.Shaw} \strong{\[1\]} 
#' \mjdeqn{HI.RTOG.max_-ref = \frac{D_{~max}}{presc.dose}}{ascii}
#' where \mjeqn{D_{max}}{ascii} is the maximum dose in the target volume.

#' \item \code{HI.RTOG.5_95} : (1) Homogeneity Index from RTOG \strong{\[1\]}
#' \mjdeqn{HI.RTOG.5_-95 = \frac{D.5pc}{D.95pc}}{ascii}
#' where \mjeqn{D.5pc}{ascii} and \mjeqn{D.95pc}{ascii} are respectively the doses 
#' at 5% and 95% of the target volume in cumulative dose-volume histogram.
  
#' \item \code{HI.ICRU.max_min} : (1) Homogeneity Index from \emph{ICRU report 62} 
#' \strong{\[14\]}
#' \mjdeqn{HI.ICRU.max_-min = \frac{D_{~max}}{D_{~min}}}{ascii}
#' where \mjeqn{D_{max}}{ascii} and \mjeqn{D_{min}}{ascii} are respectively the 
#' maximum and the minimum dose in the target volume.
 
# \item \code{HI.ICRU.max_ref} : (1) Homogeneity Index from ICRU 
# \mjdeqn{HI.ICRU.max_-ref = \frac{D_{~max}}{presc.dose}}{ascii}

#' \item \code{HI.ICRU.2.98_ref} : (1) Homogeneity Index from \emph{ICRU report 83} 
#' \strong{\[15\]}
#' \mjdeqn{HI.ICRU.2.98_-ref = 100 \frac{D.2pc - D.98pc}{presc.dose}}{ascii}
#' where \mjeqn{D.2pc}{ascii} and \mjeqn{D.98pc}{ascii} are respectively the doses 
#' at 2% and 98% of the target volume in cumulative dose-volume histogram.

#' \item \code{HI.ICRU.2.98_50 } : (1) Homogeneity Index from \emph{ICRU report 83} 
#' \strong{\[15\]}
#' \mjdeqn{HI.ICRU.2.98_-50 = 100 \frac{D.2pc - D.98pc}{D.50pc}}{ascii}
#' where \mjeqn{D.2pc}{ascii}, \mjeqn{D.98pc}{ascii} and \mjeqn{D.50pc}{ascii} are 
#' respectively the doses 
#' at 2%, 98% and 50% of the target volume in cumulative dose-volume histogram.

#' \item \code{HI.ICRU.5.95_ref} : (1) Homogeneity Index from \emph{ICRU report 83} 
#' \strong{\[15\]}
#' \mjdeqn{HI.ICRU.5.95_-ref = 100 \frac{D.5pc - D.95pc}{presc.dose}}{ascii}
#' where \mjeqn{D.5pc}{ascii} and \mjeqn{D.95pc}{ascii} are respectively the doses 
#' at 5% and 95% of the target volume in cumulative dose-volume histogram.

#' \item \code{HI.mayo2010} : (1) Homogeneity Index defined by \emph{Mayo and al} 
#' \strong{\[16\]}
#' \mjdeqn{HI.mayo2010 =\sqrt{\frac{D_{~max}}{presc.dose}~\cdot~(1 + 
#' \frac{\sigma_D}{presc.dose})}}{ascii}
#' where \mjeqn{D_{max}}{ascii} is the maximum dose in the target volume, and
#' \mjeqn{\sigma_D}{ascii} the standard deviation of the dose in the target volume.

#' \item \code{HI.heufelder} : (1) Homogeneity Index defined by \emph{Heufelder and al} 
#' \strong{\[17\]}
#' \mjdeqn{HI.heufelder = e^{-0.01~\cdot~ (1-\frac{\mu_D}{presc.dose})^2}~\cdot~
#' e^{-0.01~\cdot~ (\frac{\sigma_D}{presc.dose})^2}}{ascii}
#' where \mjeqn{\mu_D}{ascii} and \mjeqn{\sigma_D}{ascii} are
#' respectively the mean and the standard deviation of the dose in the target volume.
#' }

#' @return \mjeqn{-~gradient}{ascii} : dataframe containing 
#' \itemize{
#' \item \code{GI.ratio.50}: Gradient Index based on volumes ratio defined by 
#' \emph{Paddick and Lippitz} \strong{\[18\]}
#' \mjdeqn{GI.ratio.50 = \frac {V_{0.5~\cdot~presc.dose}}{V_{presc.dose}}}{ascii}

#' \item \code{mGI}: Modified Gradient Index defined by \emph{SABR UK Consortium 2019} 
#' \strong{\[2\]}
#' \mjdeqn{mGI = \frac{V_{0.5~\cdot~presc.dose}}{V_{target ~\ge~ presc.dose}} = 
#' \frac{V_{0.5~\cdot~presc.dose}}{V_{target} ~\cap~ V_{presc.dose}}}{ascii}


#' }


# @return \mjeqn{-~radiobiology}{ascii}: For all target and healthy structures,
# dataframe containing:
# \itemize{
# \item \code{gEUD}: Generalized  equivalent uniform dose, defined by  \emph{Niemierko} \strong{\[18\]}
# \mjdeqn{gEUD = \left(\frac{1}{N}~\cdot~ \sum_i D_i^a\right)^{\frac{1}{a}}}{ascii}
# where \mjeqn{D_i}{ascii} is the dose in each voxel of the target or healthy structure, 
# \mjeqn{N}{ascii} is the number of voxels in the target or healthy structure, 
# and \mjeqn{a}{ascii} is the parameter \code{target.a} or \code{healthy.a}.
# }
#' 

#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{SHAW19931231}{espadon}
#' @references \strong{\[2\]} \insertRef{SABR}{espadon}
#' @references \strong{\[3\]} \insertRef{LOMAX20031409}{espadon}
#' @references \strong{\[4\]} \insertRef{RIET1997731}{espadon}
#' @references \strong{\[5\]} \insertRef{Paddick2000ASS}{espadon}
#' @references \strong{\[6\]} \insertRef{Nakamura2002}{espadon}
#' @references \strong{\[7\]} \insertRef{DICE1945}{espadon}
#' @references \strong{\[8\]} \insertRef{Park2014}{espadon}
#' @references \strong{\[9\]} \insertRef{Wu2003}{espadon}
#' @references \strong{\[10\]} \insertRef{Ansari2018}{espadon}
#' @references \strong{\[11\]} \insertRef{LEFKOPOULOS2000}{espadon}
#' @references \strong{\[12\]} \insertRef{Baltas1998ACI}{espadon}
#' @references \strong{\[13\]} \insertRef{Menhel2006}{espadon}
#' @references \strong{\[14\]} \insertRef{ICRU621}{espadon}
#' @references \strong{\[15\]} \insertRef{ICRU83}{espadon}
#' @references \strong{\[16\]} \insertRef{MAYO20101457}{espadon}
#' @references \strong{\[17\]} \insertRef{HEUFELDER2003231}{espadon}
#' @references \strong{\[18\]} \insertRef{Paddick2006ASD}{espadon}

#'  
#' @references All this references are compiled by 
#' \itemize{
#' \item \insertRef{KAPLAN2021164}{espadon} and 
#' \item \insertRef{PATEL2020}{espadon}.
#' }
# @references \strong{\[18\]} \insertRef{niemierko1999generalized}{espadon}

#' @details If \code{target.roi.name}, \code{target.roi.sname}, and \code{target.roi.idx} 
#' are all set to \code{NULL},no target RoI are selected.
#' 
#' @details If \code{healthy.roi.name}, \code{healthy.roi.sname}, and 
#' \code{healthy.roi.idx} are all set to \code{NULL}, no healthy RoI are selected.
#' 
#' @seealso \link[espadon]{rt.indices.from.bin}.
#' @examples

#' # loading of toy-patient objects (decrease dxyz and increase beam.nb 
#' #  for better result)
#' step <- 5

#' patient <- toy.load.patient (modality = c("rtdose", "rtstruct"), roi.name = "eye",
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' indices <- rt.indices.from.roi (patient$rtdose[[1]],  patient$rtstruct[[1]],
#'                                 target.roi.sname = "ptv",
#'                                 healthy.roi.sname = "eye", presc.dose = 50,
#'                                 conformity.indices = c("PITV", "PDS", "CI.lomax2003", 
#'                                                        "CN", "NCI", "DSC","COIN"),
#'                                 verbose = FALSE)
#' indices



#' @import progress
#' @import mathjaxr
#' @importFrom methods is
#' @export
rt.indices.from.roi <- function (vol, struct = NULL, T.MAT = NULL, 
                                 target.roi.name = NULL, target.roi.sname = NULL, target.roi.idx = NULL,  
                                 healthy.roi.name = NULL, healthy.roi.sname = NULL, healthy.roi.idx = NULL,
                                 presc.dose = NA,
                                 healthy.tol.dose = NA,
                                 healthy.weight = 1,
                                 dosimetry = c("D.min", "D.max", "D.mean", "STD"),
                                 volume.indices = c("V.tot", "area", "V.prescdose"),
                                 conformity.indices =  c("PITV","PDS", "CI.lomax2003", "CN",
                                                         "NCI", "DSC", "CI.distance", 
                                                         "CI.abs_distance", "CDI", "CS3", 
                                                         "ULF","OHTF", "gCI", 
                                                         "COIN", "G_COSI", "COSI"),
                                 homogeneity.indices = c("HI.RTOG.max_ref", "HI.RTOG.5_95", 
                                                         "HI.ICRU.max_min", 
                                                         "HI.ICRU.2.98_ref", "HI.ICRU.2.98_50", 
                                                         "HI.ICRU.5.95_ref", "HI.mayo2010", 
                                                         "HI.heufelder"),
                                 gradient.indices = c("GI.ratio.50", "mGI"),
                                 D.xpc = NULL, D.xcc = NULL, V.xGy = NULL, 
                                 
                                 verbose = TRUE){

  if (!is (vol, "volume")) stop ("vol should be a volume class object.")
  if (is.null(vol$vol3D.data)) stop ("empty vol$vol3D.data.")
  if (vol$modality !="rtdose") warning ("vol should be of rtdose modality.")
  
  
  if (!is.null(struct)){ 
    if (!is (struct, "struct")) stop ("struct should be a struct class object.")
    if (is.null (struct$roi.data)) stop ("empty struct$roi.data.")
    struct <- struct.in.new.ref  (struct, new.ref.pseudo = vol$ref.pseudo, T.MAT =T.MAT)
    if (is.null(struct)) warning ("struct and vol do not share the same ref.pseudo. Specify T.MAT")
  }

  L.target.idx <- NULL
  L.healthy.idx <- NULL
  bin.L <- NULL
  vol.L <- NULL
  vol.names <- NULL
  
  # if (is.null(target.roi.name) &  is.null(target.roi.sname) & is.null(target.roi.idx))
  #   target.roi.sname <- "ptv"
  if (is.null(target.roi.name) &  is.null(target.roi.sname) & is.null(target.roi.idx)){
    target.roi.idx <- NULL
  } else {
    target.roi.idx <- select.names(struct$roi.info$roi.pseudo, roi.name=target.roi.name, 
                                   roi.sname=target.roi.sname, roi.idx=target.roi.idx)
  }  
  
 
  
  if (is.null(healthy.roi.name) &  is.null(healthy.roi.sname) & is.null(healthy.roi.idx)){
    healthy.roi.idx <- NULL
  } else {
    healthy.roi.idx <- select.names(struct$roi.info$roi.pseudo, roi.name=healthy.roi.name, 
                                    roi.sname=healthy.roi.sname, roi.idx=healthy.roi.idx)
  }  
  roi.idx <- c(target.roi.idx, healthy.roi.idx)
  

  
  if (!is.null(roi.idx) & !is.null(struct)) {
    if (!is.null(target.roi.idx))  L.target.idx <-match(target.roi.idx,roi.idx)
    if (!is.null(healthy.roi.idx))  L.healthy.idx <- match(healthy.roi.idx,roi.idx)
    
    margin <- 5
    if (verbose) pb <- progress_bar$new(format = " processing [:bar] :percent",
                                        total = length(roi.idx), clear = FALSE, width= 60)
    vol.L <- list ()
    D.max <- (floor(vol$max.pixel/10) + 1)*10
    for (idx in 1:length(roi.idx)) {
      vr <- nesting.roi(vol, struct, roi.idx =roi.idx[idx], xyz.margin=rep(margin,3), 
                        T.MAT = T.MAT, alias = struct$roi.info$name[roi.idx[idx]])
      b <- bin.from.roi(vr, struct, roi.idx = roi.idx[idx], T.MAT = T.MAT, 
                        alias = struct$roi.info$name[roi.idx[idx]])
      vr <- vol.from.bin(vr,b, alias = struct$roi.info$name[roi.idx[idx]], description=b$description)
      # h <- histo.vol (vr, breaks = seq (0, D.max, by = 0.001),
      #                 alias = struct$roi.info$name[roi.idx[idx]])
      vol.L[[idx]] <- list(b=b, v=vr)#, h = h)
      if (verbose) pb$tick()
    }
    vol.names <- trimws(struct$roi.info$name[roi.idx])
    names(vol.L) <- vol.names
    
    bin.L <- lapply(vol.L, function(l) l$b)
    # histo.L <- lapply(vol.L, function(l) l$h)
    vol.L <- lapply(vol.L, function(l) l$v)
    
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
                    D.xpc, D.xcc, V.xGy, 
                    healthy.weight, healthy.tol.dose
                    )
  
  
  
  
}

################################################################################
#' @importFrom Rvcg vcgArea vcgClost
.indices.from.bin <- function (vol,presc.dose, bin.L,vol.L,vol.names, 
                               L.target.idx, L.healthy.idx,
                               dosimetry,
                               volume.indices,
                               conformity.indices,
                               homogeneity.indices,
                               gradient.indices,
                               D.xpc, D.xcc, V.xGy, 
                               healthy.weight, healthy.tol.dose){#, DVH){
  
  vol.over.val.df <- function (vol,val,coln =NULL) {
    df <- as.data.frame (matrix (sapply (val, function(r) sum (vol$vol3D.data>=r, na.rm=T) * 
                                           abs (prod(vol$dxyz))/1000), ncol=length (val)))
    if (is.null(coln)) {colnames(df)   <- val} 
    else {colnames(df) <-  coln }
    return (df)
  }

  
  ###################################################################################
  #volume.indices
  
  V.xGy.b <- data.frame(matrix(ncol= 0,nrow=length(vol.names)+1,dimnames= list(c("isodose",vol.names),NULL)))
  if (!is.null (V.xGy))  {
    VxGy.names <- paste("V_",V.xGy,"Gy", sep="")
    V.xGy.b <- do.call(rbind,lapply(vol.L,function(V)  vol.over.val.df (V, V.xGy,VxGy.names)))
    isodose.xGy <- vol.over.val.df (vol, V.xGy,VxGy.names)
    V.xGy.b <- rbind (isodose =isodose.xGy,V.xGy.b)
    # dvh.interpol <-function(d , dvh){
    #   idx <- floor(d/dvh$step) + 1
    #   dvh$vol[idx] + (dvh$vol[idx+1]-dvh$vol[idx])*(d-dvh$breaks[idx])/dvh$step
    # }
    # Vx.Gy_ <- do.call( rbind.data.frame,lapply(vol.names, function(i) sapply(V.xGy, 
    #                                      function(x) dvh.interpol(x,dvh.L[[i]]))))
    # rownames (Vx.Gy_) <- vol.names
  }
  
 
  ##########################

  
  V.refdose <- NA
  V.50pc.refdose <- NA
  V.95pc.refdose <- NA
  V.105pc.refdose <- NA
  if (!is.na(presc.dose[1]))  {
    V.refdose <- as.data.frame (matrix (sapply (presc.dose, function(r) 
      sum (vol$vol3D.data>=r, na.rm=T) * abs (prod(vol$dxyz))/1000),
      ncol = length (presc.dose), dimnames=list("V.refdose",presc.dose)))
    V.50pc.refdose <- as.data.frame (matrix (sapply (0.5*presc.dose, function(r) 
      sum (vol$vol3D.data>=r, na.rm=T) * abs (prod(vol$dxyz))/1000),
      ncol = length (presc.dose), dimnames=list("V.refdose",0.5*presc.dose)))
    if ("CS3" %in% conformity.indices) 
      V.95pc.refdose <- as.data.frame (matrix (sapply (0.95*presc.dose, function(r) 
        sum (vol$vol3D.data>=r, na.rm=T) * abs (prod(vol$dxyz))/1000),
        ncol = length (presc.dose), dimnames=list("V.refdose",0.95*presc.dose)))
    if ("CS3" %in% conformity.indices) 
      V.105pc.refdose <- as.data.frame (matrix (sapply (1.05*presc.dose, function(r) 
        sum (vol$vol3D.data>=r, na.rm=T) * abs (prod(vol$dxyz))/1000),
        ncol = length (presc.dose), dimnames=list("V.refdose",1.05*presc.dose)))
  }
  
  if (is.null(bin.L)) {
    if (is.na(presc.dose[1]) & ncol(V.xGy.b)==0) return(NULL)
    if (is.na(presc.dose[1])) return(list(volume.indices=V.xGy.b))
    volume.indices.b <- V.refdose
    colnames(volume.indices.b) <- paste0("V_",presc.dose,"Gy")
    if (ncol(volume.indices.b)>0 & ncol(V.xGy.b)>0){
      m <- match(colnames(volume.indices.b), colnames(V.xGy.b))
      m <- m[!is.na(m)]
      if (length(m)>0) V.xGy.b <- V.xGy.b[-m]
    }
    
    volume.indices.b <- cbind (volume.indices.b, V.xGy.b)
    rownames(volume.indices.b) <- "isodose"
    volume.indices.b<-.reduc.tab(volume.indices.b)
    return(list(volume = volume.indices.b))
  }
  ###############################################################################
  dosimetry.names <- c("D.min", "D.max", "D.2%", "D.5%", "D.95%",
                       "D.98%","D.median", "D.mean", "STD")
  dosimetry.probs <- 1-c (1, 0, 0.02, 0.05, 0.95, 0.98, 0.5)
  
  
  Dxpc.names <- NULL
  if (!is.null (D.xpc))  {
    Dxpc.names <-  paste("D.",D.xpc,"%", sep="")
    probs.Dxpc <-  1- D.xpc/100
    probs.Dxpc [probs.Dxpc>1] <- 1
    probs.Dxpc [probs.Dxpc<0] <- 0
    f <- is.na(match(Dxpc.names,dosimetry.names))
    if (any(f)){
      dosimetry.probs <- c(dosimetry.probs,probs.Dxpc[f])
      dosimetry.names <- c(dosimetry.names[1:7],Dxpc.names[f], "D.mean", "STD")
    }
  }
  
  Dosimetry <- as.data.frame(do.call (rbind, lapply(vol.L,function(V) {
    D <- matrix(c(quantile (V$vol3D.data, probs=dosimetry.probs, na.rm=TRUE),
                  mean (V$vol3D.data, na.rm=TRUE), sd (V$vol3D.data, na.rm=TRUE)), ncol=length(dosimetry.names), byrow = T)
    dimnames(D) <- list(V$object.alias,dosimetry.names)
    D
  })))
  
  
  
  Dxcc.names <- NULL
  if (!is.null (D.xcc))  {
    Dxcc.names <-   paste("D",D.xcc,"cc", sep="")
    Dosimetry.xcc <- as.data.frame(do.call (rbind, lapply(vol.L,function(V) {
      probs.Dxcc <-  1-D.xcc/  (sum (V$vol3D.data>=0, na.rm = TRUE) * abs(prod (V$dxyz))/1000)
      probs.Dxcc [probs.Dxcc>1] <- 1
      probs.Dxcc [probs.Dxcc<0] <- 0
      D <- matrix(quantile (V$vol3D.data, probs=probs.Dxcc, na.rm=TRUE),
                  ncol=length(Dxcc.names), byrow = T)
      dimnames(D) <- list(V$object.alias,Dxcc.names)
      D
    })))  
    Dosimetry <- cbind(Dosimetry[,1:(ncol(Dosimetry)-2)], 
                       Dosimetry.xcc,Dosimetry[,(ncol(Dosimetry)-1):ncol(Dosimetry)])
  }
  
  

  
 ####################################################################################### 
  V.target <- do.call (rbind, lapply(vol.L[L.target.idx],function(V) vol.over.val.df (V, 0)))
  V.healthy <- do.call (rbind, lapply(vol.L[L.healthy.idx],function(V) vol.over.val.df (V, 0)))
  V.target.over.refdose <- NULL
  V.target.over.refdose.95pc <- NULL
  V.target.over.refdose.105pc <- NULL
  
  V.target.under.refdose <- NULL
  V.healthy.over.refdose <- NULL
  
  

  if (!is.na(presc.dose[1]) & !is.null(L.target.idx)){
    V.target.over.refdose <- do.call (rbind, lapply(vol.L[L.target.idx],function(V) vol.over.val.df (V, presc.dose)))
    # V.target.under.refdose <- do.call (rbind, lapply(vol.L[L.target.idx],function(V) vol.under.val.df (V, presc.dose)))
    V.target.under.refdose <- as.data.frame(-sweep(as.matrix(V.target.over.refdose),1,as.matrix(V.target)))
  }
  if (!is.na(presc.dose[1]) & !is.null(L.healthy.idx)){
    V.healthy.over.refdose <- do.call (rbind, lapply(vol.L[L.healthy.idx],function(V) vol.over.val.df (V, presc.dose)))
  }
  
  
  ##########################
  f.conformity <- any (!is.na (match (conformity.indices, c("CI.distance","CI.abs_distance","CDI"))))
  S.refdose <- NA
  if (f.conformity & !is.na(presc.dose[1])) {
    mesh.refdose <- lapply (presc.dose, function(r) {
      b <- bin.from.vol(vol, min=r)
      mesh.from.bin(b, alias=r, smooth.iteration = 10, smooth.type ="angWeight")
      })
    names(mesh.refdose) <- presc.dose
    
    S.refdose <- matrix(sapply ( mesh.refdose, function (m) {
      if (is.null(m)) return(0)
      vcgArea(m$mesh, perface = FALSE)}), dimnames=list(presc.dose,"area"))/100
  }
 
  if (!is.null(L.target.idx) &  (f.conformity | "area" %in% volume.indices)){   
    mesh.target <-lapply(bin.L[L.target.idx], function (b) 
      mesh.from.bin(b, alias=b$object.alias, smooth.iteration = 10, smooth.type ="angWeight")) 
    names(mesh.target) <- names (bin.L[L.target.idx])
    S.target <- matrix(sapply ( mesh.target, function (m) {
      if (is.null(m)) return(0)
      vcgArea(m$mesh, perface = FALSE)}), dimnames=list(vol.names[L.target.idx],"area"))/100
  }
  
  S.healthy <- NULL
  if ("area" %in% volume.indices  & !is.null(L.healthy.idx)){   
    mesh.healthy <-lapply(bin.L[L.healthy.idx], function (b) 
      mesh.from.bin(b, alias=b$object.alias, smooth.iteration = 10, smooth.type ="angWeight")) 
    names(mesh.healthy) <- names (bin.L[L.healthy.idx])
    S.healthy <- matrix(sapply ( mesh.healthy, function (m) {
      if (is.null(m)) return(0)
      vcgArea(m$mesh, perface = FALSE)}), dimnames=list(vol.names[L.healthy.idx],"area"))/100
  }


  if ((("CI.distance" %in% conformity.indices) |  ("CI.abs_distance" %in% conformity.indices)) 
      & !is.na(presc.dose[1]) & !is.null(L.target.idx)) {
    dist.fct<- function(vect) {sqrt(sum(vect^2))}
    G.target <- lapply(bin.L[L.target.idx], function (b)apply(get.xyz.from.index(which(b$vol3D.data), b),2, mean))
    
    unit.vector <-  matrix(.fansphereC(1),ncol=5, byrow=TRUE)[,1:3]
    
    ##### prec.dose
    Isodose.to.G.dist <- NULL
    if (length(presc.dose)>0)
      Isodose.to.G.dist <- lapply(1:length(presc.dose), function(pd.idx){
        dist <- matrix(NA, ncol=length(G.target), nrow= nrow(unit.vector), 
                       dimnames = list(NULL,names(G.target)))
        for (col in colnames(dist)){
          if (!is.null (mesh.refdose[[pd.idx]]))
          # if (get.value.from.xyz(G.target[[col]],vol,interpolate = TRUE) > presc.dose[pd.idx])
          if(vcgClost(matrix(G.target[[col]],ncol=3), mesh.refdose[[pd.idx]]$mesh)$quality>=0)
            dist[ ,col] <- .mesh.pt.distance(mesh.refdose[[pd.idx]]$mesh, 
                                           xyz.pt=G.target[[col]], unit.vector)
          # open3d()
          # par3d(windowRect=wr)
          # shade3d(mesh.refdose[[pd.idx]]$mesh)
          # points3d(matrix(G.target[[col]], ncol=3), col="red")
          # points3d(sweep(unit.vector*dist[,col],2,G.target[[col]],"+"), col="yellow")
        }
        dist
      })
    names(Isodose.to.G.dist) <- presc.dose
    
    Target.to.G.dist <- matrix(NA,ncol=length(G.target), nrow= nrow(unit.vector), 
                   dimnames = list(NULL,names(G.target)))
    for (col in colnames(Target.to.G.dist)){
      if (!is.null (mesh.target[[col]]))
      if (vcgClost(matrix(G.target[[col]],ncol=3), mesh.target[[col]]$mesh)$quality>=0)
        Target.to.G.dist[,col] <- .mesh.pt.distance(mesh.target[[col]]$mesh, 
                                       xyz.pt=G.target[[col]], unit.vector)
  
     
      # open3d()
      # par3d(windowRect=wr)
      # wire3d(mesh.target[[col]]$mesh,col="black")
      # points3d(matrix(G.target[[col]], ncol=3), col="red")
      # points3d(pt[na.idx,], col="yellow")
      # points3d(sweep(unit.vector*Target.to.G.dist[,col],2,G.target[[col]],"+"), col="yellow")
    }
  

  }


  
  ################################################
  
  if (!is.na(presc.dose[1]) & !is.null(L.target.idx)){
    tol.dose_ <- NULL
    healthy.weight_ <- NULL
    COSI.info <- NULL 
    if (any(!is.na (match(c("COSI","G_COSI"),conformity.indices))) & !is.null(L.healthy.idx)) {
      
      tol.dose_ <- rep(rev(healthy.tol.dose)[1],length(L.healthy.idx))
      tol.dose_ [1:length(healthy.tol.dose)] <- healthy.tol.dose
      healthy.weight_ <-  matrix(rep(rev(healthy.weight)[1],length(L.healthy.idx)), nrow=1)
      healthy.weight_ [1,1:length(healthy.weight)] <- healthy.weight
      
      COSI.info <- data.frame (matrix(NA, ncol =2, nrow =length(L.healthy.idx),
                                      dimnames = list(vol.names[L.healthy.idx], c("weight","tol.dose"))))
      
      COSI.info$weight <- healthy.weight_[1,]
      COSI.info$tol.dose <- tol.dose_
      
      V.healthy.over.toldose <- do.call (rbind, 
                                         lapply(1:length(L.healthy.idx),function(idx) {
                                           df <-vol.over.val.df (vol.L[[L.healthy.idx [idx]]],
                                                                 tol.dose_[idx])
                                           colnames(df) <- "toldose"
                                           df
                                         }))
      row.names(V.healthy.over.toldose) <- names(vol.L[L.healthy.idx])
      frac.healthy <- t(V.healthy.over.toldose/V.healthy)[1,]
      COSI <- as.data.frame(do.call (rbind, lapply (1:length(presc.dose), function(col.idx) {
        do.call (rbind, lapply(1:nrow(V.target), function (target.idx) 
          c (col.idx,target.idx,
             frac.healthy/(V.target.over.refdose[target.idx, col.idx]/V.target[target.idx,1]))))
      })))
      colnames(COSI) <- c("presc.dose", "target", vol.names[L.healthy.idx])
    } 
    
    m <- match(conformity.indices, c("PITV","PDS", "CI.lomax2003", "CN",
                                     "NCI", "DSC", "CI.distance", 
                                     "CI.abs_distance", "CDI", "CS3", 
                                     "ULF","OHTF", "gCI", 
                                     "COIN", "G_COSI"))
    
    cn <- conformity.indices[!is.na(m)]
    conformity.b <- data.frame (matrix(NA, ncol =length(cn), nrow =length(L.target.idx),
                                       dimnames = list(vol.names[L.target.idx], cn)))

    if (ncol(conformity.b)>0){
      
      conformity.b$target <- vol.names[L.target.idx]
      conformity.b$presc.dose <- NA
      conformity.b <- conformity.b[, c(ncol(conformity.b):(ncol(conformity.b)-1), 1:(ncol(conformity.b)-2))]
      
      conformity <- do.call (rbind, lapply(1:length(presc.dose), function(col.idx) {
        # CI.RTOG <- V.refdose[rep(1,nrow(V.target)),col.idx]/V.target
        # colnames(CI.RTOG) <- "CI.RTOG"
        df <- conformity.b
        df$presc.dose <- presc.dose[col.idx]
        if ("PITV" %in% colnames(df)) df$PITV <- as.numeric ((V.refdose[rep(1,nrow(V.target)),col.idx]/V.target)[,1])
        if ("PDS" %in% colnames(df)) df$PDS <- V.refdose[rep(1,nrow(V.target)),col.idx]/V.target.over.refdose[ ,col.idx]
        if ("CI.lomax2003" %in% colnames(df)) df$CI.lomax2003  <- as.numeric ((V.target.over.refdose[,col.idx]/
                                                                                 V.refdose[rep(1,nrow(V.target)),col.idx]))#[,1])
        CN  <- as.numeric (((V.target.over.refdose[,col.idx]^2) /
                              V.target / 
                              V.refdose[rep(1,nrow(V.target)),col.idx])[,1])
        if ("CN" %in% colnames(df)) df$CN <- CN
        if ("NCI" %in% colnames(df)) df$NCI  <- as.numeric (1/ CN)
        if ("DSC" %in% colnames(df)) df$DSC  <- as.numeric ((2 * V.target.over.refdose[,col.idx] / 
                                                               (V.target + V.refdose[rep(1,nrow(V.target)),col.idx]))[,1])
        
        
        if ("CI.distance" %in% colnames(df)) df$CI.distance  <- 100 * 
          as.numeric (apply((Isodose.to.G.dist[[col.idx]]-Target.to.G.dist)/Target.to.G.dist,2,mean, na.rm = TRUE))
        
        if ("CI.abs_distance" %in% colnames(df)) df$CI.abs_distance  <- 100 *
          as.numeric (apply(abs((Isodose.to.G.dist[[col.idx]]-Target.to.G.dist)/Target.to.G.dist),2,mean, na.rm = TRUE))
        
        
        if ("CDI" %in% colnames(df)) df$CDI <-as.numeric ((2* (V.refdose[rep(1,nrow(V.target)),col.idx] + 
                                                                 V.target - 2* V.target.over.refdose[,col.idx]) / 
                                                             (S.refdose[col.idx,rep(1,nrow(V.target))] + S.target))[,1]) 
        
        
        if ("CS3" %in% colnames(df)) df$CS3  <- as.numeric (((V.95pc.refdose[rep(1,nrow(V.target)),col.idx] + V.refdose[rep(1,nrow(V.target)),col.idx] +
                                                                V.105pc.refdose[rep(1,nrow(V.target)),col.idx])/ 
                                                               (3*V.target))[,1])
        
        
        if (any (!is.na(match (colnames(df),c("ULF","OHTF","gCI"))))){
          ULF<- as.numeric ((V.target.under.refdose[,col.idx] / V.target) [,1])
          OHTF <- rep(NA,length(L.target.idx))
          if (!is.null(V.healthy.over.refdose))
            OHTF <- as.numeric((sum(V.healthy.over.refdose[,col.idx]) / V.target)[,1])
          if ("ULF" %in% colnames(df)) df$ULF  <- ULF
          if ("OHTF" %in% colnames(df)) df$OHTF <- OHTF
          if ("gCI" %in% colnames(df)) df$gCI <- ULF + OHTF
        }
        
        if (!is.null(L.healthy.idx)){
          if ("COIN" %in% colnames(df)) 
            df$COIN <- as.numeric(((apply(1 - (V.healthy.over.refdose[,col.idx] /V.healthy),2,prod) *
                                     (V.target.over.refdose[,col.idx]^2) /V.target / 
                                     V.refdose[rep(1,nrow(V.target)),col.idx])[,1]))
          
          if ("G_COSI" %in% colnames(df) & !is.na(healthy.tol.dose[1])) {
            df$G_COSI <- 1- apply (as.matrix(COSI[COSI$presc.dose==col.idx, 3:ncol(COSI)], ncol=ncol(healthy.weight_))* 
                                    matrix(healthy.weight_[rep(1,nrow(V.target)),], ncol=ncol(healthy.weight_)),1,sum)
            }
        }
        df
        
      }))
      rownames(conformity) <- NULL
      conformity <- .reduc.tab(conformity)
    }else {
      conformity <- NULL
    }
    
    
    if ("COSI" %in% conformity.indices & !is.null(L.healthy.idx) & !is.na(healthy.tol.dose[1])) {
      COSI[ , 3:ncol(COSI)] <- 1-COSI[ , 3:ncol(COSI)]
      COSI$presc.dose <- presc.dose [COSI$presc.dose]
      COSI$target <- vol.names[COSI$target]
    }else COSI <- NULL
  } else {
    conformity <- NULL
    COSI <- NULL
    COSI.info <- NULL
  }

  ################################################ 
  homogeneity <- NULL
  if (!is.null(L.target.idx)){
    m <- match(homogeneity.indices,  c("HI.RTOG.max_ref", "HI.RTOG.5_95", 
                                       "HI.ICRU.max_min", #"HI.ICRU.max_ref", 
                                       "HI.ICRU.2.98_ref", "HI.ICRU.2.98_50", 
                                       "HI.ICRU.5.95_ref", "HI.mayo2010", 
                                       "HI.heufelder"))
    cn <- homogeneity.indices[!is.na(m)]
    homogeneity.b <- data.frame (matrix(NA, ncol =length(cn), nrow =length(L.target.idx),
                                        dimnames = list(vol.names[L.target.idx], cn)))
    
    if (ncol(homogeneity.b)>0){
      homogeneity.b$target <- vol.names[L.target.idx]
      homogeneity.b$presc.dose <- NA
      homogeneity.b <- homogeneity.b[, c(ncol(homogeneity.b):(ncol(homogeneity.b)-1), 1:(ncol(homogeneity.b)-2))]
      
      
      homogeneity <- as.data.frame(do.call (rbind, lapply (1:length(presc.dose), function(col.idx) {
        df <- homogeneity.b
        df$presc.dose <- presc.dose[col.idx]
        if ("HI.RTOG.max_ref" %in% colnames(df)) df$HI.RTOG.max_ref = Dosimetry[L.target.idx, "D.max"]/presc.dose[col.idx]
        if ("HI.RTOG.5_95" %in% colnames(df)) df$HI.RTOG.5_95 = Dosimetry[L.target.idx, "D.5%"] / Dosimetry[L.target.idx, "D.95%"] 
        if ("HI.ICRU.max_min" %in% colnames(df)) df$HI.ICRU.max_min = Dosimetry[L.target.idx, "D.max"] / Dosimetry[L.target.idx, "D.min"]
        # if ("HI.ICRU.max_ref" %in% colnames(df)) df$HI.ICRU.max_ref = Dosimetry[L.target.idx, "D.max"] / presc.dose[col.idx]
        if ("HI.ICRU.2.98_ref" %in% colnames(df)) df$HI.ICRU.2.98_ref = 100 * (Dosimetry[L.target.idx, "D.2%"] - Dosimetry[L.target.idx, "D.98%"])/presc.dose[col.idx]
        if ("HI.ICRU.2.98_50" %in% colnames(df)) df$HI.ICRU.2.98_50 = 100 * (Dosimetry[L.target.idx, "D.2%"] - Dosimetry[L.target.idx, "D.98%"])/ Dosimetry[L.target.idx, "D.median"]
        if ("HI.ICRU.5.95_ref" %in% colnames(df)) df$HI.ICRU.5.95_ref = 100 * (Dosimetry[L.target.idx, "D.5%"] - Dosimetry[L.target.idx, "D.95%"])/presc.dose[col.idx]
        if ("HI.mayo2010" %in% colnames(df)) df$HI.mayo2010 = sqrt ((Dosimetry[L.target.idx, "D.max"] / presc.dose[col.idx])*(1+ Dosimetry[L.target.idx, "STD"]/presc.dose[col.idx]))
        if ("HI.heufelder" %in% colnames(df)) df$ HI.heufelder = exp (-0.01 *  ((1 - (Dosimetry[L.target.idx, "D.mean"] / presc.dose[col.idx]))^2 + (Dosimetry[L.target.idx, "STD"]/presc.dose[col.idx])^2))
        # if ("D.STD" %in% colnames(df)) df$D.STD= Dosimetry[, "STD"]
        df
      })))
      rownames(homogeneity) <- NULL
    }
    homogeneity <- .reduc.tab(homogeneity)
  }
  ################################################ 

  ### gradient
  gradient <- NULL
  if (!is.null(L.target.idx)){
    m <- match(gradient.indices, c("GI.ratio.50", "mGI"))
    cn <- gradient.indices[!is.na(m)]
    gradient.b <- data.frame (matrix(NA, ncol =length(cn), nrow =length(L.target.idx),
                                     dimnames = list(vol.names[L.target.idx], cn)))
    
    
    if (!is.na(presc.dose[1]) & ncol(gradient.b)>0){
      gradient.b$target <- vol.names[L.target.idx]
      gradient.b$presc.dose <- NA
      gradient.b <- gradient.b[, c(ncol(gradient.b):(ncol(gradient.b)-1), 1:(ncol(gradient.b)-2))]
      gradient <- as.data.frame(do.call (rbind, lapply (1:length(presc.dose), function(col.idx) {
        df <- gradient.b
        df$presc.dose <- presc.dose[col.idx]
        if ("GI.ratio.50" %in% colnames(df)) 
          df$GI.ratio.50 = (V.50pc.refdose[rep(1,nrow(V.target)),col.idx]/ V.refdose[rep(1,nrow(V.target)),col.idx])
        
        if ("mGI" %in% colnames(df) & all(!is.na (V.target.over.refdose)))
          df$mGI = (V.50pc.refdose[rep(1,nrow(V.target)),col.idx]/ V.target.over.refdose[,col.idx])
        df
      })))
      rownames(gradient) <- NULL
      gradient <- .reduc.tab(gradient)
    } 
    }
    
 
  
  ###################################################################################
  
  # 
  # dvh.l <- NULL
  # if (DVH)
  #   dvh.l <- lapply(histo.L,function(l) histo.DVH(l, alias= l$alias))
  # 
  # m <- match(radiobiology.indices, c("gEUD", "TCP","NTCP"))
  # cn <- radiobiology.indices[!is.na(m)]
  # biology.b <- data.frame (matrix(NA, ncol =length(cn), nrow = length(bin.L),
  #                                    dimnames = list(vol.names, cn)))
  # 
  # if (ncol(biology.b)>0){
  # 
  #   if ("gEUD" %in% radiobiology.indices ){
  #     target.a_ <- rep(rev(target.a)[1],length(L.target.idx))
  #     target.a_ [1:length(target.a)] <- target.a 
  #     healthy.a_ <- NULL
  #     if (!is.null(L.healthy.idx)){
  #       healthy.a_ <- rep(rev(healthy.a)[1],length(L.healthy.idx))
  #       healthy.a_ [1:length(healthy.a)] <- healthy.a  
  #     }
  #     a <- c(target.a_, healthy.a_)
  #     biology.b$a <- a
  #     geud.col <-which(colnames(biology.b)=="gEUD")
  #     if (geud.col>1) {
  #       biology.b <- biology.b[,c(1:(geud.col-1),ncol(biology.b),geud.col:(ncol(biology.b)-1))]
  #     } else {
  #       biology.b <- biology.b[,c(ncol(biology.b),geud.col:(ncol(biology.b)-1))]
  #     }
  #     # Dosimetry$gEUD <- sapply(1:length(a), function(idx) {
  #     #   N <- sum (histo.L[[idx]]$counts)
  #     #   (sum (histo.L[[idx]]$counts * (histo.L[[idx]]$mids)^a[idx])/N)^(1/a[idx])
  #     #   })
  #     # # geud <- sapply(1:length(a), function(idx) sum (histo.L[[idx]]$counts / 
  #     #                                                  sum (histo.L[[idx]]$counts) * 
  #     #                                                  (histo.L[[idx]]$mids)^a[idx])^(1/a[idx]))
  #     
  #     biology.b$gEUD <-  sapply(1:length(a), function(idx) 
  #       (sum(vol.L[[idx]]$vol3D.data^a[idx],na.rm = T)/sum(bin.L[[idx]]$vol3D.data,na.rm = T))^(1/a[idx]))
  #     
  #     biology.b <- .reduc.tab(biology.b)
  #   }
  # 
  # } else {biology.b=NULL}
  # 
  # 

  
  #volume.indices suite
  m <- match (volume.indices, c("V.tot","area","V.prescdose"))
  volume.indices.b <- data.frame(matrix(ncol= 0,nrow=length(vol.names)+1,
                                   dimnames= list(c("isodose",vol.names),NULL)))
  if (any (!is.na(m))) {
    volume.indices.b <- lapply( volume.indices[!is.na(m)], function(i) NULL)
    names(volume.indices.b) <-  volume.indices[!is.na(m)]
    cn <-volume.indices.b
    if ("V.tot" %in% volume.indices) {
      volume.indices.b[["V.tot"]] <- rbind(isodose=NA, V.target,V.healthy)
      cn[["V.tot"]] <- "V_tot"
    }
    if ("V.prescdose" %in% volume.indices & !is.na(presc.dose[1]))  {
      volume.indices.b[["V.prescdose"]] <-rbind(isodose=V.refdose, V.target.over.refdose,V.healthy.over.refdose)
      cn[["V.prescdose"]] <- paste0("V_",colnames( volume.indices.b[["V.prescdose"]]),"Gy")
    }
    
    if ("area" %in% volume.indices) {
      volume.indices.b[["area"]] <-rbind(isodose=NA,S.target, S.healthy)
      cn[["area"]] <-"area"
    }
    
    volume.indices.b <- do.call(cbind, volume.indices.b[!sapply(volume.indices.b,is.null)])
    cn <- as.character(do.call(c, cn[!sapply(cn,is.null)]))
    colnames(volume.indices.b) <- cn
  }
  
  if (ncol(volume.indices.b)>0 & ncol(V.xGy.b)>0){
    m <- match(colnames(volume.indices.b), colnames(V.xGy.b))
    m <- m[!is.na(m)]
    if (length(m)>0) V.xGy.b <- V.xGy.b[,-m]
  }
  
  volume.indices.b <- cbind (volume.indices.b, V.xGy.b)
  if (ncol(volume.indices.b)==0) {
    volume.indices.b <- NULL
  } else {volume.indices.b <- .reduc.tab(volume.indices.b)}
  

  ###################################################################################
  #dosimetry suite
  dosimetry.name <- c(dosimetry,Dxpc.names,Dxcc.names)
  m <- match(c(dosimetry.name), colnames(Dosimetry))
  m <- m[!is.na(m)]
  dosimetry.b <- NULL
  if(length(m)>0) dosimetry.b  <- Dosimetry[,m]

  L <- list( dosimetry = dosimetry.b, volume = volume.indices.b, conformity = conformity,
             COSI = COSI, COSI.info = COSI.info, 
             homogeneity = homogeneity, gradient = gradient)#, DVH=dvh.l)
  return (L[!sapply(L,is.null)])
}




