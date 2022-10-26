#' Gamma index 2D - 3D
#' @description The \code{rt.gamma.index} function computes the local or global  
#' Gamma index from a measurement and a reference. These latter are "volume" class
#' objects containing one (2D) or several planes (3D).
#' @param vol "volume" class object, which represents the measured volume.
#' @param vol.ref "volume" class object, which represents the reference volume.
#' @param vol.max Positive number, by default equal to the maximum value of the reference volume.
#' See Details.
#' @param dose.th Positive number, in percent, used to determine the dose difference criterion. See Details.
#' @param delta.r Positive number, in mm. Distance difference criterion.
#' @param analysis.th Positive number, in percent. Only the voxels whose value is 
#' greater than or equal \code{analyse.th*vol.max} are processed.
#' @param local Boolean. If \code{local = FALSE} (default), a global Gamma index 
#' is computed, and a local Gamma index otherwise.
#' @param local.th Positive number, in percent. Local threshold, only used if 
#' \code{local = TRUE}. See Details.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to Gamma index setup.
#' \loadmathjax
#' @details The Gamma index of a voxel \mjeqn{n}{ascii} was defined by \emph{Low and al} \strong{\[1\]}. 
#' It is computed from the formulae:
#' \mjdeqn{\gamma_n = min \left( \sqrt{\frac{(D_i-Dref_n)^2}{{\Delta D}^2} + \frac{{r_i}^2}{{\Delta r}^2}}\right)}{ascii}
#' whith \mjeqn{D_i}{ascii} the measured dose at voxel \mjeqn{i}{ascii},
#' \mjeqn{Dref_n}{ascii} the reference dose at voxel \mjeqn{n}{ascii},
#' \mjeqn{r_i}{ascii} the distance between voxels \mjeqn{i}{ascii} and \mjeqn{n}{ascii},
#' \mjeqn{\Delta r}{ascii} the distance difference criterion equal to \code{delta.r},
#' \mjeqn{\Delta D}{ascii} the distance difference criterion at voxel \mjeqn{n}{ascii} defined as follows:
#' \itemize{
#' \item If \code{local = FALSE} a global Gamma index is computed and 
#' \mjeqn{\Delta D~=~dose.th~\cdot~vol.max}{ascii}. 
#' \item If \code{local = TRUE}, then \mjeqn{\Delta D~=~dose.th~\cdot~Dref_n}{ascii} when  
#' \mjeqn{Dref_n~\ge~local.th~\cdot~vol.max}{ascii}, and
#' \mjeqn{\Delta D~=~dose.th~\cdot~local.th~\cdot~vol.max}{ascii} otherwise.
#' }
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions).  The \code{$vol3D.data} field represents the Gamma index. 
#' Two fields are added: 
#' the \code{$setup} field recalls the calculation setup, and the \code{$gamma.info} field
#' details the rate of Gamma indices below 1, above 1.2 and 1.5, the max and the 
#' mean Gamma index.
#' 
#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{low1998}{espadon}
#' 
#' @seealso \link[espadon]{rt.chi.index}
#' @examples
#' # Creation of a reference volume  and measured volume
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' patient <- toy.load.patient (modality = c ("rtdose", "rtstruct"), 
#'                              roi.name = "ptv", dxyz = c (3, 3, 3))
#' D.ref <- patient$rtdose[[1]]  
#' # We will assume that the measured dose is equal to the reference dose shifted
#' # by one pixel on the x axis
#' D.meas <- vol.copy (D.ref, alias = "measured_dose")
#' D.meas$vol3D.data[1:(D.meas$n.ijk[1] - 1) ,,] <- D.ref$vol3D.data[2:D.ref$n.ijk[1],,] 
#' 
#' gamma <- rt.gamma.index (D.meas, D.ref, delta.r = 6)  
#' gamma$gamma.info  
#' 
#' # Display gamma index at isocenter
#' G.iso <- patient$rtstruct[[1]]$roi.info$Gz[
#'   patient$rtstruct[[1]]$roi.info$name == "ptv"]
#' display.plane(gamma, view.coord = G.iso, 
#'               bottom.col = c ("#00FF00", "#007F00", "#FF8000", "#FF0000", 
#'                               "#AF0000"),
#'               bottom.breaks = c (0, 0.8, 1, 1.2, 1.5, gamma$max.pixel), 
#'               bg = "blue", interpolate = FALSE)


#' @export
rt.gamma.index <- function (vol, vol.ref, vol.max = vol.ref$max.pixel, 
                         dose.th = 0.02, delta.r = 3, 
                         analysis.th = 0.05, 
                         local = FALSE, local.th = 0.3, 
                         alias = "", description = NULL){
  
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if (!is (vol.ref, "volume")) {
    warning ("vol.ref should be a volume class object.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("Check input data : empty vol$vol3D.data.")
    return (NULL)
  }
  if(is.null(vol.ref$vol3D.data)){
    warning ("Check input data : empty vol.ref$vol3D.data.")
    return (NULL)
  }
  a.th <- analysis.th * vol.max
  f.analyse <- !is.na(vol.ref$vol3D.data) & !is.na(vol$vol3D.data) & vol.ref$vol3D.data >= a.th
  
  gamma <- vol$vol3D.data - vol.ref$vol3D.data
  gamma[!f.analyse] <- NA
  if (local){
    f <- vol.ref$vol3D.data[f.analyse]<local.th*vol.max
    gamma[f.analyse][f] <- abs(gamma[f.analyse][f]/(local.th*vol.max))
    gamma[f.analyse][!f] <- abs(gamma[f.analyse][!f]/(dose.th *  vol.ref$vol3D.data[f.analyse][!f]))
  } else {
    gamma[f.analyse] <- abs(gamma[f.analyse]/(dose.th*vol.max))
  }

  max.gamma <-  max(gamma,na.rm = TRUE)
  
  dxyz_rel = vol.ref$dxyz/delta.r
  borne.ijk <- ceiling(abs(max.gamma/dxyz_rel))
  borne.ijk[is.infinite(borne.ijk)] <- 0
  
  ijk = expand.grid((-borne.ijk[1]):borne.ijk[1], 
                    (-borne.ijk[2]):borne.ijk[2], 
                    (-borne.ijk[3]):borne.ijk[3])
  s <- sqrt((ijk[,1]*dxyz_rel[1])^2 + (ijk[,2]* dxyz_rel[2])^2  +( ijk[3] * dxyz_rel[3])^2)
  f <-  s<=max.gamma
  ijk <- ijk[f,]
  s <- s[f]
  obr <- order(s)
  s <- s[obr]
  ball.ijk <- ijk[obr, ]
  ball.ijk <- ball.ijk[-1,]
  s <- s[-1]
  rownames(ball.ijk) <- NULL
  around_idx <- max(which(s[1:7]==s[1]))
  zone_s <- match(s,unique(s))
  inspect.idx <- which(f.analyse)
  desc <-  paste0(ifelse(local, "local ", "global "),
                  dose.th*100, "% ", delta.r,"mm", 
                  ifelse(local, paste0(" - local th ", 100*local.th,"%"),""))
  if (is.null(description)) description <- desc
  
  gammaindex <- vol.copy(vol.ref,alias = "gamma_index", modality ="gammaindex",
                          description = description)
  gammaindex$vol3D.data[] <- NA
  
  gammaindex$vol3D.data[inspect.idx] <- .gammaindex(vol3D = as.vector(vol$vol3D.data),
                                                     vol3D_ref = as.vector(vol.ref$vol3D.data), 
                                                     inspect_idx = inspect.idx,
                                                     n_ijk = vol.ref$n.ijk, rel_dxyz= dxyz_rel,
                                                     ball_i= as.vector(ball.ijk[,1]),ball_j= as.vector(ball.ijk[,2]),
                                                     ball_k= as.vector(ball.ijk[,3]),
                                                     around_idx = around_idx,
                                                     distance = s,
                                                     D_norm = vol.max,
                                                     local = local,
                                                     local_th_pc = local.th,
                                                     ref_pc = dose.th)
  gammaindex$max.pixel <- max(gammaindex$vol3D.data, na.rm = TRUE)
  gammaindex$min.pixel <- min(gammaindex$vol3D.data, na.rm = TRUE)
 
 
  le <- length(inspect.idx)
  gammaindex$setup <-  data.frame(label= c("Measure", "Reference", "Analysis threshold", "mode"),
                                   value = c(vol$object.alias, vol.ref$object.alias, 
                                             paste0(analysis.th * 100,"%"), desc))
  gammaindex$gamma.info <-  data.frame(label= c("<1 (%)","max","mean",">1.5 (%)",">1.2 (%)"),
                                        value = round(c(sum(gammaindex$vol3D.data[inspect.idx]<1)*100/le,
                                                        gammaindex$max.pixel,mean(gammaindex$vol3D.data[inspect.idx]),
                                                        sum(gammaindex$vol3D.data[inspect.idx]>1.5)*100/le,
                                                        sum(gammaindex$vol3D.data[inspect.idx]>1.2)*100/le),2))
  return(gammaindex)
}