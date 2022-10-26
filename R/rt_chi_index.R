#' Chi index 2D - 3D
#' @description The \code{rt.chi.index} function computes the local or global  
#' Chi index from a measurement and a reference. These latter are "volume" class
#' objects containing one (2D) or several planes (3D).
#' @param vol "volume" class object, which represents the measured volume.
#' @param vol.ref "volume" class object, which represents the reference volume.
#' @param abs Boolean. If \code{TRUE} (default), the absolute value of Chi is computed.
#' @param vol.max Positive number, by default equal to the maximum value of the reference volume.
#' See Details.
#' @param dose.th Positive number, in percent, used to determine the dose difference criterion. See Details.
#' @param delta.r Positive number, in mm. Distance difference criterion.
#' @param analysis.th Positive number, in percent. Only the voxels whose value are 
#' greater than or equal \code{analyse.th * vol.max} are processed.
#' @param local Boolean. If \code{local = FALSE} (default), a global Chi index 
#' is computed, and a local Chi index otherwise.
#' @param local.th Positive number, in percent. Local threshold, only used if 
#' \code{local = TRUE}. See Details.
#' @param alias Character string, \code{$object.alias} of the created object.
#' @param description Character string, describing the created object. If 
#' \code{description = NULL} (default value), it will be set to Chi index setup.
#' \loadmathjax
#' @details The Chi index of a voxel \mjeqn{n}{ascii} was defined by \emph{Bakai and al} \strong{\[1\]}. 
#' It is computed from the formulae:
#' \mjdeqn{\chi_n = \frac{D_i-Dref_n}{\sqrt{{\Delta D}^2 + {\Delta r}^2~\cdot~\Vert \nabla Dref_n \Vert^2}}}{ascii}
#' If \code{abs = TRUE}, the used formulae is : 
#' \mjdeqn{\chi_n = \frac{\vert D_i-Dref_n\vert}{\sqrt{{\Delta D}^2 + {\Delta r}^2~\cdot~\Vert \nabla Dref_n \Vert^2}}}{ascii}
#' with \mjeqn{D_i}{ascii} the measured dose at voxel \mjeqn{i}{ascii},
#' \mjeqn{Dref_n}{ascii} the reference dose at voxel \mjeqn{n}{ascii},
#' \mjeqn{\nabla Dref_n}{ascii} the gradient of reference dose at voxel \mjeqn{n}{ascii},
#' \mjeqn{\Delta r}{ascii} the distance difference criterion equal to \code{delta.r}, and
#' \mjeqn{\Delta D}{ascii} the distance difference criterion at voxel \mjeqn{n}{ascii} defined as follows:
#' \itemize{
#' \item If \code{local = FALSE} a global Chi index is computed and 
#' \mjeqn{\Delta D~=~dose.th~\cdot~vol.max}{ascii}. 
#' \item If \code{local = TRUE}, then \mjeqn{\Delta D~=~dose.th~\cdot~Dref_n}{ascii} when  
#' \mjeqn{Dref_n~\ge~local.th~\cdot~vol.max}{ascii}, and
#' \mjeqn{\Delta D~=~dose.th~\cdot~local.th~\cdot~vol.max}{ascii} otherwise.
#' }
#' @return Returns a "volume" class object (see \link[espadon]{espadon.class} 
#' for class definitions).  The \code{$vol3D.data} field represents the Chi index. 
#' Two fields are added: 
#' the \code{$setup} field recalls the calculation setup, and the \code{$chi.info} field
#' details the rate of absolute values of the Chi index below 1, above 1.2 and 1.5, 
#' the max and the mean Chi index.

#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{bakai2003}{espadon}

#' @seealso \link[espadon]{rt.gamma.index}
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
#' abs_chi <- rt.chi.index (D.meas, D.ref, delta.r = 6)  
#' abs_chi$chi.info  
#' 
#' # Display chi index at isocenter
#' G.iso <- patient$rtstruct[[1]]$roi.info$Gz[
#'   patient$rtstruct[[1]]$roi.info$name == "ptv"]
#' display.plane(abs_chi, view.coord = G.iso, 
#'               bottom.col = c ("#00FF00", "#007F00", "#FF8000", "#FF0000", 
#'                               "#AF0000"),
#'               bottom.breaks = c (0, 0.8, 1, 1.2, 1.5, abs_chi$max.pixel),
#'               interpolate = FALSE, bg = "blue")

#' @export
rt.chi.index <- function (vol, vol.ref, abs = TRUE, vol.max = vol.ref$max.pixel, 
                          dose.th = 0.02, delta.r = 3, analysis.th = 0.05,
                          local = FALSE, local.th = 0.3, alias = "", description = NULL) {
  
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
  
  delta.D <- abs(vol$vol3D.data - vol.ref$vol3D.data)
  chi <- delta.D
  grad.D <-vol.gradient(vol.ref)
  
  if (local) {
    th.loc <- local.th*vol.max
    f <- chi>th.loc
    chi[f] <-  delta.D / sqrt(((dose.th * vol.ref$vol3D.data[f])^2) +
                                (delta.r^2) * (grad.D$vol3D.data[f])^2)
    chi[!f] <-  delta.D / sqrt(((dose.th * th.loc)^2) +
                                 (delta.r^2) * (grad.D$vol3D.data[!f])^2)
  } else {
    th <- vol.max * dose.th
    chi <- delta.D / sqrt((th^2) + (delta.r^2) * (grad.D$vol3D.data)^2)
  }
  chi[!f.analyse] <- NA
  desc <-  paste0(ifelse(local, "local ", "global "),
                  dose.th*100, "% ", delta.r,"mm", 
                  ifelse(local, paste0(" - local th ", 100*local.th,"%"),""))
  if (is.null(description)) description <- desc
  
  chiindex <- vol.copy(vol.ref,alias = "chi_index", modality ="chiindex",
                       description = description)
  if (abs) {
    chiindex$vol3D.data <- abs(chi)
  } else {
    chiindex$vol3D.data <- chi
  }
  chiindex$max.pixel <- max(chiindex$vol3D.data, na.rm = TRUE)
  chiindex$min.pixel <- min(chiindex$vol3D.data, na.rm = TRUE)
  
  le <- sum(f.analyse)
  chiindex$setup <-  data.frame(label= c("Measure", "Reference", "Analysis threshold", "mode"),
                                value = c(vol$object.alias, vol.ref$object.alias, 
                                          paste0(analysis.th * 100,"%"), desc))
  
  chiindex$chi.info <-  data.frame(label= c("<1 (%)","max","mean",">1.5 (%)",">1.2 (%)"),
                                   value = round(c(sum(abs(chi[f.analyse])<1)*100/le,
                                                   chiindex$max.pixel, mean(abs(chi[f.analyse])),
                                                   sum(abs(chi[f.analyse])>1.5)*100/le,
                                                   sum(abs(chi[f.analyse])>1.2)*100/le),2))
  
  return (chiindex)
  
}