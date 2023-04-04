#' Gamma index 2D - 3D
#' @description The \code{rt.gamma.index} function computes the local or global  
#' Gamma index from a measurement and a reference. These latter are "volume" class
#' objects containing one (2D) or several planes (3D).
#' @param vol "volume" class object, which represents the measured volume.
#' @param vol.ref "volume" class object, which represents the reference volume.
#' @param over.sampling.factor Strictly positive integer, or a vector of 3 strictly 
#' positive integers, default to 1. Defined to oversample grids of \code{vol} and \code{vol.ref}.
#' Oversampling can be very time consuming.
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
#' @param project.to.isocenter Boolean. If \code{TRUE}, and if \code{vol} and 
#' \code{vol.ref} are of modality "rtimage", the size of the pixels is corrected 
#' to correspond to that found if the sensor was at the isocenter.
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
#' details the number of dose points, the number of evaluated dose points,the rate 
#' of evaluated dose points, the rate of Gamma indices below 1, above 1.2 and 1.5, 
#' the max and the mean Gamma index.
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
#' # by 3 pixels on the x axis
#' D.meas <- vol.copy (D.ref, alias = "measured_dose")
#' D.meas$vol3D.data[1:(D.meas$n.ijk[1] - 3) ,,] <- D.ref$vol3D.data[4:D.ref$n.ijk[1],,] 
#' D.max <- as.numeric(quantile(as.numeric(D.ref$vol3D.data), 
#'                              probs = 99.99/100, na.rm = TRUE))
#' gamma <- rt.gamma.index (D.meas, D.ref, delta.r = 6, vol.max = D.max)  
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
rt.gamma.index <- function (vol, vol.ref, 
                            over.sampling.factor = 1,
                            vol.max = vol.ref$max.pixel, 
                            dose.th = 0.02, delta.r = 3, 
                            analysis.th = 0.05, 
                            local = FALSE, local.th = 0.3, 
                            project.to.isocenter = TRUE,
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
  
  if (!((length(over.sampling.factor)==1 | length(over.sampling.factor)==3)) | 
      !all(abs(as.integer(over.sampling.factor)) == over.sampling.factor) | 
      any(over.sampling.factor < 0)){
    warning ("over.sampling.factor should be an integer >0 or vector of length 3 of integers >0")
    return (NULL)
  }
  
  if (length(over.sampling.factor)==1) over.sampling.factor=rep(over.sampling.factor,3)
  desc <-  paste0(ifelse(local, "local ", "global "),
                  dose.th*100, "% ", delta.r,"mm", 
                  ifelse(local, paste0(" - local th ", 100*local.th,"%"),""))
  if (is.null(description)) description <- desc
  
  if (project.to.isocenter & vol$modality == "rtimage" & vol.ref$modality == "rtimage") {
    vol <- .im.projection(vol)
    vol.ref <- .im.projection(vol.ref)
  }
  
  gammaindex <- vol.copy(vol.ref,alias = "gamma_index", modality ="gammaindex",
                         description = description)
  gammaindex$vol3D.data[] <- NA
  gammaindex$max.pixel  <-  gammaindex$min.pixel  <- NA
  
  
  a.th <- analysis.th * vol.max
  f.analyse <- !is.na(vol.ref$vol3D.data) & !is.na(vol$vol3D.data) & vol.ref$vol3D.data >= a.th
  
  inspect.idx <- which(f.analyse)
  le <- length(inspect.idx)
  
  gammaindex$setup <-  data.frame(label= c("Measure", "Reference", "Analysis threshold", "mode"),
                                  value = c(vol$object.alias, vol.ref$object.alias, 
                                            paste0(analysis.th * 100,"%"), desc))
  nb.pt <- prod(gammaindex$n.ijk)
  gammaindex$gamma.info <-  data.frame(label= c("nb of pts","evaluated pts","evaluated pts (%)",
                                                "<1 (%)","max", "mean",
                                                ">1.5 (%)",">1.2 (%)"),
                                       value = round(c(nb.pt, le, 100*le/nb.pt,
                                                       0,NA,NA,0,0),2))
  
  
  if (le==0) return(gammaindex)
  
  # search for max gamma
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
  
  t.mat <- ref.cutplane.add(gammaindex,ref.cutplane="rcp")
  # restrict volume
  bin.g <- vol.copy(vol.in.new.ref(gammaindex,"rcp",t.mat))
  bin.g$vol3D.data[] <- bin.g$min.pixel <-FALSE
  bin.g$vol3D.data[inspect.idx] <- bin.g$max.pixel <-TRUE
  bin.g$modality <- "binary"
  xyz.margin <- c(rep(max.gamma,3) * delta.r, 0)
  
  
  idx.c <- which(apply(abs(bin.g$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
  idx.r <-  which(apply(abs(bin.g$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    #2D
    u <- bin.g$xyz.from.ijk 
    u[idx.r,idx.c]<- 1
    ijk.from.xyz <- solve(u)
    ijk.from.xyz[idx.r,idx.c] <- 0
    ijk.margin <- xyz.margin %*% t(ijk.from.xyz)
  } else {
    ijk.margin <- xyz.margin %*% t(solve(bin.g$xyz.from.ijk))
  }
  
  ijk.margin[1:3][bin.g$n.ijk < 2] <- 0
  xyz.margin <- (ijk.margin %*% t(bin.g$xyz.from.ijk))[1:3]
  
  vol_ <- nesting.bin(vol.in.new.ref(vol,"rcp",t.mat),bin.g, 
                      xyz.margin = xyz.margin, vol.restrict = TRUE)
  vol.ref_<- nesting.bin(vol.in.new.ref(vol.ref,"rcp",t.mat),bin.g, 
                         xyz.margin = xyz.margin, vol.restrict = TRUE)
  
  f.analyse_ <- !is.na(vol.ref_$vol3D.data) & !is.na(vol_$vol3D.data) & vol.ref_$vol3D.data >= a.th
  inspect.idx_ <- which(f.analyse_)
  
  vol.ref_o <- vol.oversampling (vol.ref_, fact.ijk = over.sampling.factor)
  vol_o <- vol.oversampling (vol_, fact.ijk = over.sampling.factor)
  
  # display.kplane(vol.ref_,k=0, bg="green", interpolate =F)
  pt_ <- .get.ijkt.from.index(inspect.idx_,vol.ref_)
  # points(as.numeric(pt_[,1]),as.numeric(pt_[,2]),pch=".", col="red")
  
  # display.kplane(vol.ref_o,k=0, bg="green", interpolate =F)
  pt_o <- sweep(pt_[, 1:3],2,over.sampling.factor,"*")
  # points(as.numeric(pt_o[,1]),as.numeric(pt_o[,2]),pch=".", col="red")
  inspect.idx_o <-pt_o[,1] + vol.ref_o$n.ijk[1]*pt_o[,2] +
    prod(vol.ref_o$n.ijk[1:2])*pt_o[,3] + 1
  
  
  #ball definition
  #---------------
  
  idx.c <- which(apply(abs(vol.ref_o$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
  idx.r <-  which(apply(abs(vol.ref_o$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  
  if (length(idx.c)>0) {
    #2D
    u <- vol.ref_o$xyz.from.ijk 
    u[idx.r,idx.c]<- 1
    ijk.from.xyz <- solve(u)
    ijk.from.xyz[idx.r,idx.c] <- 0
    borne.ijk <- c(xyz.margin,0) %*% t(ijk.from.xyz)
  } else {
    borne.ijk <- c(xyz.margin,0) %*% t(solve(vol.ref_o$xyz.from.ijk))
  }
  borne.ijk[1:3][vol.ref_o$n.ijk < 2] <- 0
  borne.ijk <- ceiling(abs(borne.ijk)[1:3])
  
  
  ijkt = expand.grid((-borne.ijk[1]):borne.ijk[1], 
                     (-borne.ijk[2]):borne.ijk[2], 
                     (-borne.ijk[3]):borne.ijk[3],0)
  xyzt_rel <- (t(vol.ref_o$xyz.from.ijk %*% t (ijkt)))[,1:3]/delta.r
  s <- sqrt(xyzt_rel[,1]^2 + xyzt_rel[,2]^2  + xyzt_rel[,3]^2)
  
  f <-  s<=max.gamma
  ijk <- ijkt[f,1:3]
  s <- s[f]
  obr <- order(s)
  s <- s[obr]
  ball.ijk <- ijk[obr, ]
  ball.ijk <- ball.ijk[-1,]
  s <- s[-1]
  rownames(ball.ijk) <- NULL
  zone_s <- match(s,unique(s))
  around_idx <- max(which(zone_s<3))

  gammaindex$vol3D.data[inspect.idx] <- .gammaindex(vol3D = as.vector(vol_o$vol3D.data),
                                                    vol3D_ref = as.vector(vol.ref_o$vol3D.data), 
                                                    inspect_idx = inspect.idx_o - 1,
                                                    n_ijk = vol.ref_o$n.ijk, 
                                                    rel_dxyz= vol.ref_o$dxyz/delta.r,
                                                    ball_i= as.vector(ball.ijk[,1]),
                                                    ball_j= as.vector(ball.ijk[,2]),
                                                    ball_k= as.vector(ball.ijk[,3]),
                                                    around_idx = around_idx,
                                                    distance = s,
                                                    D_norm = vol.max,
                                                    local = local,
                                                    local_th_pc = local.th,
                                                    ref_pc = dose.th)
  gammaindex$max.pixel <- max(gammaindex$vol3D.data, na.rm = TRUE)
  gammaindex$min.pixel <- min(gammaindex$vol3D.data, na.rm = TRUE)
  
  
  gammaindex$gamma.info <-  data.frame(label= c("nb of pts","evaluated pts","evaluated pts (%)",
                                                "<1 (%)","max", "mean",
                                                ">1.5 (%)",">1.2 (%)"),
                                       value = round(c(nb.pt, le, 
                                                       100*le/nb.pt,
                                                       sum(gammaindex$vol3D.data[inspect.idx]<1)*100/le,
                                                       gammaindex$max.pixel,mean(gammaindex$vol3D.data[inspect.idx]),
                                                       sum(gammaindex$vol3D.data[inspect.idx]>1.5)*100/le,
                                                       sum(gammaindex$vol3D.data[inspect.idx]>1.2)*100/le),2))
  return(gammaindex)
}
