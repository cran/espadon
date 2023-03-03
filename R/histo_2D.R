#' 2D histograms of 2 volumes
#' @description The \code{histo.2D} function creates a "histo2D" class object, 
#' containing the two-dimensional array of histograms of two "volume" class 
#' objects  that have the same grid.
#' @param x.vol,y.vol "volume" class objects. The 2 volumes must have the
#' grid (i.e. share the same voxels location).
#' @param x.breaks,y.breaks Vectors giving the breakpoints of x and y axes. See 
#' Details.
#' @param include.outer Boolean. If \code{include.outer = TRUE}, the values out
#' the \code{x.breaks} and \code{y.breaks} of each volume are counted in the first
#' and the last cell of the histograms. They are not taken into account otherwise.
#' @param alias Character string, \code{$alias} of the created object
#' @param description Character string, describing the created object.
#' @details The arguments \code{x.breaks} and \code{y.breaks} represent the scales
#' of the x and y axes of 2D-histogram graph. If they are \code{NULL}, the 
#' \code{histo.2D} function defaults to 256 cells between the values 
#' \code{vol$min.pixel} and \code{vol$max.pixel}.
#' @return Returns a "histo2D" class object. This is a list including:
#' \itemize{
#' \item \code{$patient}: set to \code{x.vol$patient}.
#' \item \code{$patient.name}: set to \code{x.vol$patient.name}.
#' \item \code{$patient.bd}: set to \code{x.vol$patient.bd}.
#' \item \code{$patient.sex}: set to \code{x.vol$patient.sex}.
#' \item \code{$file.basename}: set to "".
#' \item \code{$file.dirname}: set to "".
#' \item \code{$object.name}: set to "".
#' \item \code{$object.alias}: alias of the histo2D object.
#' \item \code{$frame.of.reference}: set to \code{x.vol$frame.of.reference}.
#' \item \code{$ref.pseudo} : set to \code{x.vol$ref.pseudo}.
#' \item \code{$modality} : set to \code{"histo2D"}.
#' \item \code{$description}: description of the histo2D object.
#' \item \code{$creation.date}: set to \code{Sys.Date}.
#' \item \code{$nb.pixels}: number of elements in the \code{density.map}.
#' \item \code{$x.file.src}: set to x.vol$object.alias
#' \item \code{$y.file.src}: set to y.vol$object.alias
#' \item \code{x.breaks}: vector of x-axis breakpoints.
#' \item \code{y.breaks}: vector of y-axis breakpoints.
#' \item \code{x.mids}: vector of x-axis cell centers.
#' \item \code{y.mids}: vector of y-axis cell centers.
#' \item \code{density.map}: array of densities.
#' \item \code{total.counts}: number of counted voxels.
#' }
#' @seealso \link[espadon]{display.2D.histo}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "mr", "rtstruct"), 
#'                              roi.name =  "brain", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' MR <- patient$mr[[1]]
#' S <- patient$rtstruct[[1]]
#' T.MAT <- patient$T.MAT
#'
#' # restriction of the volume around the RoI
#' CT.on.roi <- nesting.roi (CT, S, roi.name = "brain", vol.restrict = TRUE,
#'                           xyz.margin = c (1, 1, 1), alias = CT$description)
#' MR.on.CT <- vol.regrid (vol = MR, back.vol = CT.on.roi, interpolate = TRUE,
#'                         T.MAT = T.MAT, alias = CT$description,
#'                         description = NULL)

#' # selection of voxels included in the RoI.
#' roi.bin <- bin.from.roi (vol = CT.on.roi, struct = S, roi.sname = "brain")
#' MR.select <- vol.from.bin (MR.on.CT, roi.bin, alias = MR$description)
#' CT.select <- vol.from.bin (CT.on.roi, roi.bin, alias = CT$description)

#' # 2D histogram
#' H2D <- histo.2D (MR.select, CT.select, x.breaks = seq (50, 400, 10),
#' 			  y.breaks = seq (50, 400, 10), alias = "H2D MR CT")
#' str (H2D)

#' @export
#' @importFrom methods is
histo.2D  <- function (x.vol, y.vol, x.breaks = NULL, y.breaks = NULL,
                       include.outer = TRUE, alias = "", description =""){
  
  
  .find.scale <- function(vol, breaks){
    if (all (is.na (vol$vol3D.data))) return(NULL)
    if (is.null (breaks)) {
      nb.mids <- 256
      min.vol <- min(vol$vol3D.data,na.rm = TRUE)
      max.vol <- max(vol$vol3D.data,na.rm = TRUE)
      
      step <- ifelse (max.vol-min.vol+1 < nb.mids, 1,((max.vol-min.vol)/(nb.mids-1)))
      nb.mids <- 1+(max.vol-min.vol)/step
      first <- min.vol
      last <- first+(nb.mids-1)*step
      breaks  <- seq (first- 0.5*step, last+0.5*step,step)
    } else {
      if (length(breaks)<2)return (NULL)
      stepu<-unique(diff(breaks))
      stepu <- stepu[!is.na(stepu)]
      if (length(stepu) == 0)return (NULL)
      step<-mean(stepu)
      if (any(abs(stepu-step)>1e-3))return (NULL)
      
      
      nb.mids <- length(breaks) - 1
      first <- mean(breaks[1:2])
      last <-  mean(rev(breaks[1:2]))
    }
    
    list(step=step,nb.mids=nb.mids,first=first, last=last,breaks=breaks)
  }
  
  if (!is (x.vol, "volume")){
    warning ("x.vol should be a volume class object.")
    return (NULL)
  }
  if (!is (y.vol, "volume")){
    warning ("y.vol should be a volume class object.")
    return (NULL)
  }
  
  if (is.null(x.vol$vol3D.data)| is.null(y.vol$vol3D.data)){
    warning ("no data in vol.")
    return (NULL)
  }
  
  if (!grid.equal (x.vol, y.vol)){
    warning ("both volumes must share same grid.")
    return (NULL)
  }
  param1 <- .find.scale (x.vol,x.breaks)
  if (is.null(param1)){
    warning ("undefined step for x.vol.")
    return (NULL)
  }
  param2 <- .find.scale (y.vol,y.breaks)
  if (is.null(param2)){
    warning ("undefined step for y.vol.")
    return (NULL)
  }
  
  scale1 <- c (0,(param1$nb.mids-1)) * param1$step + param1$first
  scale2 <- c(0,(param2$nb.mids-1)) * param2$step + param2$first
  
  H <- list()
  H$patient <- x.vol$patient
  H$patient.name <- x.vol$patient.name
  H$patient.bd <- x.vol$patient.bd
  H$patient.sex <- x.vol$patient.sex
  
  H$file.basename <- ""
  H$file.dirname <- ""
  H$object.name <- alias
  H$object.alias <- alias
  H$frame.of.reference <- x.vol$frame.of.reference
  H$ref.pseudo <- x.vol$ref.pseudo
  H$modality <- "histo2D"
  H$description <- description
  
  H$creation.date <- format(Sys.Date(), "%Y%m%d")
  
  H$x.file.src <- x.vol$object.alias
  H$y.file.src <- y.vol$object.alias
  H$x.breaks <- param1$breaks
  H$y.breaks <- param2$breaks
  H$x.mids <- param1$breaks[-1] - param1$step/2
  H$y.mids <- param2$breaks[-1] - param2$step/2
  H$nb.pixels <- 0
  
  vol1 <- floor((x.vol$vol3D.data-param1$first)/param1$step)
  vol2 <- floor((y.vol$vol3D.data-param2$first)/param2$step)
  
  if(include.outer) {
    vol1[vol1<0] <- 0
    vol1[vol1>param1$nb.mids -1] <- param1$nb.mids -1
    vol2[vol2<0] <- 0
    vol2[vol2>param2$nb.mids -1] <- param2$nb.mids -1
  } else {
    vol1[vol1<0] <- NA
    vol1[vol1>param1$nb.mids -1] <- NA
    vol2[vol2<0] <- NA
    vol2[vol2>param2$nb.mids -1] <- NA
  }
  H$density.map <- array(NA,dim=c(param1$nb.mids, param2$nb.mids))
  H$nb.pixels <- param1$nb.mids * param2$nb.mids
  pt <- data.frame(cbind (vol1[!is.na(vol1) & !is.na(vol2)],
                          vol2[!is.na(vol1) & !is.na(vol2)]))
  tab<-table(pt, useNA = "no")
  ctab<- lapply(dimnames(tab),function(v) as.numeric(v))
  H$density.map[as.matrix(expand.grid(ctab[[1]]+1,ctab[[2]]+1))] <- tab
  H$density.map[H$density.map==0] <- NA
  H$total.counts <- sum(H$density.map,na.rm = TRUE)
  H$density.map <- H$density.map/H$total.counts
  class (H) <- "histo2D"
  if (alias=="") return(H)
  return(.set.ref.obj(H,list(x.vol,y.vol)))
  
}