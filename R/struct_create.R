#' Struct creating from contours list
#' @description The \code{struct.create} function creates a struct object from a 
#' list of polygons, representing the contours of a shape.
#' @param contours.list List of 3-column matrices or dataframes providing contour 
#' coordinates in z steps.
#' @param roi.name Character string, representing the name of created RoI.
#' @param roi.nb Positive integer, representing the number of created RoI.
#' @param roi.color Color of the created RoI, in hex code format ("#RRGGBB").
#' @param roi.type Type of RoI, from among "", "EXTERNAL", "PTV", "CTV", "GTV", 
#' "TREATED_VOLUME", "IRRAD_VOLUME", "OAR", "BOLUS", "AVOIDANCE", "ORGAN", "MARKER", 
#' "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", "CAVITY", "BRACHY_CHANNEL", 
#' "BRACHY_ACCESSORY", "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", "SUPPORT", "FIXATION", 
#' "DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT".
#' @param ref.pseudo Character string, frame of reference pseudonym of the 
#' created object.By defaukt equal to "ref1"
#' @param frame.of.reference Character string, frame of reference of the 
#' created object.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' @return Returns a "struct" class object (see \link[espadon]{espadon.class}
#' for class definition), including the unique \code{roi.name} as region of interest.
#' @seealso \link[espadon]{struct.from.mesh}.
#' @examples
#' contours.z <- -50:50
#' theta <- seq(0,2*pi, length.out = 100)
#' contours <- lapply(contours.z,function(z){
#'   if (z<(-25)) return(data.frame(x = (50 + z) * cos(theta),
#'                                  y = (50 + z) * sin(theta),
#'                                  z = z))
#'   if (z>25) return(data.frame(x = (50 - z) * cos(theta),
#'                               y = (50 - z) * sin(theta),
#'                               z = z))
#'   return(data.frame(x = 25 * cos(theta),
#'                     y = 25 * sin(theta),
#'                     z = z))
#' })
#' 
#' contours <- contours[!sapply(contours, is.null)]
#' S <- struct.create(contours, roi.name="myshape",
#'                    roi.nb = 1,
#'                    roi.color = "#ff0000",
#'                    roi.type = "",
#'                    ref.pseudo = "ref1", 
#'                   alias="", description = NULL)
#' display.3D.contour(S)

#' @export

struct.create <- function (contours.list, roi.name,
                           roi.nb = 1,
                           roi.color = "#ff0000",
                           roi.type = "",
                           ref.pseudo = "ref1", 
                           frame.of.reference = "",
                           alias="", description = NULL) {
  
  
  defined.type <- c ("","EXTERNAL", "PTV", "CTV", "GTV",
                     "TREATED_VOLUME", "IRRAD_VOLUME", "OAR",
                     "BOLUS", "AVOIDANCE", "ORGAN", "MARKER",
                     "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT",
                     "CAVITY", "BRACHY_CHANNEL", "BRACHY_ACCESSORY",
                     "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD",
                     "SUPPORT", "FIXATION", "DOSE_REGION",
                     "CONTROL", "DOSE_MEASUREMENT")
  roi.type <- roi.type[!is.na(match(roi.type,defined.type))][1]
  if (is.na(roi.type)) roi.type <- ""
  if (roi.name=="") roi.name <- "my_ROI"
  if (is.null(description)) description <- "ROI from my contours"
  roi.all.z <- sapply(contours.list, function(M) M[1,3])
  z <- sort(unique(roi.all.z))
  thickness <- 0
  if (length(z)>1) thickness <- min(diff(z))
  
  obj <- obj.create(class="struct")
  obj$description <- description
  obj$frame.of.reference <-frame.of.reference
  obj$ref.pseudo <- ref.pseudo
  obj$nb.of.roi <- 0
  obj$thickness <- thickness
  obj$roi.info[1,] <- c(roi.nb[1], roi.name,"","AUTOMATIC",roi.color,
                        castlow.str(roi.name))
  obj$object.alias <- alias[1]
  obj$description <- description
  obj$roi.obs[1,] <- c(1,roi.nb[1],"","","","","",roi.type,"")
  obj$roi.data <- list()
  obj$roi.data[[1]] <- lapply(contours.list,function(M) {
    M <- round(M,6)
    
    slope.x <- diff(M[,1])
    slope.y <- diff(M[,2])
    d <-  sqrt(slope.x^2 + slope.y^2)
    pt <- M[c(TRUE, !(diff(round(slope.x/d,6))==0 & diff(round(slope.y/d,6))==0), TRUE),]
    colnames(pt) <- c("x","y","z")
    row.names(pt) <- NULL
    list(type=ifelse(nrow(pt)==1,"POINT ",
                     ifelse(all(pt[1,1:2]==pt[nrow(pt),1:2]),"CLOSED_PLANAR ","OPEN_PLANAR ")),
         pt =as.data.frame(pt), level =0)})
  
  # on vérifie si les contours sont des contours inscrits ou non.
  
  ord <- order(roi.all.z)
  obj$roi.data[[1]] <- obj$roi.data[[1]][ord]
  roi.all.z <- roi.all.z[ord]
  if (length(roi.all.z)>0) {     
    kz <- rep(0,length(roi.all.z))
    if (obj$thickness!=0) kz <- round((roi.all.z -roi.all.z[1])/obj$thickness)
    
    for (k.value in unique(kz)){
      same.k.roi <- which (kz ==k.value)
      if (length(same.k.roi)>1) {
        for (j in same.k.roi){
          ptj<- obj$roi.data[[1]][[j]]$pt
          roi.index.k <-same.k.roi[same.k.roi!=j]
          # if (length(roi.index.z)!=0) {
          r <- unique (sapply (roi.index.k, function (k) {
            if (castlow.str (obj$roi.data[[1]][[k]]$type) != "closedplanar") return(NA)
            ptk <- obj$roi.data[[1]][[k]]$pt
            keep <- .pt.in.polygon (ptj[ ,1], ptj[ ,2],
                                              ptk[ ,1], ptk[ ,2]) > 0.5
            return (ifelse (all(keep), k,NA))}))
          r <- r[!is.na (r)]
          obj$roi.data[[1]][[j]]$level <- ifelse (length(r)!=0, length(r), 0)
        } #else obj$roi.data[[1]][[j]]$level <- 0
      } else obj$roi.data[[1]][[same.k.roi]]$level <- 0
    }
  }
  
  
  
  db <- .struct.moreinfo (obj$roi.data, obj$ref.from.contour,obj$thickness)
  obj$roi.info  <- cbind (obj$roi.info, db)
  
  class (obj) <- "struct"
  return(obj)
}