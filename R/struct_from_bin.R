#' Creation of struct class object from a binary volume
#' @description The \code{struct.from.bin} function creates a struct object with 
#' a unique RoI, defined by the contours of binary volume.
#' @param vol "volume" class object, of binary modality.
#' @param roi.name Character string, representing the name of created RoI.
#' @param roi.nb Positive integer, representing the number of created RoI.
#' @param roi.color Color of the created RoI, in hex code format ("#RRGGBB").
#' @param roi.type Type of RoI, from among "", "EXTERNAL", "PTV", "CTV", "GTV", 
#' "TREATED_VOLUME", "IRRAD_VOLUME", "OAR", "BOLUS", "AVOIDANCE", "ORGAN", "MARKER", 
#' "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", "CAVITY", "BRACHY_CHANNEL", 
#' "BRACHY_ACCESSORY", "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", "SUPPORT", "FIXATION", 
#' "DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT"
#' @param external.only Boolean. If \code{TRUE}, only external contours are kept.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the created object.
#' @return Returns a "struct" class object (see \link[espadon]{espadon.class}
#' for class definition), including the unique \code{roi.name} as region of interest.
#' @examples
#' # Contours of a sphere of 10 mm radius
#' b.sphere <- vol.create (n.ijk = c (40, 40, 40), dxyz = c(1,1,1), 
#'                         mid.pt = c (0, 0, 0), modality = "binary", 
#'                         default.value = FALSE)
#' xyz <- expand.grid (-20:19, -20:19, -20:19)
#' R <- 10
#' Sphere.flag <- (xyz[, 1]^2 + xyz[, 2]^2 + xyz[, 3]^2) <= R^2
#' b.sphere$vol3D.data[Sphere.flag] <- TRUE
#' b.sphere$max.pixel <- TRUE
#' S.sphere <- struct.from.bin (b.sphere, roi.name = "sphere", external.only = TRUE)
#' str (S.sphere$roi.info)                 

## @importFrom sp point.in.polygon
#' @importFrom methods is
#' @export
struct.from.bin <- function (vol, roi.name = vol$description, roi.nb = 1,
                             roi.color = "#379DA2", 
                             roi.type = c ("","EXTERNAL", "PTV", "CTV", "GTV", 
                                           "TREATED_VOLUME", "IRRAD_VOLUME", "OAR", 
                                           "BOLUS", "AVOIDANCE", "ORGAN", "MARKER", 
                                           "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", 
                                           "CAVITY", "BRACHY_CHANNEL", "BRACHY_ACCESSORY", 
                                           "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", 
                                           "SUPPORT", "FIXATION", "DOSE_REGION", 
                                           "CONTROL", "DOSE_MEASUREMENT"),
                             external.only = FALSE, alias = "", 
                             description = paste ("RoI from", vol$object.alias)) {
  
  
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if ((vol$modality!="binary")) {
    stop ("vol must be of binary modality.")
  }
  
  t.mat <- ref.cutplane.add(vol,ref.cutplane = "intern", origin = c(0,0,0))
  
  
  
  bin_ <- vol.in.new.ref(vol, new.ref.pseudo="intern", t.mat)
  roi.name <- roi.name[1]
  roi.color <- roi.color[1]
  alias <- alias[1]
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
  
  struct <- list()
  
  struct$patient <-vol$patient
  struct$patient.name <- vol$patient.name
  struct$patient.bd <- vol$patient.bd
  struct$patient.sex <- vol$patient.sex
  struct$file.basename <- ""
  struct$file.dirname <- ""
  struct$object.name <- alias
  struct$object.alias <- alias
  struct$frame.of.reference <-  vol$frame.of.reference
  struct$ref.pseudo <- vol$ref.pseudo
  struct$modality <- "rtstruct"
  struct$description <- description
  d <- format (Sys.time(),"%Y%m%d")
  struct$acq.date <- vol$acq.date
  struct$study.date <- d
  struct$creation.date <- d
  struct$study.time <- ""
  
  struct$number <- 1
  struct$nb.of.roi <- 1
  if (vol$min.pixel==vol$max.pixel) struct$nb.of.roi <- 0
  
  struct$thickness <-round(vol$slice.thickness,6)
  struct$ref.from.contour <- get.rigid.M(t.mat,src.ref="intern",
                                         dest.ref = vol$ref.pseudo)
  # struct$ref.from.contour[,4] <- c(0,0,0,1)
  class(struct) <- "struct"
  if (struct$nb.of.roi == 0)  return(struct)
  
  struct$roi.info <- data.frame (number = roi.nb[1], name=roi.name, description="", 
                                 generation.algorithm = "automatic",
                                 color=roi.color, 
                                 roi.pseudo=tolower (gsub("[[:space:],_]", "", 
                                                          iconv (roi.name,  to="ASCII//TRANSLIT"))))
  
  struct$roi.obs <- data.frame (nb = roi.nb[1], roi.nb=roi.nb[1], label=roi.name, 
                                     code.value = "", code.scheme="",code.scheme.v="",
                                     code.meaning = "", type=roi.type, interpreter="")
                                     
  struct$roi.data <- list()
  # struct$roi.data[[1]] <- .display.roi.data.from.bin(bin_, bin_$xyz0[1,])
  struct$roi.data[[1]] <- .display.roi.data.from.bin(bin_)
  #contour fermé
  
  
  
  # on vérifie si les contours sont des contours inscrits ou non.
  roi.all.z<- sapply(struct$roi.data[[1]], function(li)  li$pt[1,3])
  if (length(roi.all.z)>0) {     
    kz <- rep(0,length(roi.all.z))
    if (struct$thickness>0) kz <- round((roi.all.z -roi.all.z[1])/struct$thickness)
    
    for (k.value in unique(kz)){
      same.k.roi <- which (kz ==k.value)
      if (length(same.k.roi)>1) {
        for (j in same.k.roi){
          ptj<- struct$roi.data[[1]][[j]]$pt
          roi.index.k <-same.k.roi[same.k.roi!=j]
          # if (length(roi.index.z)!=0) {
          r <- unique (sapply (roi.index.k, function (k) {
            ptk <- struct$roi.data[[1]][[k]]$pt
            keep <- .pt.in.polygon (ptj[ ,1], ptj[ ,2],
                                      ptk[ ,1], ptk[ ,2]) > 0.5
            return (ifelse (any(keep), k,NA))}))
          r <- r[!is.na (r)]
          struct$roi.data[[1]][[j]]$level <- ifelse (length(r)!=0, length(r), 0)
        } 
      } else struct$roi.data[[1]][[same.k.roi]]$level <- 0
    }
  }
  
  if (external.only) {
    level.f  <- sapply(struct$roi.data[[1]], function(l) l$level==0)
    struct$roi.data[[1]] <- struct$roi.data[[1]][level.f]
  }
  # names (struct$roi.data) <- roi.name
  
  db <- .struct.moreinfo (struct$roi.data, struct$ref.from.contour, struct$thickness)
  struct$roi.info  <- cbind (struct$roi.info, db)
  if (alias=="") return(struct)
  return(.set.ref.obj(struct,list(vol)))
}
