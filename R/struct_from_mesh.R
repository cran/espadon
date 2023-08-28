#' Creation of struct class object from an espadon mesh
#' @description The \code{struct.from.mesh} function creates a struct object with 
#' a unique RoI, defined by the contours of a mesh.
#' @param mesh  espadon mesh class object.
#' @param z z-coordinate vector where mesh contours are computed.
#' @param thickness struct thickness between 2 adjacent contours. If NULL (default),
#' it is deduced from \code{z}.
#' @param roi.name Character string, representing the name of created RoI.
#' @param roi.nb Positive integer, representing the number of created RoI.
#' @param roi.color Color of the created RoI, in hex code format ("#RRGGBB").
#' @param roi.type Type of RoI, from among "", "EXTERNAL", "PTV", "CTV", "GTV", 
#' "TREATED_VOLUME", "IRRAD_VOLUME", "OAR", "BOLUS", "AVOIDANCE", "ORGAN", "MARKER", 
#' "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", "CAVITY", "BRACHY_CHANNEL", 
#' "BRACHY_ACCESSORY", "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", "SUPPORT", "FIXATION", 
#' "DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT"
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' @return Returns a "struct" class object (see \link[espadon]{espadon.class}
#' for class definition), including the unique \code{roi.name} as region of interest.

#' @examples
#' # Creation of an espadon mesh of a cube
#' M <- obj.create (class = "mesh")
#' M$mesh <- Rvcg::vcgIsotropicRemeshing (Rvcg::vcgBox(),0.5) 
#' M$nb.faces <- ncol (M$mesh$it)
#'
#' S <- struct.from.mesh (M, z = seq(-1,1,0.5))
#' display.3D.contour(S)
#' @importFrom Rvcg vcgFaceNormals
#' @importFrom methods is
#' @export
struct.from.mesh <- function(mesh, z, thickness = NULL, 
                             roi.name = mesh$object.alias,
                             roi.nb = 1,
                             roi.color = "#ff0000", 
                             roi.type = c ("","EXTERNAL", "PTV", "CTV", "GTV", 
                                           "TREATED_VOLUME", "IRRAD_VOLUME", "OAR", 
                                           "BOLUS", "AVOIDANCE", "ORGAN", "MARKER", 
                                           "REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", 
                                           "CAVITY", "BRACHY_CHANNEL", "BRACHY_ACCESSORY", 
                                           "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", 
                                           "SUPPORT", "FIXATION", "DOSE_REGION", 
                                           "CONTROL", "DOSE_MEASUREMENT"),
                        
                             alias="", description = NULL){
  
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
  if (roi.name=="") roi.name <- "ROI_from_mesh"
  if (is.null(description)) description <- paste("ROI from mesh",mesh$object.alias)
  if (is.null(thickness)) {
    thickness <- 0
    if (length(z)>1) thickness <- min(diff(sort(z)))
  } 
  
  obj <- obj.create(class="struct")
  obj$patient<-mesh$patient
  obj$patient.name<-mesh$patient.name
  obj$patient.bd<-mesh$patient.bd
  obj$patient.sex <- mesh$patient.sex
  obj$description <- description
  obj$frame.of.reference <- mesh$frame.of.reference
  obj$ref.pseudo <- mesh$ref.pseudo
  obj$nb.of.roi <- 1
  obj$thickness <- thickness
  obj$roi.info[1,] <- c(roi.nb[1], roi.name,"","AUTOMATIC",roi.color, castlow.str(roi.name))
  obj$object.alias <- alias[1]
  obj$description <- description
  obj$roi.obs[1,] <- c(1,roi.nb[1],roi.name,"","","","",roi.type,"")
  
  pt <- t(mesh$mesh$vb[1:3,])
  z <- sort(unique(z))
  # all.edges <- vcgGetEdge(mesh$mesh,unique=FALSE)
  it <- cbind(t(mesh$mesh$it), 1:ncol(mesh$mesh$it))
  all.edges <- as.data.frame(rbind(it[,c(1,2,4)],it[,c(2,3,4)],it[,c(3,1,4)]))
  colnames(all.edges) <- c("vert1","vert2","facept")
  vd <- pt[all.edges$vert2,] - pt[all.edges$vert1,]
  L.ct <- list()
  #repérer les points vert1 qui sont au bord
  all.edges$border <-is.na(match(paste(all.edges$vert1, all.edges$vert2),paste(all.edges$vert2, all.edges$vert1)))
  # all.edges$border  <-  vcgBorder(mesh$mesh)$bordervb[all.edges$vert1]
  planz <- abs(abs(vcgFaceNormals(mesh$mesh)[3,]) - 1) <1e-6
  all.edges$planz <- planz[all.edges$facept]  
  
  for (z_ in z){
    first.in.plane <- TRUE
    k <- rep(+Inf,nrow(vd))
    f <- vd[,3]!=0
    k[f] <- (z_-pt[all.edges$vert1[f],3])/vd[f,3]
    f <- vd[,3]==0 & pt[all.edges$vert1,3] == z_
    k[f] <- 0
    
    #on élimine tout ceux qui sont dans le plan z qui qui ne sont pas des bords
    k.flag <- k>=0 & k<=1 & 
      (!all.edges$planz | (all.edges$planz & all.edges$border))
    
    if (any(k.flag)){ 
      edges_ <- all.edges[k.flag, ]
      f <-edges_$border & edges_$planz
      le <- sum(f)
      ptz1 <- data.frame(x=rep(NA,le),y=rep(NA,le),z=rep(NA,le),ch=rep(NA,le))
      ptz2 <- data.frame(x=rep(NA,le),y=rep(NA,le),z=rep(NA,le),ch=rep(NA,le))
      if (le>0){
        ptz1[,c("x","y","z")] <- round(pt[edges_$vert1[f],],4)
        ptz1$ch <- paste(ptz1[,1],ptz1[,2])
        ptz2[,c("x","y","z")] <- round(pt[edges_$vert2[f],],4)
        ptz2$ch <- paste(ptz2[,1],ptz2[,2])
      }
      ptz_ <- round(as.data.frame(vd[k.flag,][!f,]*k[k.flag][!f] + 
                                    pt[edges_$vert1,][!f,]),4)
      if (nrow(ptz_)>0){
        colnames(ptz_) <- c("x","y","z")
        ptz_$face <- edges_$facept[!f]
        ptz_$planz<- edges_$planz[!f]
        ptz_$border <- edges_$border[!f]
        
        #on élimine les points dupliqués dans la même face
        ptz_$ch <- paste(ptz_[,1],ptz_[,2])
        f <- !duplicated(ptz_)
        # f <- !duplicated(paste(ptz_$ch , ptz_$face))
        ptz_ <- ptz_[f,]
        ptz_ <- ptz_[order(ptz_$x,ptz_$y),]
        
        face <- unique(ptz_$face)
        m <- match(face,ptz_$face)
        dum <-ptz_[m,c("x","y","z","ch")]
        # dum <-ptz_[m,c("x","y","z")] 
        ptz1 <- rbind(ptz1,dum)
        
        ptz_ <- ptz_[-m, ]
        m <- match(face,ptz_$face)
        f <- !is.na(m)
        dum[] <- NA
        dum[f,] <- ptz_[m[f],c("x","y","z","ch")]
        # dum[f,] <- ptz_[m[f],c("x","y","z")]
        ptz2 <- rbind(ptz2,dum)
      }
      
      vb <- unique(rbind(ptz1[,c("x","y","z")],ptz2[!is.na(ptz2$ch),c("x","y","z")]))
      vb <- vb[order(vb$x,vb$y), ]
      
      it <- matrix(nrow(ptz1)+1, nrow=nrow(ptz1), ncol =2)
      it[,1] <- match(paste(ptz1$x,ptz1$y),paste(vb$x,vb$y))
      it[,2] <- match(paste(ptz2$x,ptz2$y),paste(vb$x,vb$y))
      
      
      #on trie pour que it[,1]<it[,2] 
      f <- it[,1]>it[,2]  & !is.na(it[,2])
      it[f,1:2] <- it[f,2:1]
      
      
      # on nettoie les points dupliqué (sur des faces différentes) 
      
      it <- matrix(it[!duplicated(it),],ncol=2)
      
      
      # on cherche le points uniques --
      f <- is.na(it[,2])
      vb_u <- it[f,1]
      it <- it[!f,]
      vb_u <- vb_u[is.na(match(vb_u,c(it[,1],it[,2])))]
      
      
      #on traite les  faces qui n'ont qu'un point d'occurence 1
      L.ct_ <- lapply(vb_u,function(pt.idx) list(type="POINT",pt = vb[pt.idx, ],
                                                 level =0, A = 0 , G=c(vb[pt.idx, ],1)))
      
      # On traite lescontours
      if (nrow(it)>0){
        h <- hist(c(it[,1],it[,2]),seq(0.5,nrow(vb)+0.5,1), plot=FALSE)
        occ <- matrix(c(h$counts[it[,1]],  h$counts[it[,2]]),byrow=FALSE, ncol=2)
        f <- occ[,2] ==1
        it[f,1:2] <- it[f,2:1]
        occ[f,1:2] <- occ[f,2:1]
        ord <-order(occ[,1],vb[it[,1],1],vb[it[,1],2])
        it <- it[ord,]
        occ <- occ[ord,]
        
        # l.idx0 <- l.idx 
        
        
        
        it[is.na(it[,2]),2] <- -1
        tab <- as.data.frame(matrix(.isoclineC(as.numeric(it[,1]),as.numeric(it[,2])),
                                    byrow=TRUE,ncol=4,
                                    dimnames = list(NULL, c("dispo.it","dispo.pt1","dispo.pt2", "swap"))))
        f <- tab$swap==1
        it[f,1:2] <- it[f,2:1]
        l.idx <- max(tab$dispo.it)
        
        
        L.ct_ <- c(L.ct_, lapply((0:l.idx)[-1], function(idx){
          f <- tab$dispo.it==idx
          pt <- data.frame(ord=c(tab$dispo.pt1[f],tab$dispo.pt2[f]),pt =c(it[f,1],it[f,2]))
          pt <- pt[!duplicated(pt$ord),]
          pt  <- pt[order(pt$ord),]
          le <- nrow(pt)
          op <- (pt$pt[1]!=pt$pt[le])
          pt <- vb[pt$pt,]
          row.names(pt) <- NULL
          
          if (op) return(list(type="OPEN_PLANAR", pt = pt, level =0, A=0, G =c(NA,NA,NA,1)))
          i <- 2:le
          A <- sum (pt[i-1,1] * pt[i,2] - pt[i,1] * pt[i-1,2]) / 2
          G <- c (sum ((pt[i-1,1] + pt[i,1]) * (pt[i-1,1] * pt[i,2] - pt[i,1] * pt[i-1,2])) / (6 * A),
                  sum ((pt[i-1,2] + pt[i,2]) * (pt[i-1,1] * pt[i,2] - pt[i,1] * pt[i-1,2])) / (6 * A),
                  pt[1,3], 1 )
          
          return(list(type="CLOSED_PLANAR", pt = pt, level =0, A = A))
        }))
      }
      
      L.ct_ <- L.ct_[order(sapply(L.ct_,function(l) l$A))]
      same.k.roi <- (0:length(L.ct_))[-1]
      if (length(same.k.roi)>1) {
        for (j in same.k.roi){
          ptj<- L.ct_[[j]]$pt
          roi.index.k <-same.k.roi[same.k.roi!=j]
          # if (length(roi.index.z)!=0) {
          r <- unique (sapply (roi.index.k, function (k) {
            ptk <- L.ct_[[k]]$pt
            if (L.ct_[[k]]$type != "CLOSED_PLANAR") return(NA)
            keep <- .pt.in.polygon (ptj[1 ,1], ptj[1 ,2], ptk[ ,1], ptk[ ,2]) > 0.5
            return (ifelse (any(keep), k,NA))}))
          r <- r[!is.na (r)]
          L.ct_[[j]]$level <- ifelse (length(r)!=0, length(r), 0)
        } 
      } else L.ct_[[same.k.roi]]$level <- 0
      L.ct <-  c(L.ct,L.ct_)
      
    }
  }
  
  obj$roi.data[[1]] <- L.ct
  obj$roi.info <- cbind(obj$roi.info, 
                        .struct.moreinfo (obj$roi.data,obj$ref.from.contour,obj$thickness))
  return(obj)
}
