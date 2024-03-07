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
#' "DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT".
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object.
#' @param force.closed Boolean. Set to \code{TRUE} when the mesh represents the 
#' surface of a closed volume. 
#' @param verbose Boolean. If \code{TRUE} (default), a progress bar indicates 
#' the state of calculation.
#' @param ... Others parameters
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
#' @importFrom Rvcg vcgFaceNormals vcgBorder vcgIsolated checkFaceOrientation vcgGetEdge vcgNonBorderEdge
#' @importFrom Morpho mergeMeshes 
#' @importFrom rgl addNormals
#' @importFrom methods is
#' @export
struct.from.mesh <- function(mesh, z, thickness = NULL,
                             roi.name = mesh$object.alias,
                             roi.nb = 1,
                             roi.color = "#ff0000",
                             roi.type = "",
                             alias="", description = NULL,
                             force.closed = TRUE,
                             verbose = TRUE,...){
  
  
  args <- tryCatch(list(...), error = function(e)list())
  clean <- args[["clean"]]; if (is.null(clean)) clean <- TRUE; 
  
  
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
  ext.pt <- get.extreme.pt(mesh)
  obj <- obj.create(class="struct")
  obj$patient<-mesh$patient
  obj$patient.name<-mesh$patient.name
  obj$patient.bd<-mesh$patient.bd
  obj$patient.sex <- mesh$patient.sex
  obj$description <- description
  obj$frame.of.reference <- mesh$frame.of.reference
  obj$ref.pseudo <- mesh$ref.pseudo
  obj$nb.of.roi <- 0
  obj$thickness <- thickness
  obj$roi.info[1,] <- c(roi.nb[1], roi.name,"","AUTOMATIC",roi.color,
                        castlow.str(roi.name))
  obj$object.alias <- alias[1]
  obj$description <- description
  obj$roi.obs[1,] <- c(1,roi.nb[1],"","","","","",roi.type,"")
  
  nb.pt.max <- 65000
  
  z <- round(sort(unique(z)),6)
  dz <- max (max (abs (mesh$mesh$vb[3, mesh$mesh$it[1,]] - mesh$mesh$vb[3, mesh$mesh$it[2,]])),
             max (abs (mesh$mesh$vb[3, mesh$mesh$it[2,]] - mesh$mesh$vb[3, mesh$mesh$it[3,]])),
             max (abs (mesh$mesh$vb[3, mesh$mesh$it[3,]] - mesh$mesh$vb[3, mesh$mesh$it[1,]])))
  
  
  dxy <- max (max (abs (mesh$mesh$vb[1, mesh$mesh$it[1,]] - mesh$mesh$vb[1, mesh$mesh$it[2,]])),
              max (abs (mesh$mesh$vb[1, mesh$mesh$it[2,]] - mesh$mesh$vb[1, mesh$mesh$it[3,]])),
              max (abs (mesh$mesh$vb[1, mesh$mesh$it[3,]] - mesh$mesh$vb[1, mesh$mesh$it[1,]])),
              max (abs (mesh$mesh$vb[2, mesh$mesh$it[1,]] - mesh$mesh$vb[2, mesh$mesh$it[2,]])),
              max (abs (mesh$mesh$vb[2, mesh$mesh$it[2,]] - mesh$mesh$vb[2, mesh$mesh$it[3,]])),
              max (abs (mesh$mesh$vb[2, mesh$mesh$it[3,]] - mesh$mesh$vb[2, mesh$mesh$it[1,]]))
  )
  
  
  cut.range <- as.numeric(cut(z, c(z[ c(TRUE,(z - 2*dz)[-1] > (z + 2*dz)[-length(z)])]-dz*0.5, 
                                   z[length(z)]+0.5*dz), include.lowest = T))
  cut.range <- lapply(unique(cut.range), function(i) range(z[cut.range==i]))
  br <- unique(sort(c (seq(min(z) - 2*dz,max(z) + 2*dz, 4*dz), range(mesh$mesh$vb[3, ]) + c(-0.5,0.5))))
  
  
  if (is.null(mesh$mesh$normals)) mesh$mesh <- addNormals(mesh$mesh )
  M.L <- vcgIsolated(mesh$mesh, split=TRUE, silent =TRUE)
  
  keep <- sapply(M.L,function(l){
    rg <- range(l$vb[3,])
    any (sapply(cut.range, function(cr) length(c(rg[(rg >= cr[1]) & (rg <= cr[2])], cr[(cr >= rg[1]) & (cr <= rg[2])]))>0))
  })
  M.L <- M.L[keep]
  orientation <-  sapply(M.L,function(l) ifelse(checkFaceOrientation(l),1,-1))
  
  if (force.closed) {
    closed <- rep(TRUE,length(M.L))
  } else {
    closed <-  sapply(M.L,function(l) {all((vcgBorder(l)$border)==0)})
  }
  
  M.L <- list((M.L[orientation==1 & closed]), (M.L[orientation==-1 & closed]),
              (M.L[orientation==1 & !closed]), (M.L[orientation==-1 & !closed]))
  f <- sapply(M.L,function(l) length(l)>0)
  M.L <- M.L[f]
  orientation <- c(1, -1, 1, -1)
  closed <- c(TRUE, TRUE, FALSE, FALSE)
  orientation <- orientation[f]
  closed <- closed[f]
  M.L <- lapply(M.L,function(l) {if (length(l)>1) return(mergeMeshes(l)); return(l[[1]])})
  if (length(M.L)==0) return(obj)
  
  info.L <- list()
  M.idx <- 0
  Mesh.L <- list()
  as.mesh <- function(l){class(l)<- "mesh";return(l)}
  
  for(idx in 1:length(orientation)){
    z_ <- z
    histo <- hist(as.vector(M.L[[idx]]$vb[3, ]),  breaks = br, plot=FALSE)
    f <-histo$mids >= min(z_) - 2*dz & histo$mids <= max(z_) + 2*dz
    histo <- data.frame(counts = histo$counts[f], mids = histo$mids[f])
    while(nrow(histo)!=0){
      M.idx <- M.idx + 1
      f.jmax <- cumsum(histo$counts)<=nb.pt.max
      if (all(!f.jmax)) {jmax <- 1
      } else {jmax <- max(which(f.jmax))}
      info.L[[M.idx]] <- list()
      fz <- (histo$mids[1] - 2*dz <= z_ ) & (histo$mids[jmax] + 2*dz >= z_)
      info.L[[M.idx]]$rangez <- z_ [fz]
      info.L[[M.idx]]$orientation <- orientation[idx]
      info.L[[M.idx]]$closed <- closed[idx]
      Mesh.L[[M.idx]] <- nesting.cube(as.mesh(list(mesh = M.L[[idx]])), pt.min = c(ext.pt[1:2,1]- 2*dxy,min(info.L[[M.idx]]$rangez) - 2*dz),
                                      pt.max = c(ext.pt[1:2,2] + 2* dxy,max(info.L[[M.idx]]$rangez) + 2*dz))$mesh
      
      histo <- histo[-(1:jmax),]
      z_ <- z_ [!fz]
      
    }
  }
  rm(M.L)
  
  L.ct <- list()
  total <- length(unlist(lapply(info.L, function(i) i$rangez)))
  if (verbose)   pb <- progress_bar$new(format ="[:bar] :percent",total = total, width= 60)
  pb.sens <- FALSE
  if (length(Mesh.L)>0){
    for (M.idx in 1:length(Mesh.L)){
      
      pt <- round(t(Mesh.L[[M.idx]]$vb[1:3,]),6)
      
      edges_ <- cbind(vcgGetEdge(Mesh.L[[M.idx]], unique =TRUE)[,c("vert1", "vert2","border", "facept")],face2=NA,pbsens=FALSE,plan.z=NA,bump=1)
      colnames(edges_) <- c("vert1", "vert2","border", "face1","face2", "pbsens", "plan.z","bump")
      f <- edges_$border==1
      edges_ <- edges_[f, ]
      
      normals <- t(vcgFaceNormals(Mesh.L[[M.idx]]))
      # if (nrow(edges_)!=0) edges_$plan.z <-  round(normals[edges_$face1,3],6) == 1
      
      edges <- vcgNonBorderEdge(Mesh.L[[M.idx]],silent=TRUE)
      edges$pbsens <- FALSE
      v1.idx <- (Mesh.L[[M.idx]]$it[1,edges$face1]==edges$vert1) +
        2*(Mesh.L[[M.idx]]$it[2,edges$face1]==edges$vert1) +
        3*(Mesh.L[[M.idx]]$it[3,edges$face1]==edges$vert1)
      v2.idx <- (Mesh.L[[M.idx]]$it[1,edges$face1]==edges$vert2) + 
        2*(Mesh.L[[M.idx]]$it[2,edges$face1]==edges$vert2) + 
        3*(Mesh.L[[M.idx]]$it[3,edges$face1]==edges$vert2)
      
      edges$sens1 <- ((((v2.idx-v1.idx + 2) %% 3 ) == 0) * 2) - 1
      edges$sens2 <- (((((Mesh.L[[M.idx]]$it[1,edges$face2]==edges$vert2) + 
                           2*(Mesh.L[[M.idx]]$it[2,edges$face2]==edges$vert2) + 
                           3*(Mesh.L[[M.idx]]$it[3,edges$face2]==edges$vert2) -
                           (Mesh.L[[M.idx]]$it[1,edges$face2]==edges$vert1) -
                           2*(Mesh.L[[M.idx]]$it[2,edges$face2]==edges$vert1) -
                           3*(Mesh.L[[M.idx]]$it[3,edges$face2]==edges$vert1) + 2) %% 3 ) == 0) * 2 )- 1
      
      
      edges$pbsens <- !(edges$sens1 == -edges$sens2)
      pb.sens <- pb.sens | any(edges$pbsens)
      # f <- which(edges$sens1 == (!edges$sens2))
      # bary <- vcgBary(Mesh.L[[M.idx]]); f.i <-c(1774,1781,1782,1785)
      
      # for (fi in f.i) segments3d(rbind(bary[fi,],bary[fi,] + 0.2*normals[fi,]))
      edges$plan.z <- abs(round(normals[edges$face1,3],6)) == 1
      edges$plan2.z  <- abs(round(normals[edges$face2,3],6))== 1
      
      # f <- edges$plan.z & edges$plan2.z
      # edges <- edges[!f,]
      f <- (!edges$plan.z) & edges$plan2.z
      edges[f, c("face1","face2","sens1","sens2","plan.z","plan2.z")] <- edges[f, c("face2","face1","sens2","sens1","plan2.z","plan.z")]
      
      edges$plan2.z <- NULL
      
      
      
      vp <- cbind(
        normals[edges$face1,2] * normals[edges$face2,3] - normals[edges$face2,2] * normals[edges$face1,3],
        normals[edges$face1,3] * normals[edges$face2,1] - normals[edges$face2,3] * normals[edges$face1,1],
        normals[edges$face1,1] * normals[edges$face2,2] - normals[edges$face2,1] * normals[edges$face1,2]) * 
        (pt[edges$vert2,]- pt[edges$vert1,])
      edges$bump  <- sign(round(vp[,1] + vp[,2] + vp[,3],6)*edges$sens1 * info.L[[M.idx]]$orientation)  # * sign(normals[all.edges$face1[f],3]) < 0
      edges$sens1 <- edges$sens2 <- NULL
      edges$bump[!edges$plan.z] <- 1
      
      edges <- rbind(edges,edges_)
      edges <- edges[order(edges$vert1, edges$vert2), ]
      rm(edges_)
      
      v3.idx <- rbind((Mesh.L[[M.idx]]$it[,edges$face1][1,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face1][1,] != edges$vert2), 
                      (Mesh.L[[M.idx]]$it[,edges$face1][2,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face1][2,] != edges$vert2) ,
                      (Mesh.L[[M.idx]]$it[,edges$face1][3,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face1][3,] != edges$vert2))  * 
        Mesh.L[[M.idx]]$it[,edges$face1]
      edges$vert3.1 <- v3.idx[1,] + v3.idx[2,] + v3.idx[3,]
      v3.idx <- rbind((Mesh.L[[M.idx]]$it[,edges$face2][1,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face2][1,] != edges$vert2), 
                      (Mesh.L[[M.idx]]$it[,edges$face2][2,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face2][2,] != edges$vert2) ,
                      (Mesh.L[[M.idx]]$it[,edges$face2][3,] != edges$vert1  &  Mesh.L[[M.idx]]$it[,edges$face2][3,] != edges$vert2))  * 
        Mesh.L[[M.idx]]$it[,edges$face2]
      edges$vert3.2 <- v3.idx[1,] + v3.idx[2,] + v3.idx[3,]
      edges$vert3.2[is.na(edges$vert3.2)] <- +Inf
      
      vd <- round(pt[edges$vert2,] - pt[edges$vert1,],6)
      edges$plan.z <- abs(vd[,3]) == 0
      L.ct[[M.idx]] <- list()
      
      for (z_ in info.L[[M.idx]]$rangez){
        k <- rep(+Inf,nrow(vd))
        f <- vd[,3]!=0
        k[f] <- (z_-pt[edges$vert1[f],3])/vd[f,3] # pas dans le plan z_ étudié
        f <- vd[,3]==0 & pt[edges$vert1,3] == z_ # dans le plan z_ étudié
        k[f] <- 0
        
        # Le contour peut traverser une face :
        #    - par 2 arrêtes  (transition par face1, face2)
        #    - par un vertex et une arrête (transition par vertex vert1, face.opp)
        #    - par 2 vertex (=le long d'une arrete)
        
        #    - par 2 arrêtes
        
        k.flag <- k>=0 & k<=1  #& ((edges$bump > 0) | (open *edges$bump > 0))#&    (!all.edges$planz | (all.edges$planz & (all.edges$bump>0 | all.edges$border)))
        
        if (any(k.flag)){
          edges_ <- edges[k.flag, , drop=FALSE]
          
          row.names(edges_) <- NULL
          vb_ <- round(vd[k.flag, , drop=FALSE] * k[k.flag] +  pt[edges_$vert1, , drop=FALSE],6)
          
          f.planz <- edges_$plan.z
          
          vb <- as.data.frame(unique (rbind (vb_, 
                                             pt[edges_$vert1[f.planz], , drop=FALSE],
                                             pt[edges_$vert2[f.planz], , drop=FALSE])))
          colnames(vb) <- c("x","y","z")
          row.names(vb) <- NULL
          # wire3d(Mesh.L[[M.idx]], col ="black")
          # bary <- vcgBary(Mesh.L[[M.idx]])
          # text3d(bary,texts=1:nrow(bary), size=0.5)
          # text3d(t(Mesh.L[[M.idx]]$vb[1:3,]),texts=1:ncol(Mesh.L[[M.idx]]$vb),col="brown", size=0.5, adj=c(-0.2,-0.2,0))
          # for(idx in 1:nrow(edges_)) segments3d(pt[as.numeric(edges_[idx,1:2]),], col="red", lwd=2)
          # points3d(vb, col="green", size=5)
          # # text3d(vb,texts=1:nrow(vb),col="red", size=0.7)
          fact.pt <- vb$x + 1i*vb$y
          edges_$pt.nb <- match(vb_[,1] + 1i*vb_[,2],fact.pt)
          
          it.L <- list()
          
          #on traite les plans z 
          
          edges_z <- edges_[f.planz, ,drop=FALSE]
          edges_ <- edges_[!f.planz, ,drop=FALSE]
          
          it.L[[1]] <- as.data.frame(matrix(
            as.numeric(factor(c(pt[edges_z$vert1,1] + 1i* pt[edges_z$vert1,2],
                                pt[edges_z$vert2,1] + 1i* pt[edges_z$vert2,2]),
                              levels = fact.pt)),ncol=2))
          
          
          
          ptz <- unique(as.vector(it.L[[1]]))
          # ptz <- ptz[!is.na(ptz)]
          if (info.L[[M.idx]]$closed & !any(edges_z$pbsens)) it.L[[1]]<- it.L[[1]][(edges_z$bump > 0), , drop=FALSE]
          
          # it.L[[1]][is.na(it.L[[1]][,2]),2] <- +Inf
          f <- it.L[[1]][,1]>it.L[[1]][,2]
          it.L[[1]][f,1:2] <-  it.L[[1]][f,2:1]
          
          ## les autres points 
          row.names(edges_) <- edges_$vert1 + 1i*edges_$vert2
          dum <- edges_[,c("vert1", "vert3.1")]
          f <- dum$vert1 > dum$vert3.1
          dum[f,] <- dum[f,2:1]
          edges_$pt.nb1 <- edges_[as.character(dum[,1] + 1i*dum[,2]),]$pt.nb
          
          dum <- edges_[,c("vert1", "vert3.2")]
          f <- dum$vert1 > dum$vert3.2
          dum[f,] <- dum[f,2:1]
          edges_$pt.nb2 <- edges_[as.character(dum[,1] + 1i*dum[,2]),]$pt.nb
          
          dum <- edges_[,c("vert2", "vert3.1")]
          f <- dum$vert2 > dum$vert3.1
          dum[f,] <- dum[f,2:1]
          edges_$pt.nb3 <- edges_[as.character(dum[,1] + 1i*dum[,2]),]$pt.nb
          
          dum <- edges_[,c("vert2", "vert3.2")]
          f <- dum$vert2 > dum$vert3.2
          dum[f,] <- dum[f,2:1]
          edges_$pt.nb4 <- edges_[as.character(dum[,1] + 1i*dum[,2]),]$pt.nb
          
          # on enlève tous les ptz dans pt.nb
          edges_ <- edges_[is.na(match(edges_$pt.nb,ptz)), ]
          
          
          
          it.L[[2]] <- as.data.frame(rbind(as.matrix(edges_[,c("pt.nb","pt.nb1")]), 
                                           as.matrix(edges_[,c("pt.nb","pt.nb2"),]),
                                           as.matrix(edges_[,c("pt.nb","pt.nb3")]), 
                                           as.matrix(edges_[,c("pt.nb","pt.nb4"),])))
          it.L[[2]][is.na(it.L[[2]][,2]),2] <- +Inf 
          f <- it.L[[2]][,1] > it.L[[2]][,2]
          it.L[[2]][f,] <- it.L[[2]][f,2:1, drop=FALSE]
          it.L[[2]] <-  it.L[[2]][!duplicated(it.L[[2]]),, drop=FALSE]
          it.L[[2]] <-  it.L[[2]][ it.L[[2]][,1]!=  it.L[[2]][,2] ,, drop=FALSE]
          colnames(it.L[[2]] ) <- colnames(it.L[[1]] )
          # on cherche le points uniques --
          f <- is.infinite(it.L[[2]][,2])
          vb_u <- it.L[[2]][f,1]
          it.L[[2]] <- it.L[[2]][!f, , drop=FALSE]
          vb_u <- vb_u[is.na(match(vb_u,c(it.L[[2]][,1],it.L[[2]][,2])))]
          # points3d(vb[vb_u,], col="red", size=20)
          
          #et on les mets dans un contours
          L.ct_ <- lapply(vb_u,function(pt.idx) 
            list(type="POINT",pt = vb[pt.idx, c("x","y","z")], level =0))
          
          # on rassemble les connexions en cas de contour étudié comme fermé
          # if (info.L[[M.idx]]$closed) {
          it.L[[1]] <- rbind(it.L[[1]],it.L[[2]])
          it.L[[2]]  <- NULL
          f <- it.L[[1]][,1]>it.L[[1]][,2]
          it.L[[1]][f,1:2] <-  it.L[[1]][f,2:1]
          it.L[[1]] <- unique(it.L[[1]])
          row.names( it.L[[1]]) <-NULL
          # }
          
          # On traite lescontours
          for(it in it.L){
            if (nrow(it)>0){
              h <- hist(c(it[,1],it[,2]),seq(0.5,nrow(vb)+0.5,1), plot=FALSE)
              occ <- matrix(c(h$counts[it[,1]],
                              h$counts[it[,2]]),byrow=FALSE, ncol=2)
              f <- occ[,2] ==1
              it[f,1:2] <- it[f,2:1]
              occ[f,1:2] <- occ[f,2:1]
              ord <-order(occ[,1],vb[it[,1],1],vb[it[,1],2])
              it <- it[ord,]
              occ <- occ[ord,]
              
              it[is.infinite(it[,2]),2] <- -1
              tab <-
                as.data.frame(matrix(.isoclineC(as.numeric(it[,1]),as.numeric(it[,2])),
                                     byrow=TRUE,ncol=4,
                                     dimnames = list(NULL,
                                                     c("dispo.it","dispo.pt1","dispo.pt2", "swap"))))
              f <- tab$swap==1
              it[f,1:2] <- it[f,2:1]
              l.idx <- max(tab$dispo.it)
              
              
              L.ct_ <- c(L.ct_, lapply((0:l.idx)[-1], function(idx){
                f <- tab$dispo.it==idx
                pt <- data.frame(ord=c(tab$dispo.pt1[f],tab$dispo.pt2[f]),pt
                                 =c(it[f,1],it[f,2]))
                pt <- pt[!duplicated(pt$ord),]
                pt  <- pt[order(pt$ord),]
                
                op <- (pt$pt[1]!=pt$pt[nrow(pt)])
                pt <- vb[pt$pt,c("x","y","z")]
                if (clean & nrow(pt)>1) {
                  pt <- pt[c(TRUE, !(round(diff(pt[,1]),6)==0 & round(diff(pt[,2]),6)==0)),]
                  if(nrow(pt)>2) {
                    slope.x <- diff(pt[,1])
                    slope.y <- diff(pt[,2])
                    d <-  sqrt(slope.x^2 + slope.y^2)
                    pt <- pt[c(TRUE, !(diff(round(slope.x/d,6))==0 & diff(round(slope.y/d,6))==0), TRUE),]}
                }
                
                row.names(pt) <- NULL
                if (op) return(list(type="OPEN_PLANAR", pt = pt, level =0))
                return(list(type="CLOSED_PLANAR", pt = pt, level =0))
              }))
              
              # L.ct_ <- L.ct_[!sapply(L.ct_, is.null)]
              # length(sapply(L.ct_,function(l) l$type))
              if (length(L.ct_)>0){
                same.k.roi <- 1:length(L.ct_)
                if (length(same.k.roi)>1) {
                  same.k.roi <-same.k.roi[-1]
                  for (j in same.k.roi){
                    ptj<- L.ct_[[j]]$pt
                    roi.index.k <-same.k.roi[same.k.roi!=j]
                    
                    r <- unique (sapply (roi.index.k, function (k) {
                      ptk <- L.ct_[[k]]$pt
                      if (L.ct_[[k]]$type != "CLOSED_PLANAR") return(NA)
                      keep <- .pt.in.polygon (ptj[,1], ptj[,2], ptk[ ,1],
                                                        ptk[ ,2]) > 0.5
                      return (ifelse (any(keep), k,NA))}))
                    r <- r[!is.na (r)]
                    L.ct_[[j]]$level <- length(r)
                  }
                } else L.ct_[[same.k.roi]]$level <- 0
              }
              L.ct[[M.idx]] <-  c(L.ct[[M.idx]],L.ct_)
            }
          }
        }
        if (verbose)  pb$tick() 
      }
    }
  }
  obj$roi.data[[1]] <- do.call("c", L.ct)
  zm.order <- order(sapply(obj$roi.data[[1]],function(l) l$pt[1,3]))
  obj$roi.data[[1]] <- obj$roi.data[[1]][zm.order]
  obj$nb.of.roi <- 1
  obj$roi.info <- cbind(obj$roi.info,
                        .struct.moreinfo (obj$roi.data, obj$ref.from.contour, obj$thickness))
  # display.3D.contour(obj, roi.lwd = 3)
  if (pb.sens) warning("Some faces have an inconsistent orientation. Some contours will be considered open.")
  return(obj)
}