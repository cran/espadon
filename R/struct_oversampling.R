struct.oversampling <- function (struct, fact.z = 2, 
                                 roi.name = NULL, roi.sname = NULL, roi.idx = NULL,
                                 vol.ref = NULL, alias = "", 
                                 description = NULL, verbose = TRUE, ...){
  args <- tryCatch(list(...), error = function(e)list())
  args <- args[is.na(match(names(args),
                           c("smooth.iteration", "smooth.type", "smooth.lambda", 
                             "smooth.mu", "smooth.delta")))]
  
  if (!is (struct, "struct")) {
    stop ("struct should be a struct class object.")
  }
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (length(roi.idx)==0) stop ("no roi selected.")
  struct_  <- obj.create (class = "struct", alias = alias)
  if (is.null(description)) {struct_$description <- paste0 ("oversampling /",fact.z," of ", struct$object.alias)
  } else {struct_$description <- description}
  
  struct_$patient <- struct$patient
  struct_$patient.name <- struct$patient.name
  struct_$patient.bd <- struct$patient.bd
  struct_$patient.sex <- struct$patient.sex
  struct_$frame.of.reference <- struct$frame.of.reference
  struct_$ref.pseudo <- struct$ref.pseudo
  struct_$thickness <- struct$thickness / fact.z
  struct_$nb.of.roi<- length(roi.idx)
  struct_$roi.info <- struct$roi.info[roi.idx,]
  
  row.names(struct_$roi.info) <- NULL
  struct_$roi.obs <- struct$roi.obs[roi.idx,]
  row.names(struct_$roi.obs) <- NULL
  struct_$roi.data <- list()
  roi.z <- sort (unique (unlist (lapply (struct$roi.data[roi.idx], function(roi.l) unique (sapply(roi.l,function(l) l$pt[1,3]))))))
  M <- solve(struct$ref.from.contour)
  
  if (is.null(vol.ref)){
    mid.z <-roi.z[length(roi.z)/2]
    ext.pt <- get.extreme.pt(struct,roi.idx = roi.idx)
    # on le met dans le référentiel du contour :
    ext.pt <- M %*% as.matrix(rbind(ext.pt,c(1,1)))
    mid.pt <- round(apply(ext.pt[1:3,], 1, mean),3)
    mid.pt[3] <- mid.z
    dxyz <- struct$thickness * c(0.5,0.5,1)
    n.ijk <- ceiling(2*round(ext.pt[1:3,2] - mid.pt,3)/dxyz)+8
    f <- n.ijk%%2 == 0 ## on veut un nombre impair pour que le centre tombe sur un nombre impair
    n.ijk[f] <- n.ijk[f]+1
    vol.ref <- vol.create(n.ijk,dxyz,mid.pt, ref.pseudo = struct$ref.pseudo,
                          frame.of.reference = struct$frame.of.reference)
    vol.ref$orientation <- c(struct$ref.from.contour[1:3,1],struct$ref.from.contour[1:3,2])
    vol.ref$xyz.from.ijk <- struct$ref.from.contour %*% vol.ref$xyz.from.ijk
    vol.ref$xyz0 <- (cbind(vol.ref$xyz0,1) %*% t(struct$ref.from.contour))[,1:3]
  } else{
    if (!is (vol.ref, "volume")) {
      stop ("vol.ref should be a volume class object.")
    }
    if ((vol.ref$ref.pseudo!= struct$ref.pseudo) & 
        !all(vol.ref$orientation==c(struct$ref.from.contour[1:3,1],struct$ref.from.contour[1:3,2]))) {
      stop ("struct and vol.ref do not have the same orientation")
    }
    
  }
  
  for (j in 1:length(roi.idx)){
    roi.data <- struct$roi.data[[roi.idx[j]]]
    back <- nesting.roi(vol.ref, struct, roi.idx = roi.idx[j], xyz.margin = abs(vol.ref$dxyz)*2)
    args[["bin"]]  <- bin.from.roi(back, struct, roi.idx = roi.idx[j],verbose)
   
    # changement de referentiel
    mesh <- do.call(mesh.from.bin,args)
    mesh$mesh$vb <- M %*% mesh$mesh$vb
    mesh$mesh$normals <-  M %*% mesh$mesh$normals

    ctz<-sapply(roi.data, function(l) l$pt[1,3])
    uctz <- sort(unique(ctz))
    z.study <- seq(uctz[1], uctz[length(uctz)], struct$thickness/fact.z)
    z.study <- z.study[-seq(1, length(z.study), struct$thickness)]
    s <- struct.from.mesh(mesh, z.study, verbose)
    s$roi.data[[1]] <- lapply (1:length (s$roi.data[[1]]), function(i){
      if (nrow(s$roi.data[[1]][[i]]$pt)<3) return (NULL)
      pt <- polyg.simplify(s$roi.data[[1]][[i]]$pt, tol = 0.1)
      if ((castlow.str(s$roi.data[[1]][[i]]$type) == "closedplanar") & (polyg.area (pt) == 0)) return (NULL)
      return (list(type = s$roi.data[[1]][[i]]$type, pt = pt, level = s$roi.data[[1]][[i]]$level))
    })
    s$roi.data[[1]] <- s$roi.data[[1]][!sapply(s$roi.data[[1]], is.null)]
    roi.data <- c(roi.data, s$roi.data[[1]])
    
    struct_$roi.data[[j]] <- roi.data[order(sapply(roi.data , function(l) l$pt[1,3]))]
  }
  struct_$roi.info[,7:17] <- .struct.moreinfo (struct_$roi.data, struct_$ref.from.contour, struct_$thickness)
  return(struct_)
}