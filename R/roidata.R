roidata.xor.to.keyhole <- function(roi.data, thickness = NULL) {
  if (length(roi.data)==0) return(roi.data)
  level <- sapply(roi.data, function(li)  li$level)
  # if (all(level==0)) return(roi.data)
  
  roi.all.z<- sapply(roi.data, function(li)  li$pt[1,3])
  type  <- sapply(roi.data, function(li)  castlow.str(li$type))
  if (is.null(thickness)){
    thickness <- 0  
    z.diff <- diff(sort(unique(roi.all.z)))
    if (length(z.diff)>0) {
      thickness <- unlist(z.diff)
      thickness <- tryCatch(round(median(thickness [thickness >1e-2]),3), error = function (e) 0)
    }
  }
  
  kz <- rep(0,length(roi.all.z))
  if (thickness>0) kz <- round((roi.all.z -roi.all.z[1])/thickness)
  
  
  made <- rep(FALSE, length(roi.data))
  
  for(slevel in rev(sort(unique(level)))){
    idx.level = which(level==slevel & type =="closedplanar")
    if(slevel%%2 !=0){
      for(ct.idx in idx.level){
        # cat(ct.idx,"\n")
        potentiel.parent <- which(kz==kz[ct.idx] &  level == slevel-1 & type =="closedplanar")
        parent <- potentiel.parent[sapply(potentiel.parent, function(parent.idx)
          any(.pt.in.polygon(roi.data[[ct.idx]]$pt[,1],roi.data[[ct.idx]]$pt[,2],
                                       roi.data[[parent.idx]]$pt[,1],
                                       roi.data[[parent.idx]]$pt[,2], 1e-20) == 1))]
        
        roi.data[[ct.idx]]$pt <- polyg.sort (roi.data[[ct.idx]]$pt, clockwise =  FALSE)
        roi.data[[parent]]$pt <- polyg.sort (roi.data[[parent]]$pt, clockwise = TRUE)
        le.ct <- nrow(roi.data[[ct.idx]]$pt)
        le.p <- nrow(roi.data[[parent]]$pt)
        
        #on met en premier le plus petit en y
        
        pt.link <- pt.nearest(roi.data[[ct.idx]]$pt[1:(le.ct -1),],
                              roi.data[[parent]]$pt[1:(le.p -1),], full.info = TRUE)
        
        pt.link <- pt.link$pt.index
        
        roi.data[[ct.idx]]$pt <- rbind(roi.data[[ct.idx]]$pt[ pt.link[1]:(le.ct-1),],
                                       roi.data[[ct.idx]]$pt[-(pt.link[1]:le.ct),],
                                       roi.data[[ct.idx]]$pt[ pt.link[1],])
        roi.data[[parent]]$pt <- rbind(roi.data[[parent]]$pt[1:pt.link[2],],
                                       roi.data[[ct.idx]]$pt,
                                       roi.data[[parent]]$pt[pt.link[2]:le.p,])#,
        # roi.data[[parent]]$pt[1,])
        row.names(roi.data[[parent]]$pt) <- NULL 
        row.names(roi.data[[ct.idx]]$pt) <- NULL 
        made[c(parent,ct.idx)] <- TRUE
        
        
        roi.data[[ct.idx]] <- list()                               
      }
    } else {
      for(ct.idx in idx.level) {
        #on met les level pair à 0 s'ils sont different de 0
        roi.data[[ct.idx]]$level <- 0
        if (type[ct.idx] =="closedplanar" & !made[ct.idx]) {
          roi.data[[ct.idx]]$pt <- polyg.sort (roi.data[[ct.idx]]$pt, clockwise = TRUE)
        }
        
      }    
    }
  }
  
  roi.data <- roi.data[sapply(roi.data, function(li) length(li)!=0)] 
  return(roi.data)
}

###################################################################################
roidata.keyhole.to.xor <- function(roi.data, thickness = NULL, clockwise = TRUE){
  if (length(roi.data)==0) return(roi.data)
  level <- sapply(roi.data, function(li)  li$level)
  
  roi.all.z<- sapply(roi.data, function(li)  li$pt[1,3])
  
  new.roi.data <- list()
  for(i in 1:length(roi.data)){
    roi.data[[i]]$pt <- round(roi.data[[i]]$pt, 6)
    if (castlow.str(roi.data[[i]]$type)=="closedplanar"){
      pt <- roi.data[[i]]$pt
      le <- nrow(pt)
      while (all(pt[le,]==pt[1,])){pt <- pt[-le,]; le <- le-1}
      end.idx <- duplicated.data.frame(pt)
      
      
      while (any(end.idx)){
        end.idx <- which(end.idx)[1]
        
        first.idx <- sapply(end.idx,function(j) 
          which(pt[,1]==pt[j,1] & pt[,2]==pt[j,2] & (1:le)!=j))
        if (end.idx == first.idx + 1){ 
          pt <- pt[-end.idx,]
        } else {
          new.roi.data <- c(list(list(type = roi.data[[i]]$type, 
                                      pt = pt[first.idx:end.idx,], level=0)),
                            new.roi.data)
          if (all(pt[first.idx-1,] == pt[end.idx+1,])){
            pt <-  pt[c(1:(first.idx-2),(end.idx+1):le),]
          } else {
            pt <- pt[c(1:first.idx,(end.idx+1):le),] 
          }
        }
        le <- nrow(pt)
        row.names(pt) <- NULL
        
        end.idx <- duplicated.data.frame(pt)
      }
      new.roi.data <- c(list(list(type = roi.data[[i]]$type, pt = pt[c(1:le,1),], level=0)),
                        new.roi.data)
    }
    
  }
  
  return(roidata.sort(roi.data = new.roi.data, thickness = thickness, clockwise = clockwise))
  
}

################################################################################
roidata.sort <- function(roi.data, thickness = NULL, clockwise = TRUE){
  if (length(roi.data)==0) return(roi.data)
  roi.all.z<- sapply(roi.data, function(li)  li$pt[1,3])
  
  if (is.null(thickness)){
    thickness <- 0  
    z.diff <- diff(unique(roi.all.z))
    if (length(z.diff)>0) {
      thickness <- unlist(z.diff)
      thickness <- tryCatch(round(median(thickness [thickness >1e-2]),3), error = function (e) 0)
    }
  }
  
  ord <- order(roi.all.z,decreasing = sign(thickness)==0)
  roi.data <- roi.data[ord]
  roi.all.z <- roi.all.z[ord]
  kz <- rep(0,length(roi.all.z))
  
  if (thickness!=0) kz <- round((roi.all.z -roi.all.z[1])/thickness)
  type  <- sapply(roi.data, function(li)  castlow.str(li$type))
  
  for (k.value in unique(kz)){
    same.k.roi <- which (kz ==k.value)
    
    if (length(same.k.roi)>1) {
      for (j in same.k.roi){
        ptj<- roi.data[[j]]$pt
        roi.index.k <-same.k.roi[same.k.roi!=j]
        r <- unique (sapply (roi.index.k, function (k) {
          if ( type[k] != "closedplanar") return(NA)
          ptk <- roi.data[[k]]$pt
          keep <- .pt.in.polygon (ptj[ ,1], ptj[ ,2], ptk[ ,1], ptk[ ,2], 1e-20) ==1 
          return (ifelse (any(keep), k,NA))}))
        r <- r[!is.na (r)]
        
        roi.data[[j]]$level <- ifelse (length(r)!=0, length(r), 0)
        
        if (type[j] == "closedplanar") 
          # pt <- polyg.sort(l$pt)
          
          roi.data[[j]]$pt  <- polyg.sort (roi.data[[j]]$pt, 
                                              clockwise = (roi.data[[j]]$level %%2) == as.numeric(!clockwise))
      } 
      

      
      A <- sapply(same.k.roi,function(j) polyg.area(roi.data[[j]]$pt))
      level <- sapply(same.k.roi,function(j) roi.data[[j]]$level)
      Y <-sapply(same.k.roi,function(j) roi.data[[j]]$pt[1,2]) 
      X <-sapply(same.k.roi,function(j) roi.data[[j]]$pt[1,1]) 
      ord <- same.k.roi[order(level,max(A)-A,Y,X)]
      roi.data[same.k.roi] <- roi.data[ord]

    } else {
      roi.data[[same.k.roi]]$level <- 0
      roi.data[[same.k.roi]]$pt  <- polyg.sort (roi.data[[same.k.roi]]$pt, 
                                                   clockwise = clockwise)
    }
  }
  return(roi.data)
  
}

##################################################################################
pt.nearest <- function(pt1, pt2, full.info = FALSE){
  
  if (length(pt1)==0) return(numeric(0))
  result <- rep(NA,length(pt1))
  if (length(pt2)==0) return(result)
  res <- .ptnearestC(pt1[,1],pt1[,2],pt1[,3],pt2[,1],pt2[,2],pt2[,3], full.info)
  return(res)
}