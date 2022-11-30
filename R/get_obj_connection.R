#' List of connections between objects
#' @description The \code{get.obj.connection} function describes with 4 matrices 
#' the different connections between the DICOM objects of the patient.
#' @param pat "patient" class object, as loaded using \link[espadon]{load.patient.from.dicom}, 
#'  \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{toy.load.patient}. 
#' @return Returns a list of 4 named matrices:
#' \itemize{
#' \item the \code{$adjacency} matrix matrix specifies the source objects that 
#' generated the destination objects: the column names correspond to the 
#' destinations, and the row names to the sources.
#' \item the \code{$same.object} matrix specifies the elements belonging to the same 
#' DICOM object.
#' \item the \code{$components} matrix specifies the objects belonging to the same study.
#' \item the \code{$same.ref} matrix specifies the objects that share the same frame of 
#' reference, or with frames of reference linked in T.MAT (by a DICOM reg file 
#' for instance)
#' }
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (dxyz = c (5, 5, 5), beam.nb = 1)
#' get.obj.connection (patient)
#' display.obj.links (patient)
#' @seealso \link[espadon]{display.obj.links}
#' @export
get.obj.connection <- function (pat) {
  
  if (is.null(pat)) stop("pat is NULL")
  if (any(is.na(match(c("patient", "pat.pseudo","description","T.MAT"),names(pat))))) stop("pat is not a patient")
  obj.list <- obj.alias  <- obj.name <- c()
  obj.type <- obj.idx <- c()
  obj.con <- list ()
  end.idx <- length (pat)
  if (names(pat)[end.idx]=="dicom.dvh") end.idx <- end.idx - 1
  for (i in (which(names(pat)=="T.MAT")+1):end.idx) {
    for (j in 1:length (pat[[i]])) {
      obj.idx <- c(obj.idx, j)
      obj.list <- c (obj.list, paste(names (pat)[i], j))
      obj.alias <- c (obj.alias, pat[[i]][[j]]$object.alias)
      obj.name <- c (obj.name, pat[[i]][[j]]$object.name)
      if (is.null (pat[[i]][[j]]$ref.object.alias)) obj.con <- c(obj.con, "NA")
      else obj.con <- c(obj.con, list (pat[[i]][[j]]$ref.object.alias))
      names(obj.con) <- obj.alias
    }
  }
  
  M <- matrix (0, nrow=length (obj.list), ncol=length (obj.list))
  rownames (M) <- obj.list
  colnames (M) <- obj.list
  
  for (i in 1:length (obj.list)) {
    idx <- which (sapply (obj.con, function (V) obj.alias[i] %in% V))
    M[i, idx] <- 1
  }
 
  M.same.obj <- matrix (0, nrow=length (obj.list), ncol=length (obj.list))
  rownames (M.same.obj) <- obj.list
  colnames (M.same.obj) <- obj.list
  for (i in 1:length (obj.list)) {
    idx <- which (sapply (obj.name, function (V) obj.name[i] %in% V))
    M.same.obj[i, idx] <- 1
  }
  
  M_ <- (M+ diag(ncol(M)) + t(M)) == 1
  M_ [] <- sapply(M_,as.numeric)
  pere <- obj.list[apply(M,2,sum)==0]
  pere.idx <- match(pere, obj.list)
  M0 <- M_
  M0[] <- 0
  L <- lapply(1: length(pere), function(idx) {
    A <- M0
    A[pere[idx], ] <- M_[pere[idx], ]
    A[ ,pere[idx]] <- M_[pere[idx], ]
    A
  })
  
  for(idx in 1:length(L)){
    old.vect <-rep(FALSE, ncol(M_))
    new.vect <-  L[[idx]][pere.idx[idx],]>0
    while(!all(new.vect==old.vect)){
      old.vect <- L[[idx]][pere.idx[idx],]>0
      L[[idx]] <- L[[idx]] %*% M_
      new.vect <- L[[idx]][pere.idx[idx],]>0
    }
    L[[idx]] <- L[[idx]]>0
    L[[idx]][new.vect,] <- L[[idx]][rep(pere.idx[idx],sum(new.vect)),]
    L[[idx]] [] <- sapply(L[[idx]],as.numeric)
  }
  
  components <- M0
  for(idx in 1:length(L)) components <- components + L[[idx]]
  components <- (components>0) + 0

  same.ref <- M0
  for (idx in 1:length(pat$description.by.reg)){
    reg <- unlist(strsplit(pat$description.by.reg[[idx]]$object.alias,";"))
    ma <- match( obj.alias, reg)
    vect <- !is.na(ma) 
    same.ref[vect,] <- matrix (rep(vect,sum(vect)), ncol = length(vect), byrow=T)
  }
  
  return (list (adjacency = M, same.object = M.same.obj,
                components = components, same.ref=same.ref))
}