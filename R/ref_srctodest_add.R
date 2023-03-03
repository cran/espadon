#' Linking two existing frames of reference in T.MAT
#' @description The \code{ref.srctodest.add} function links the source frame of 
#' reference with the destination frame of reference. 
#' @param src.ref Character string, pseudonym of the source frame of reference.
#' @param dest.ref Character string, pseudonym of the destination frame of reference.
#' @param TM 4x4 tansfert matrix for moving from \code{src.ref} to {dest.ref}.
#' @param T.MAT "t.mat" class object created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, then only the link 
#' between \code{src.ref} and \code{dest.ref} is established.

#' @return Returns a "t.mat" class object, which contains the transfer 
#' matrices from or to \code{dest.ref} pseudonym of the new frame
#' of reference. If the \code{T.MAT} is \code{NULL}, then the returned object 
#' will contain only 4 matrices: "src.ref<-src.ref",
#' "src.ref<-dest.ref", "dest.ref<- dest.ref", "dest.ref<-src.ref".
#' @seealso \link[espadon]{ref.add}, \link[espadon]{ref.cutplane.add}, 
#' \link[espadon]{ref.remove}.
#' @examples
#' local.Tmat <-  ref.srctodest.add ("ref1","ref2", 
#'                                   TM = matrix(c (0.5, -sin (pi / 3), 0, 0,
#'                                                  sin (pi / 3), 0.5, 0, 0,
#'                                                  0, 0, 1, 0, 0, 0, 0, 1),
#'                                               ncol = 4))
#' str (local.Tmat)   
                                           

#' @export
#' @importFrom methods is
ref.srctodest.add <- function (src.ref, dest.ref, TM = diag (4), T.MAT = NULL){
  
 
  if (as.character(src.ref)==as.character(dest.ref) |
      as.character(src.ref)=="" | as.character(dest.ref)=="") stop ("Redefine src.ref or dest.ref.")
  
  
  new.list.matrix<- list ()
  new.list.matrix[[1]] <- diag (4)
  new.list.matrix[[2]] <- TM
  new.list.matrix[[3]] <- solve (TM) 
  new.list.matrix[[4]] <- diag (4)
  names (new.list.matrix) <- c (paste(src.ref, src.ref,sep="<-"),
                                paste(dest.ref,src.ref,sep="<-"),
                                paste(src.ref,dest.ref,sep="<-"),
                                paste(dest.ref,dest.ref,sep="<-"))
  new.matrix.description <- data.frame( t=names (new.list.matrix), src=c (src.ref, src.ref, dest.ref, dest.ref),
                                        dest=c (src.ref,dest.ref,src.ref,dest.ref), 
                                        type=rep("RIGID",4), stringsAsFactors = FALSE)
  
  
  db.file <- data.frame(t=paste(dest.ref,src.ref,sep="<-"),path = "local")
  ref.info <- data.frame( ref.pseudo=c(src.ref,dest.ref),ref=c("",""), stringsAsFactors = FALSE)
  if (is.null(T.MAT)){
    r <-(list (ref.info = ref.info,
               reg.info = list (patient = data.frame(patient=character(0),
                                                     patient.name=character(0),
                                                     patient.bd=character(0),
                                                     patient.sex=character(0)),
                                file = db.file),
               matrix.description = new.matrix.description, matrix.list = new.list.matrix))
    class (r) <- "t.mat"
    return(r)
  }
  
  if (!is (T.MAT, "t.mat")) stop ("T.MAT should be a t.mat class object.")
  if (!(src.ref %in% T.MAT$ref.info$ref.pseudo)) stop ("src.ref should be in T.MAT.")
  if (!(dest.ref %in% T.MAT$ref.info$ref.pseudo)) stop ("dest.ref should be in T.MAT.")  
  
  
  if (!is.null (T.MAT$matrix.list[[paste0 (dest.ref, "<-", src.ref)]])) 
    stop (paste ("The link between", src.ref, "and", dest.ref, "already exists."))
  
  
  ref.info <- T.MAT$ref.info
  ref.info_ <- ref.info[ ,1:2]
  reg.info <- T.MAT$reg.info
  reg.info$file <- rbind(reg.info$file, db.file)
  matrix.description <- unlist(lapply(ref.info_$ref.pseudo,
                                      function(r) 
                                        paste(paste(ref.info_$ref.pseudo,r,sep="<-"),r,ref.info_$ref.pseudo, sep=";")))
  matrix.description <- do.call(rbind.data.frame,strsplit(matrix.description,";"))
  names(matrix.description)  <- c("t","src","dest")
  matrix.description$type=""
  matrix.description [] <- data.frame (lapply(matrix.description [], as.character), stringsAsFactors=FALSE)
  
  list.matrix <-lapply(matrix.description$t,function (r) return(NULL))
  names(list.matrix) <- matrix.description$t
  for(i in 1:nrow(matrix.description)) {
    if (!is.null(T.MAT$matrix.list[[matrix.description$t[i]]])){
      list.matrix[[i]] <- T.MAT$matrix.list[[matrix.description$t[i]]]
      matrix.description[i,4] <-T.MAT$matrix.description[T.MAT$matrix.description$t==matrix.description$t[i],4]
    } else if (!is.null(new.list.matrix[[matrix.description$t[i]]])){
      list.matrix[[i]] <- new.list.matrix[[matrix.description$t[i]]]
      matrix.description[i,4] <-new.matrix.description[new.matrix.description$t==matrix.description$t[i],4]
    }
  }
  change <- TRUE
  while(change){
    idx.left <- which(matrix.description$type=="RIGID" & matrix.description$dest != matrix.description$src)
    change <- FALSE
    for (idx in idx.left){
      new.link <- which(matrix.description$src == matrix.description[idx,]$dest & matrix.description$type=="RIGID" &
                          matrix.description$dest != matrix.description[idx,]$src & matrix.description$dest != matrix.description$src )
      if (length(new.link)>0) for  (new.link.idx in 1:length(new.link)){
        to.complete.idx <- which(matrix.description$src == matrix.description[idx,]$src & matrix.description$dest==matrix.description[new.link[new.link.idx],]$dest)
        if (matrix.description[to.complete.idx, ]$type == ""){
          
          matrix.description[to.complete.idx, ]$type <- "RIGID"
          list.matrix[[matrix.description[to.complete.idx, ]$t]] <- list.matrix[[matrix.description[new.link[new.link.idx],]$t]] %*% list.matrix[[matrix.description[idx, ]$t]]
          change <- TRUE
        }
      }
    }
  }
  r <- (list (ref.info = ref.info, reg.info=reg.info, matrix.description = matrix.description, matrix.list = list.matrix))
  class (r) <- "t.mat"
  
  return(r)
}