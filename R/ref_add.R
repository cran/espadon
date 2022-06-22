#' Adding a frame of reference in T.MAT
#' @description The \code{ref.add} function adds the transfer matrices from or 
#' to a new frame of reference defined from 2 unit vectors and an origin point.
#' @param src.ref Character string, pseudonym of the frame of reference in which 
#' the \code{orientation} vector and the origin point \code{origin} are defined.
#' @param orientation Vector of 6 or 9 elements, composed of the coordinates of the 2 orthonormal vectors (i, j),
#' or of the 3 orthonormal vectors (i, j, k) of the new coordinate system, 
#' in the \code{src.ref} frame of reference.
#' @param origin Vector of the x, y, z coordinates of the origin point of the 
#' new frame of reference in the \code{src.ref} frame of reference. Default to 
#' c (0, 0, 0).
#' @param new.ref.pseudo Character string, pseudonym of the new frame of 
#' reference to add.
#' @param T.MAT "t.mat" class object created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, then only the link 
#' between \code{src.ref} and \code{new.ref.pseudo} is computed.

#' @return Returns a "t.mat" class object, which contains the transfer 
#' matrices from or to \code{new.ref.pseudo} pseudonym of the new frame
#' of reference. If the \code{T.MAT} is \code{NULL}, then the returned object 
#' will contain only 4 matrices: "src.ref<-src.ref", "src.ref<-new.ref.pseudo", 
#' "new.ref.pseudo<- new.ref.pseudo", "new.ref.pseudo<-src.ref".
#' @return Returns a \code{NULL} if  \code{orientation} is not well defined.
#' 
#' @seealso \link[espadon]{ref.cutplane.add}, \link[espadon]{ref.remove}.
#' @examples
#' # Adding of the reference frame "ref1_60", which is a 60 degree rotation of 
#' # reference frame "ref1".
#' orientation <- c (cos (pi / 3), sin (pi / 3), 0, 
#'                   -sin (pi / 3), cos (pi / 3), 0)
#'
#' local.Tmat <- ref.add (src.ref = "ref1", orientation = orientation,
#'                   new.ref.pseudo = "ref1_60")
#'
#' str(local.Tmat)
#' @export
#' @importFrom methods is
ref.add <- function (src.ref, orientation = c (1, 0, 0, 0, 1, 0),
                     origin = c (0, 0, 0), 
                     new.ref.pseudo = "newref",
                     T.MAT = NULL){
  
  if (length (orientation)!=6 & length (orientation)!=9) {
    warning ("bad orientation length.")
    return (NULL)
  }
  
  if (length (orientation)>=6 & (round(as.numeric(orientation[1:3]%*%orientation[1:3])-1,9)> 1e-6 |
                                 round(as.numeric(orientation[4:6]%*%orientation[4:6])-1,9)> 1e-6 |
                                 round(as.numeric(orientation[1:3]%*%orientation[4:6]),9) > 1e-6)) {
    warning ("orientation is not an orthonormal basis.")
    return (NULL)
  }
  
  if (length (orientation)==9 & (round(as.numeric(orientation[7:9]%*%orientation[7:9])-1, 9) > 1e-6 |
                                 round(as.numeric(orientation[4:6]%*%orientation[7:9]), 9) > 1e-6 |
                                 round(as.numeric(orientation[1:3]%*%orientation[4:6]), 9) > 1e-6)) {
    warning ("plane.orientation is not an orthonormal basis.")
    return (NULL)
  }
  
  
  TM <- .ref.create (orientation, origin)
  
  new.list.matrix<- list ()
  new.list.matrix[[1]] <- diag (4)
  new.list.matrix[[2]] <- TM
  new.list.matrix[[3]] <- solve (TM)
  new.list.matrix[[4]] <- diag (4)
  names (new.list.matrix) <- c (paste(src.ref, src.ref,sep="<-"),
                                paste(new.ref.pseudo,src.ref,sep="<-"),
                                paste(src.ref,new.ref.pseudo,sep="<-"),
                                paste(new.ref.pseudo,new.ref.pseudo,sep="<-"))
  new.matrix.description <- data.frame( t=names (new.list.matrix), src=c (src.ref, src.ref, new.ref.pseudo, new.ref.pseudo),
                                        dest=c (src.ref,new.ref.pseudo,src.ref,new.ref.pseudo), 
                                        type=rep("RIGID",4), stringsAsFactors = FALSE)

  db.file <- data.frame(t=paste(new.ref.pseudo,src.ref,sep="<-"),path = "local")
  if (is.null(T.MAT)) {
    ref.info <- data.frame( ref.pseudo=c(src.ref,new.ref.pseudo),ref=c("",""), stringsAsFactors = FALSE)
    r <-(list (ref.info = ref.info,
               reg.info = list (patient = data.frame(patient=character(0),
                                                     patient.bd=character(0),
                                                     patient.sex=character(0)),
                                file = db.file),
               matrix.description = new.matrix.description, matrix.list = new.list.matrix))
    class (r) <- "t.mat"
    return(r)
  }
  
  if (!is (T.MAT, "t.mat")) {
    stop ("T.MAT should be a t.mat class object.")
  }
  
  
  if (new.ref.pseudo %in% T.MAT$ref.info$ref.pseudo) {
    stop ("added ref is already in T.MAT. Remove it first with ref.remove")

  }
  
  ref <-paste (sample(1:9, 20, replace=TRUE),collapse ="")
  
  while (ref %in% T.MAT$ref.info$ref) ref <-paste (sample(1:9, 20, replace=TRUE),collapse ="")
  
  ref.info <- data.frame( ref.pseudo=c (T.MAT$ref.info$ref.pseudo,new.ref.pseudo),ref=c(T.MAT$ref.info$ref,ref), stringsAsFactors = FALSE)
  reg.info <- T.MAT$reg.info
  reg.info$file <- rbind(reg.info$file, db.file)
  matrix.description <- unlist(lapply(ref.info$ref.pseudo,
                                      function(r) 
                                        paste(paste(ref.info$ref.pseudo,r,sep="<-"),r,ref.info$ref.pseudo, sep=";")))
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