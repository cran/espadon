#' Deletion of a frame of reference in T.MAT
#' @description The \code{ref.remove} function removes the management of a frame
#' of reference in T.MAT.
#' @param ref.name Character string, pseudonym of the frame of reference to delete.
#' @param T.MAT "t.mat" class object in which the \code{ref.name} frame of 
#' reference is to be deleted.
#' @return Returns a "t.mat" class object, which no longer contains transfer 
#' matrices from or to the ref.pseudo \code{ref.name}.
#' \link[espadon]{ref.cutplane.add}.
#' @examples
#' # Adding of the reference frame "ref1_60", which is a 60 degree rotation of 
#' # reference frame "ref1".
#' orientation <- c (cos (pi / 3), sin (pi / 3), 0, 
#'                   -sin (pi / 3), cos (pi / 3), 0)
#'
#' local.Tmat <- ref.add (src.ref = "ref1", orientation = orientation,
#'                   new.ref.pseudo = "ref1_60")
#' str(local.Tmat)
#' 
#' # Removal of  "ref1_60"
#' local.Tmat <- ref.remove (ref.name =  "ref1_60", T.MAT = local.Tmat)
#' str(local.Tmat)

#' @export
#' @importFrom methods is

ref.remove <- function (ref.name, T.MAT) {
  
  if (!is (T.MAT, "t.mat")){
    warning ("T.MAT should be a t.mat class object.")
    return (NULL)
  }
  
  ref.info <- T.MAT$ref.info[T.MAT$ref.info$ref.pseudo!=ref.name,]
  row.names(ref.info) <- NULL
  
  
  reg.info <-T.MAT$reg.info 
  if (nrow(T.MAT$reg.info$file)!=0)
    reg.info$file <- T.MAT$reg.info$file[sapply(T.MAT$reg.info$file$t, 
                                                function(str) !(ref.name %in%unlist(strsplit(str,"<-")))),]
  
  matrix.description <- unlist(lapply(ref.info$ref.pseudo,function(r) paste(paste(ref.info$ref.pseudo,r,sep="<-"),r,ref.info$ref.pseudo, sep=";")))
  matrix.description <- do.call(rbind.data.frame,strsplit(matrix.description,";"))
  names(matrix.description)  <- c("t","src","dest")
  matrix.description$type=""
  matrix.description [] <- data.frame (lapply(matrix.description [], as.character), stringsAsFactors=FALSE)
  
  list.matrix <-lapply(matrix.description$t,function (r) return(NULL))
  names(list.matrix) <- matrix.description$t
  for(i in 1:nrow(matrix.description)) {
    if (!is.null(T.MAT$matrix.list[[matrix.description$t[i]]])){
      list.matrix[[i]] <- T.MAT$matrix.list[[matrix.description$t[i]]]
      matrix.description[i,]$type <-T.MAT$matrix.description[T.MAT$matrix.description$t==matrix.description$t[i], ]$type
    }
  }
  r <- list (ref.info = ref.info, reg.info = reg.info, matrix.description = matrix.description, matrix.list = list.matrix)
  class(r) <- "t.mat"
  return(r)
}