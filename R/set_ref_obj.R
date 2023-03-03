
#' Set the reference objects of a espadon object
#' @description The function \code{set.reference.obj} adds to an espadon object 
#' the information identifying the espadon objects from which it derives.
#' @param obj espadon object of class "dvh", "fan", "histo", "histo2D", "mesh", 
#' "rtplan", "struct", "undef" or "volume".
#' @param ref.obj espadon object of class "dvh", "fan", "histo", "histo2D", "mesh", 
#' "rtplan", "struct", "undef" or "volume". List of espadon objects.
#' @param add Boolean. If TRUE, the reference objects are added to those already 
#' contained by \code{obj}.
#' @return Returns the espadon object \code{obj}, containing the ref.object.alias 
#' and ref.object.info fields identifying its reference objects
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 6
#' pat <- toy.load.patient (modality = c("ct", "rtdose", "rtstruct"),
#'                          roi.name = c("eye"), dxyz = rep (step, 3),
#'                          beam.nb = 3)
#' display.obj.links(pat)                          
#' pat$rtstruct[[1]] <- set.reference.obj(pat$rtstruct[[1]],pat$ct[[1]])  
#' display.obj.links(pat)                          
#' @export
set.reference.obj <- function (obj, ref.obj, add =TRUE){
  if (is.na(match(class(obj),espadon.class()[-c(6,7,10)]))) stop ("obj must be a espadon object")
  if (!is.na(match(class(ref.obj),espadon.class()[-c(6,7,10)]))) ref.obj <-list(ref.obj)
  class.ref <- sapply(ref.obj,class)
  m <- match(class.ref,espadon.class()[-c(6,7,10)])
  ref.obj <- ref.obj[!is.na(m)]
  return(.set.ref.obj(obj, ref.obj, add=add))
}


.set.ref.obj<- function(obj,ref,add=TRUE){
  
  if (!is.na(match(class(ref),espadon.class()[-c(6,7,10)]))) ref <-list(ref)
  
  if (!add){
    obj$ref.object.alias <- character(0)
    obj$ref.object.info <- list(SOP.ID = character(0), SOP.label = character(0))
  } else {
    obj$ref.object.alias  <- c(obj$ref.object.alias ,character(0))
    obj$ref.object.info <- list(SOP.ID = c(obj$ref.object.info$SOP.ID,character(0)),
                                SOP.label =  c(obj$ref.object.info$SOP.label,character(0)))
  }


  for(ref.obj in ref){
    obj$ref.object.alias <- unique(c(obj$ref.object.alias , ref.obj$object.alias))
    if (!is.null(ref.obj$object.info)){
      obj$ref.object.info$SOP.ID <- unique(c(obj$ref.object.info$SOP.ID,ref.obj$object.info$SOP.ID))
      obj$ref.object.info$SOP.label <- unique(c(obj$ref.object.info$SOP.label,ref.obj$object.info$SOP.label))
    } 
  }
  if(length(obj$ref.object.info$SOP.ID)==0 & length(obj$ref.object.info$SOP.label)==0) 
    obj$ref.object.info <- NULL
  
  idx <- match(obj$ref.object.alias,"")
  obj$ref.object.alias <- obj$ref.object.alias[is.na(idx)] 
  
  if (length(obj$ref.object.alias)==0) obj$ref.object.alias <- NULL
  
  cl <- class(obj)
  n <- names(obj)
  idx <- match(c("object.alias","ref.object.alias"), n)
  if (!(is.na(idx[2]))){
    obj <- obj[c(1:idx[1], idx[2],(1:length(n))[-c(1:idx[1],idx[2])])] 
    n <- names(obj)
    idx <- match(c("ref.object.alias","object.info", "ref.object.info"), n)
    idx1 <-max(idx[1:2], na.rm =TRUE)
    if (!(is.na(idx[3])))
      obj <- obj[c(1:idx1, idx[3],(1:length(n))[-c(1:idx1,idx[3])])] 
    class(obj) <- cl
  }
  return(obj)
}

.resort.obj<-function(obj){
  
  if (length(obj$ref.object.alias)>1) {
    idx <- match(obj$ref.object.alias,"")
    obj$ref.object.alias <- obj$ref.object.alias[is.na(idx)] 
  }
  
  cl <- class(obj)
  n <- names(obj)
  idx <- match(c("object.alias","ref.object.alias"), n)
  if (!(is.na(idx[2]))){
    obj <- obj[c(1:idx[1], idx[2],(1:length(n))[-c(1:idx[1],idx[2])])] 
    n <- names(obj)
    idx <- match(c("ref.object.alias","object.info", "ref.object.info"), n)
    idx1 <-max(idx[1:2], na.rm =TRUE)
    if (!(is.na(idx[3])))
      obj <- obj[c(1:idx1, idx[3],(1:length(n))[-c(1:idx1,idx[3])])] 
    class(obj) <- cl
  }
  return(obj)
}