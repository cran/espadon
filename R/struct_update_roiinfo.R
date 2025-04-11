####################################################################
#' Update roi.info 
#' @description The \code{struct.update_roiinfo} function updates the roi.info 
#' element in the swordfish object of the struct class when point coordinates 
#' have been modified in the roi.data element.
#' @param struct espadon object of class struct.
#' @return Returns espadon object of class struct with updated roi.info element. 
#' @export
struct.update_roiinfo<- function(struct){
  if (is.null(struct$roi.data)) {
    struct$roi.info <- struct$roi.info [,1:6]
    return(struct)
  }
  struct$roi.info <- cbind(struct$roi.info [,1:6],
                       .struct.moreinfo (struct$roi.data, struct$ref.from.contour, struct$thickness))
  return(struct)
}