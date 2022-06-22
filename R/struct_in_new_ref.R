#' Change of frame of reference of a "struct" class object.
#' @description The \code{struct.in.new.ref} function allows you to change the 
#' frame of reference of a struct.
#' @param struct "struct" class object.
#' @param new.ref.pseudo pseudonym of the frame of reference in which the struct 
#' should be located. This \code{new.ref.pseudo} must exist in the \code{T.MAT} list.
#' @param T.MAT "t.mat" class object, created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
#' \link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.
#' @param alias Character string, \code{$alias} of the created object.
#' @return Returns "struct" class object in the new frame of reference 
#' \code{new.ref.pseudo}.
#' @seealso \link[espadon]{vol.in.new.ref}
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("mr", "rtstruct"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#' S <- patient$rtstruct[[1]]
#' S.in.new.ref <- struct.in.new.ref (S, patient$mr[[1]]$ref.pseudo, patient$T.MAT)

#' @export
#' @importFrom methods is

struct.in.new.ref <- function (struct, new.ref.pseudo, T.MAT, alias="") {
  
  if (!is (struct, "struct")) {
    warning ("struct should be a struct class object.")
    return (NULL)
  }
  if (new.ref.pseudo!=struct$ref.pseudo){
    if (!is (T.MAT, "t.mat")) {
      warning ("T.MAT should be a t.mat class object.")
      return (NULL)
    }}
  
  if (is.null(struct$roi.data))  {
    warning ("struct$roi.data is NULL.")
    
  }
  M <- get.rigid.M (T.MAT, struct$ref.pseudo, new.ref.pseudo)
  if (is.null (M)) return (NULL)
  
  struct$object.alias <- alias
  struct$ref.pseudo <- new.ref.pseudo
  struct$frame.of.reference <- T.MAT$ref.info[T.MAT$ref.info$ref.pseudo==new.ref.pseudo, ]$ref
  
  struct$ref.from.contour <- M %*% struct$ref.from.contour
  struct$ref.from.contour [abs (struct$ref.from.contour) < 1.0e-6] <- 0
  if (!is.null(struct$roi.data))
    struct$roi.info[ , 7:17] <- .struct.moreinfo(struct$roi.data,
                                                 struct$ref.from.contour,
                                                 struct$thickness)
  return (struct)
}