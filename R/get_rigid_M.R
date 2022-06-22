#' Transfer matrix between two frames of reference
#' @description The function \code{get.rigid.M} provides, from the T.MAT list created
#' by the functions \link[espadon]{load.patient.from.Rdcm},  
#' \link[espadon]{load.patient.from.dicom} or \link[espadon]{load.T.MAT}, the 4x4 
#' transfer matrix from the FoR (frame o reference) pseudonym \code{src.ref} to 
#' the FoR pseudonym \code{dest.ref}.
#' @param T.MAT "t.mat" class object, created by the functions 
#' \link[espadon]{load.patient.from.Rdcm}, \link[espadon]{load.patient.from.dicom} 
#'or \link[espadon]{load.T.MAT}
#' @param src.ref Pseudonym of the source frame of reference
#' @param dest.ref Pseudonym of the destination frame of reference
#' @return Returns the 4x4 transfer matrix \code{dest.ref} from \code{src.ref}.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = c ("ct", "mr"), roi.name = "",
#'                              dxyz = c(5, 5, 5))
#' get.rigid.M (patient$T.MAT, "ref1", "ref2")

#' @export
#' @importFrom methods is

get.rigid.M <- function (T.MAT, src.ref, dest.ref) {
  
  
  if (is.null (T.MAT)) {
    if (src.ref != dest.ref) return (NULL)
    return (matrix (c (rep ( c (1.0, 0.0, 0.0, 0.0, 0.0), 3), 1.0), nrow=4))
  }
  
  if (!is(T.MAT,"t.mat")) stop("T.MAT should be a t.mat class object.")
  T.name <- paste (dest.ref, src.ref, sep = "<-")
  T.info <- T.MAT$matrix.description[T.MAT$matrix.description == T.name, ]
  if (nrow(T.info)==1) if (castlow.str(T.info$type)== "rigid")
    return (T.MAT$matrix.list[[T.name]])
  return (NULL)
  
}