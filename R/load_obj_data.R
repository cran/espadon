#' Load data of an \pkg{espadon} class object
#' @description The \code{load.obj.data} function loads all the data of an \pkg{espadon} 
#' object of class '\code{struct}' or '\code{volume}'.
#' @param obj  \code{struct} or "volume" class object
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, whose structure it must keep. 
#' This dataframe is used to parse DICOM files in case obj was extracted from DICOM files.

#' @return Returns the \pkg{espadon} object with data \code{$vol3D.data} or {$roi.data}
#' @seealso \link[espadon]{load.obj.from.dicom} and \link[espadon]{load.obj.from.Rdcm}

#' @examples
#' # First, save toy patient objects to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = c (4, 4, 4))
#' save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
#' rm( patient)
#' 
#' 
#' patient <- load.patient.from.Rdcm (pat.dir, data = FALSE)
#' CT <- load.obj.data (patient$ct[[1]])
#' str (CT, max.level = 2)

#' @export
#' @importFrom methods is
load.obj.data <- function (obj, tag.dictionary = dicom.tag.dictionary ()){

  if (!(is (obj, "volume") | is (obj, "struct"))) stop ("obj must be of class 'struct' or 'volume'")
  fp <- file.path (obj$file.dirname, obj$file.basename)
  ref.pseudo <- obj$ref.pseudo
  if (grepl("[.]Rdcm$",fp[1])) return(load.obj.from.Rdcm(fp[1]))
  obj <- load.obj.from.dicom(fp,ref.pseudo = ref.pseudo, 
                             tag.dictionary=tag.dictionary,
                             verbose=F)
  return (obj)  
}