#' Loading an \pkg{espadon} object from *.Rdcm file
#' @description The \code{load.obj.from.Rdcm} function loads a DICOM object into 
#' memory, creating a list containing the information necessary for its subsequent 
#' use with the \pkg{espadon} package.
#' @param Rdcm.filename Character string, representing the full name of a *.Rdcm 
#' file created by \link[espadon]{dicom.to.Rdcm.converter}.
#' @param data Boolean. Only works for objects usable by the \pkg{espadon} package, namely
#' ct, mr, rtdose, rtstruct, pt... If \code{data = TRUE}, either the values of the 
#' voxels when modality is (ct, mr, rtdose), or the coordinates of the RoI when modality 
#' is rtstruct, are loaded into memory.
#' @param nb Vector of integers, active only if \code{data = TRUE}, and only
#' operating on rtstruct. If \code{nb = NULL}, all the RoI of rtstruct are
#' loaded into memory. Otherwise only data of the RoI indices defined by the vector 
#' \code{nb} are loaded.
#' @param upgrade.to.latest.version Boolean. If \code{TRUE}, the function attempts 
#' to upgrade to the latest version, parsing the DICOM data. It may take longer 
#' to load the data. Consider using the \link[espadon]{Rdcm.upgrade} function.
#' @return Returns an \pkg{espadon} object of class "dvh","histo","histo2D","mesh", 
#' "rtplan","struct", "undef" or "volume" depending on the object modality. See 
#' \link[espadon]{espadon.class} for class definitions.

#' @seealso \link[espadon]{load.obj.data} and \link[espadon]{load.obj.from.dicom}
#' @examples
#' # First, save toy patient objects to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = c (4, 4, 4))
#' save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
#' save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
#' save.T.MAT (patient$T.MAT, dirname = pat.dir)
#' # Rdcm files in pat.dir
#' list.files(pat.dir)
#' 

#' CT <- load.obj.from.Rdcm (file.path (pat.dir, 
#'                                      list.files(pat.dir, pattern="ct1[.]Rdcm")[1]),
#'                           data=TRUE)
#' MR <- load.obj.from.Rdcm (file.path (pat.dir, 
#'                                      list.files(pat.dir, pattern="mr1[.]Rdcm")[1]),
#'                           data=TRUE)
#' Reg <-load.obj.from.Rdcm (file.path (pat.dir,"ref1_from_ref2.Rdcm"), data=TRUE)     
#' str(Reg)
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)  
                  
#' @export
load.obj.from.Rdcm <- function (Rdcm.filename, data = TRUE, nb = NULL,
                                upgrade.to.latest.version = FALSE){
  if (length(Rdcm.filename)>1){
    message  ("only first file is used.")
    Rdcm.filename <- Rdcm.filename[1]
  }
  if ((length(Rdcm.filename)==0) | !file.exists(Rdcm.filename)){
    warning ("no file to load.")
    return (NULL)
  }
  
  if (grepl ("[.]Rdcm$", Rdcm.filename)){
    Lobj <- load.Rdcm.raw.data (Rdcm.filename,  data=data, 
                                upgrade.to.latest.version = upgrade.to.latest.version)
    if (Lobj$update.needed){

      if (upgrade.to.latest.version) {
        warning(paste(Lobj$header$file.basename, "version has been upgraded, consider using Rdcm.upgrade() for faster loading."))
      } else {
        warning(paste(Lobj$header$file.basename, "version is not up to date, consider using Rdcm.upgrade()."))  
      }
    }
    return (.load.object (Lobj, data=data, nb=nb))
  }
  return (NULL)
}