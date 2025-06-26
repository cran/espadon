#' Updating Rdcm files.
#' @description The \code{Rdcm.upgrade} function updates Rdcm files that were 
#' created with a previous version.
#' @param Rdcm.files String vector, representing the list of the full names of the  
#' Rdcm files, or its directories.
#' @return Saves the updated Rdcm files. If the Rdcm files were generated from 
#' the dicom files, the data is updated from the DICOM fields.
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
#' # test of Rdcm.upgrade
#' 
#' Rdcm.upgrade (pat.dir)
#' # or
#' Rdcm.upgrade (list.files (pat.dir, full.names = TRUE))
#' 
#' # Cleaning  temporary directories
#' unlink (pat.dir, recursive = TRUE)
#' 
#' @export
#'

#' @import progress
Rdcm.upgrade <- function(Rdcm.files){
  if (length(Rdcm.files)==0) {
    warning ("no files to load.")
    return (NULL)
  } 
  flag <- dir.exists(Rdcm.files)
  Rdcm.dir <- Rdcm.files[flag]
  Rdcm.filenames1 <-  list.files(Rdcm.dir,recursive = TRUE,full.names = TRUE)
  Rdcm.filenames2 <- Rdcm.files[!flag]
  Rdcm.filenames2 <- Rdcm.filenames2[file.exists(Rdcm.filenames2)]
  Rdcm.filenames <- c(Rdcm.filenames1,Rdcm.filenames2)
  
  if (length(Rdcm.filenames)==0) {
    warning ("no files to load.")
    return (NULL)
  } 
  pb <- progress_bar$new(format = " objects creation [:bar] :percent",
                         total = length(Rdcm.filenames),  width= 60)
  for(fn in Rdcm.filenames) {
    obj <-suppressWarnings(tryCatch(load.Rdcm.raw.data(fn, address= TRUE, data=TRUE,
                             upgrade.to.latest.version = TRUE),
                   error = function (e) list(stopmessage= e$message)))
    if  (!is.null(obj$stopmessage)){
      warning(paste(obj$stopmessage, basename(fn)," can not be upgraded."))
    } else {
      if (obj$update.needed) .save.dicom.raw.data.to.Rdcm(obj, fn)
    }
    pb$tick ()
  }
}