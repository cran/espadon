#' Conversion of DICOM object into files that can be interpreted by the \pkg{espadon} 
#' package
#' @description The \code{dicom.to.Rdcm.converter} function creates, for each DICOM
#' object, a *.Rdcm file usefull for using \pkg{espadon} package. Each Rdcm file created is
#' referenced by the date of acquisition of the object (if it is not available, 
#' its creation date), the patient's PIN, a reference number, an object number 
#' in this reference system, and the object modality (mr, ct, rtstruct...).
#' @param dcm.files String vector, representing the list of the full names of 
#' the DICOM files of the same patient, or its directory.
#' @param pat.dest.dir Character string, representing the full name of patient 
#' directory, which will contain files converted \pkg{espadon}. 
#' @param update Boolean. If set to \code{TRUE}, and if \code{pat.dest.dir} 
#' contains previously converted files, these files are updated,even if they are 
#' duplicated. They retain the same \pkg{espadon} reference frame assignment.
#' @param ignore.duplicates Boolean. If \code{TRUE}, the function ignores duplicated objects.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @param verbose Boolean. If \code{TRUE}, a progress bar indicates the progress 
#' of the conversion.
#' @note For each DICOM object, \code{dicom.to.Rdcm.converter} 
#' creates a *.Rdcm file whose basename is made up of the date of the acquisition 
#' (or creation date if previous not found), the patient's PIN, the pseudonym of 
#' the frame of reference ("ref1", "ref2"...), the number of the volume object 
#' in the directory in this frame of reference ("do1", "do2"...), and the object 
#' modality ("mr", "ct", "rtdose", "rtstruct"...). 
#' @note For example: \code{BASE = "20160514_a008e9ac_ref2_do1_mr"}
#' 
#' @return Returns the list of basenames of the created files. 
#' @return Returns \code{NULL} if there are no DICOM files in \code{dcm.files}
#' @export
#'
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file for testing.
#' pat.src.dir <- file.path (tempdir(), "PM_dcm") 
#' dir.create (pat.src.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "PM_rtplan", tmpdir = pat.src.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' 
#' # Create a temporary destination directory where the *.Rdcm file will be saved
#' pat.dest.dir <- file.path (tempdir(), "PM_Rdcm")
#' 
#' dicom.to.Rdcm.converter (pat.src.dir, pat.dest.dir, update = TRUE)
#' # or
#' dicom.to.Rdcm.converter (dcm.filename, pat.dest.dir, update = TRUE)
#' 
#' list.files (pat.dest.dir)
#' 
#' # Cleaning  temporary directories
#' unlink (pat.src.dir, recursive = TRUE)
#' unlink (pat.dest.dir, recursive = TRUE)


dicom.to.Rdcm.converter <- function (dcm.files, pat.dest.dir, update = TRUE, 
                                     ignore.duplicates = FALSE,
                                     tag.dictionary = dicom.tag.dictionary (),
                                     verbose = TRUE){
  
  if (length(dcm.files)==0) {
    warning ("no files to format.")
    return (NULL)
  } 
  flag <- dir.exists(dcm.files)
  dcm.dir <- dcm.files[flag]
  dcm.filenames1 <-  list.files(dcm.dir,recursive = TRUE,full.names = TRUE)
  dcm.filenames2 <- dcm.files[!flag]
  dcm.filenames2 <- dcm.filenames2[file.exists(dcm.filenames2)]
  dcm.filenames <- c(dcm.filenames1,dcm.filenames2)
  if (length(dcm.filenames)==0) {
    warning ("no files to format.")
    return (NULL)
  } 
  
  pat.dest.dir <- pat.dest.dir[1]
  pat.dest.dir <- gsub("\\\\|/", .Platform$file.sep, pat.dest.dir)
  flag <- file.exists(dcm.filenames) & !dir.exists(dcm.filenames) &
    !grepl("[.]doc$|[.]docx$|[.]dot$|[.]dotx$|[.]odt$|[.]xml$|[.]txt$",dcm.filenames) &
    !grepl("[.]pdf$|[.]ppt$|[.]pptx$|[.]tif$|[.]png$|[.]jpg$|[.]mp3$",dcm.filenames) &
    !grepl("[.]mp4$|[.]wmv$|[.]bmp$|[.]gif$|[.]svg$|[.]exe$|[.]bat$",dcm.filenames) &
    !grepl("[.]hex$|[.]csv$|[.]xls$|[.]xlsx$|[.]qs$|[.]Rdcm$",dcm.filenames)
  dcm.filenames <- dcm.filenames[flag]
  
  if (length (dcm.filenames)==0) {
    warning ("no files to format.")
    return (NULL)
    } 
   
      
  if (!dir.exists(pat.dest.dir)) dir.create(pat.dest.dir,recursive = TRUE)
  
  if (!dir.exists(pat.dest.dir)) {
    warning (paste0("can not create ",pat.dest.dir,"."))
    return (NULL)
    } 
        
        
  L.list <- .create.Rdcm.object (dcm.filenames = dcm.filenames, 
                                 existing.obj = Rdcm.inventory (pat.dest.dir), 
                                 verbose = verbose, update = update, 
                                 save.flag = TRUE, save.dir = pat.dest.dir, 
                                 only.header = FALSE, ignore.duplicates =  ignore.duplicates,
                                 tag.dictionary = tag.dictionary)
        
      
  
    
  
  return (L.list)
}