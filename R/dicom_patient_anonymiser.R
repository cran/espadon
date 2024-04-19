#' Anonymisation of a patient's DICOM files
#' @description the \code{dicom.patient.anonymiser} function anonymises all DICOM 
#' files in a patient's directory.
#' @param dcm.files String vector, representing the list of the full names of the  
#' DICOM files of the same patient, or its directories.
#' @param pat.dest.dir Character string,representing the full name of the patient's 
#' directory, which will contain the patient's anonymized files.
#' @param offset Integer, default to 0. Each date of the DICOM will be shifted 
#' by this offset expressed in days.
#' @param new.PIN Character string, representing the PIN remplacing the old one.
#' @param reset.private.tag Boolean, if \code{TRUE}, the value of tags that are 
#' not in the \code{tag.dictionary} is removed.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, whose structure it must keep. This 
#' dataframe is used to parse DICOM files.
#' @param verbose Boolean. If \code{TRUE}, a progress bar indicates loading progress.

#' @return Creation of the \code{pat.dest.dir} directory, with anonymous DICOM files
#' @note The files are anonymized as follows:  
#' \itemize{
#' \item Each date of the DICOM file will be shifted by \code{offset} expressed in days.
#' \item Each patient's name, and patient'ID are remplaced by \code{new.PIN}
#' \item All other patient data are deleted, except age, weight, height, gender 
#' and shifted birthday.
#' \item All address, phone, physician, operator, author, reviewer, service.
#' \item If \code{reset.private.tag = TRUE}, the values of the tags not contained 
#' in the \code{tag.dictionary} are deleted.
#' }
#' File names are composed of their modality and the SOP UID.
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file pat.dir for testing.
#' temp <- tempdir()
#' pat.dir <- file.path (temp, "toy_dcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = pat.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' 
#' # Files anonymisation
#' anonymous.pat.dir <- file.path (temp, "anonymous") 
#' dicom.patient.anonymiser (dcm.files = pat.dir, pat.dest.dir = anonymous.pat.dir,
#'                           offset = 0, new.PIN = "Anonymous patient",
#'                           reset.private.tag = TRUE)
#' lf <- list.files(anonymous.pat.dir, full.names = TRUE)  
#' dp <- dicom.parser(lf[1])   
#' dp[grep("^[(]0008|^[(]0010", dp$TAG),]                   
#' @export
#' @import progress
dicom.patient.anonymiser <- function(dcm.files, pat.dest.dir, offset = 0,
                                     new.PIN = "Anonymous patient",
                                     reset.private.tag = FALSE,
                                     tag.dictionary = dicom.tag.dictionary(),
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
  
  if (!dir.exists(pat.dest.dir)) dir.create(pat.dest.dir, recursive =TRUE)
  if (verbose) pb <- progress_bar$new(format = " downloading [:bar] :percent",
                                      total = length(dcm.filenames), width= 60)
  
  for(fn in dcm.filenames) {
    raw.data <- dicom.raw.data.loader(fn)
    dicom.df <- dicom.browser(raw.data, stop.tag = "(0008,0060)", full.info = TRUE)
    m <- match(c("(0008,0060)","(0008,0018)"),dicom.df$tag)
    new.fn <- paste0(paste(sapply(m, function(idx) {
      dicom.tag.parser (dicom.df$start[idx], dicom.df$stop[idx], 
                        dicom.df$VR[idx], dicom.df$endian[idx],
                        raw.data, try.parse= FALSE)}), collapse ="_"),".dcm")
    
    an.raw.data <- dicom.raw.data.anonymizer(
      raw.data,
      offset = offset,
      new.PIN = new.PIN,
      reset.private.tag = reset.private.tag,
      tag.dictionary = tag.dictionary)
    zz <- file (file.path(pat.dest.dir, new.fn), "wb")
    writeBin (an.raw.data, zz, size = 1)
    close(zz)
    if (verbose) pb$tick()
  }
  
}