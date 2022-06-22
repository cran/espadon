#' Converting DICOM files to .xlsx files
#' @description The \code{xlsx.from.dcm} function creates an Excel file from 
#' DICOM files.
#' @param dcm.filenames String vector, representing the list of full names of 
#' DICOM files.
#' @param xlsx.filenames String vector, representing the list of full names of 
#' created *.xlsx files. If \code{multipage = TRUE}, only the \code{xlsx.filenames[1]} 
#' is used. 
#' @param multipage Boolean. If TRUE, all \code{dcm.filenames} are converted into 
#' multiple pages of the same *.xlsx file.
#' @param txt.sep String. Used if \code{as.txt = TRUE}. Separator of the tag value elements.
#' @param txt.length Positive integer. Used if \code{as.txt = TRUE}. Maximum number 
#' of letters in the representation of the TAG value.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @return Returns a boolean vector, establishing the existence of the created 
#' Excel files.
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_dcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "PMrtplan", tmpdir = pat.dir, fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' list.files (pat.dir)
#' 

#' # Creating an Excel file
#' xlsx.fnames <- file.path (pat.dir, 
#'                           paste (basename (dcm.filename),"xlsx", sep = "."))
#' xlsx.from.dcm (dcm.filename, xlsx.fnames)
#' list.files (pat.dir)
#' 
#' # Cleaning temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @import openxlsx
#' @export
xlsx.from.dcm <- function(dcm.filenames, xlsx.filenames , multipage  = TRUE,
                          txt.sep = "\\", txt.length = 100, 
                          tag.dictionary = dicom.tag.dictionary ()) {
  
  if (!multipage & (length(dcm.filenames)!=length(xlsx.filenames))) {
    warning ("dcm.filenames and xlsx.filenames must have the same length.")
    return (rep (FALSE,length (dcm.filenames)))
  }
  flag <-file.exists(dcm.filenames) & !dir.exists(dcm.filenames) &
    !grepl("[.]doc$|[.]docx$|[.]dot$|[.]dotx$|[.]odt$|[.]xml$|[.]txt$",dcm.filenames) &
    !grepl("[.]pdf$|[.]ppt$|[.]pptx$|[.]tif$|[.]png$|[.]jpg$|[.]mp3$",dcm.filenames) &
    !grepl("[.]mp4$|[.]wmv$|[.]bmp$|[.]gif$|[.]svg$|[.]exe$|[.]bat$",dcm.filenames) &
    !grepl("[.]hex$|[.]csv$|[.]xls$|[.]xlsx$|[.]qs$|[.]Rdcm$",dcm.filenames)
  dcm.filenames <- dcm.filenames[flag]
  xlsx.filenames <- xlsx.filenames[flag]
  
  if (length(dcm.filenames)==0) return(NULL)
  if (multipage) wb <- createWorkbook ()
  page.idx <- 1
  for (f.idx in 1:length(dcm.filenames)){
    m <- dicom.raw.data.loader(dcm.filenames[f.idx])
    # dicom.df <- dicom.browser(m ,tag.dictionary = tag.dictionary)
    db <- dicom.parser(m, try.parse = TRUE,
                       txt.sep = txt.sep, txt.length = txt.length, 
                       tag.dictionary = tag.dictionary)
    
    if (!is.null(db)) {
      
      if ("(0008,0005)" %in% db$TAG) Encoding(db$Value) <- db$Value[db$TAG==("(0008,0005)")]#"UTF-8"
      if (!multipage) {
        wb <- createWorkbook ()
        openxlsx::addWorksheet (wb, as.character (page.idx))
        openxlsx::setColWidths (wb, as.character (page.idx), cols = 1, widths = 30)
        openxlsx::setColWidths (wb, sheet = as.character (page.idx), cols=2:ncol(db), widths = "auto")
        openxlsx::writeDataTable (wb, as.character(page.idx), x=db)
        
        if (!dir.exists(dirname(xlsx.filenames[f.idx]))) dir.create(dirname(xlsx.filenames[f.idx]),recursive = TRUE)
        openxlsx::saveWorkbook(wb,xlsx.filenames[f.idx], overwrite = TRUE)
      } else {
        openxlsx::addWorksheet (wb, as.character (page.idx))
        openxlsx::setColWidths (wb, as.character (page.idx), cols = 1, widths = 30)
        openxlsx::setColWidths (wb, sheet = as.character (page.idx), cols=2:ncol(db), widths = "auto")
        openxlsx::writeDataTable (wb, as.character(page.idx), x=db)
        page.idx <- page.idx + 1
      }
    
    }
  }
  if (multipage) {
    if (!dir.exists(dirname(xlsx.filenames[1]))) dir.create(dirname(xlsx.filenames[1]),recursive = TRUE)
    openxlsx::saveWorkbook(wb, xlsx.filenames[1], overwrite = TRUE) 
    return (file.exists(xlsx.filenames[1]))
  } 
  flag[flag] <- file.exists(xlsx.filenames)
  return (flag)
}