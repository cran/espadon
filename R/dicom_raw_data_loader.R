#' DICOM file loading in raw data
#' @description the \code{dicom.raw.data.loader} function loads a DICOM file as 
#' raw data.
#' @param dcm.filename Character string, representing the full name of a DICOM file.
#'
#' @return Returns a vector of raw data from \code{dcm.filename}.
#' @seealso \link[espadon]{dicom.browser}, \link[espadon]{dicom.tag.parser} 
#' @export
#'
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file for testing.
#' pat.src.dir <- file.path (tempdir(), "toy_dccm")
#' dir.create (pat.src.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = pat.src.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' 
#' # loading of file
#' dicom.raw.data <- dicom.raw.data.loader (dcm.filename)
#' 
#' # checks if it is consistent with the original raw data
#' all ( dicom.raw.data == toy.dicom.raw () )
#' 
#' # Cleaning  temporary directory
#' unlink (pat.src.dir, recursive = TRUE)

dicom.raw.data.loader <- function (dcm.filename){
  dcm.filename <- dcm.filename[1]
  return (readBin (dcm.filename, raw(), n=file.size (dcm.filename), size = 1, signed=FALSE))
}
