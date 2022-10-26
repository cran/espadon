############################################################################
#' Loading a *.Rdcm file
#' @description the \code{load.Rdcm.raw.data} function loads the content of a *.Rdcm file.
#' @param Rdcm.filename Character string, representing the full name of a *.Rdcm 
#' file created by \link[espadon]{dicom.to.Rdcm.converter}.
#' @param address boolean. If TRUE, a dataframe with the address of the tags in 
#' the raw DICOM data is returned.
#' @param data boolean. If TRUE, the DICOM information are returned as an R list.
#' @return Returns a list containing the information, converted by \pkg{espadon}, of a 
#' DICOM object..
#' @seealso \link[espadon]{dicom.to.Rdcm.converter}, \link[espadon]{load.obj.from.Rdcm}. 
#' @import qs
#' @examples
#' # For testing, save first toy.dicom.raw () raw data to a temporary file, and
#' # convert it in Rdcm fie
#' pat.src.dir <- file.path (tempdir(), "PM_dcm") 
#' dir.create (pat.src.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "PM_rtplan", tmpdir = pat.src.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' pat.dir <- file.path (tempdir(), "PM_Rdcm")
#' dicom.to.Rdcm.converter (pat.src.dir, pat.dir, update = TRUE)
#' lf <- list.files (pat.dir, pattern = "[.]Rdcm$", full.names = TRUE)
#' lf
#' 
#' # Inspect Rdcm raw data
#' L <- load.Rdcm.raw.data (lf[1])
#' str (L, max.level =3)

#' @export
load.Rdcm.raw.data <- function (Rdcm.filename, address= TRUE, data=TRUE) {
  if (!file.exists(Rdcm.filename)) return (NULL)
  zz <-  file(Rdcm.filename, "rb")
  l <- readBin(zz,what="int",size=4, n=3, endian="little")
  h <- qdeserialize (readBin(zz,what="raw", n=l[1]))
  a <- NULL
  if (is.null(h$espadon.version)){
    espadon.version <- character(0)
  } else {
    espadon.version <- h$espadon.version
    h$espadon.version <- NULL
  }
  if ((address | data) & l[2]>0) a <-  qdeserialize (readBin(zz,what="raw", n=l[2]))
  if (data) d <-  qdeserialize (readBin(zz,what="raw", n=l[3]))
  close (zz)
  h$file.dirname <- dirname (Rdcm.filename)
  h$file.basename <- basename (Rdcm.filename)
  from.dcm <- l[2]>0
  
  if (length(espadon.version)==0) {
    n <- names(h)
    idx <- which(n=="ref.object.name")
    if (length(idx)!=0) n[idx] <- "ref.object.alias"
    idx <- which(n=="patient.xyz0")  
    if (length(idx)!=0) n[idx] <- "xyz0"  
    idx <- which(n=="patient.orientation")  
    if (length(idx)!=0) n[idx] <- "orientation"  

    names(h) <- n
  }
  
  if (!address & !data) return(list(header=h, from.dcm=from.dcm))
  if (address & !data) return(list(header=h, address=a, from.dcm=from.dcm))
  if (!address & data) return(list(header=h, data=d, from.dcm=from.dcm))
  if (address & data) return(list(header=h, address=a, data=d, from.dcm=from.dcm,
                                  espadon.version=espadon.version))
}