#' DICOM raw data browser
#' @description the \code{dicom.browser} function creates a dataframe describing 
#' the tags contained in the raw data of a DICOM file, as well as the information 
#' to access them.
#' @param dicom.raw.data Raw vector, representing the binary extraction of the 
#' DICOM file.
#' @param nbTAG Integer. If \code{nbTAG = 0} (default), and \code{stop.tag = ""}, 
#' all the DICOM raw data is browsed. Otherwise, the function only browses the 
#' first \code{nbTAG} tags.
#' @param stop.tag Character string, representing the tag that stops the browse 
#' of the \code{dicom.raw.data}.
#' @param stop.level Positive integer, specifying the encapsulation level of the 
#' \code{stop.tag} in \code{dicom.raw.data}.
#' @param full.info Boolean. If \code{TRUE}, more information about the DICOM 
#' data is returned.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @return Returns a dataframe if \code{dicom.raw.data} is DICOM raw data, 
#' \code{NULL} otherwise.
#' @return If \code{full.info = FALSE}, dataframe columns are 
#' \itemize{
#' \item tag : the tags contained in \code{dicom.raw.data},
#' \item VR : value representation of the content of the tag,
#' \item endian : the endianness of the tag content,
#' \item start : the start address in \code{dicom.raw.data} of the tag content.
#' \item stop : the stop address  in \code{dicom.raw.data} of the tag content.
#' }
#' @return If \code{full.info = TRUE}, the following columns are added :
#' \itemize{
#' \item encaps.load : If the tag contains encapsulated data, this column gives 
#' the number of bytes remaining until the end of the encapsulation. If there 
#' are several levels of encapsulation, these numbers are collapsed and separated 
#' by a space.
#' \item load.start : the start address in \code{dicom.raw.data} of the tag load size.
#' \item load.stop : the stop address in \code{dicom.raw.data}of the tag load size.
#' \item tag.start : the start address in \code{dicom.raw.data} of the tag.
#' } 
#' @seealso \link[espadon]{dicom.raw.data.loader}, \link[espadon]{dicom.tag.parser} 
#' @export
#' @examples
 
#' # DICOM information dataframe of the dummy raw data toy.dicom.raw ()
#' df <- dicom.browser (toy.dicom.raw (), full.info = TRUE)
#' str (df)

dicom.browser <- function(dicom.raw.data, nbTAG = 0, stop.tag = "", stop.level= 0, 
                          full.info = FALSE, tag.dictionary = dicom.tag.dictionary ()){

  L <- .dicombrowser (dicom.raw.data, tag.dictionary, nbTAG = nbTAG,
                         stop_tag=stop.tag, stop_level= stop.level, 
                         full_info = full.info, verbose = FALSE)
  if (any(L=="not dicom compliant")) return(NULL)
  VRness <- L[length(L)]
  dicom.df <- do.call(rbind.data.frame, strsplit(L[-length(L)],";"))
  
  if (full.info) {
    colnames(dicom.df) <- c("tag", "VR", "endian","start", "stop","encaps.load",
                            "load.start","load.stop","tag.start")
    dicom.df [ , 7:9] <- suppressWarnings(lapply(dicom.df [ , 7:9] ,as.numeric ))
  } else {
    colnames(dicom.df) <- c("tag", "VR", "endian","start", "stop")
  }
  dicom.df [ , 4:5] <- suppressWarnings(lapply(dicom.df [ , 4:5] ,as.numeric ))
  lt.idx <- nrow(dicom.df)
  if (is.na(dicom.df[lt.idx, 5]) & dicom.df[lt.idx, 2] == "UN" & VRness=="1") 
    warning ( paste("Last decoded TAG",dicom.df[lt.idx, 1], 
                    "is unknown with unspecified length. DICOM parsing is unpredictable."))
  return(dicom.df)
}