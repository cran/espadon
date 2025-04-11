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
#' \item encaps.load : If the tag contains nested data, this column gives 
#' the number of bytes remaining until the end of the nesting. If there 
#' are several levels of nesting, these numbers are collapsed and separated 
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
  if (!is.raw(dicom.raw.data)) stop ("dicom.raw.data should be raw data :).")

  L <- .dicombrowser (dicom.raw.data, tag.dictionary, nbTAG = nbTAG,
                         stop_tag=stop.tag, stop_level= stop.level, 
                         full_info = full.info, verbose = FALSE)
  if (any(L$error!="")) {
    warnings(L$error)
    return(NULL)
  }

  lt.idx <- nrow(L$db)
  if (is.na(L$db[lt.idx, 5]) & L$db[lt.idx, 2] == "UN" & L$VRness=="1") 
    warning ( paste("Last decoded TAG",L$db[lt.idx, 1], 
                    "is unknown with unspecified length. DICOM parsing is unpredictable."))
  return(L$db)
}