#' DICOM TAG dictionary
#' @description The \code{dicom.tag.dictionary} function gives the dictionary of 
#' tags used by default in the \pkg{espadon} package.
#' @param add.dict Vector of the list of additional dictionaries. Put to NULL, if 
#' no additional dictionary is requested.
#' @return Returns a 3-column dataframe, describing the VR (value representation)
#' and the name of each DICOM TAG. 
#' @return This dataframe is the fusion of the "nema.tag" dictionary, provided 
#' by \emph{nema} \strong{\[1\]}, with the dictionaries defined in the 
#' \code{add.dict} vector:
#' \itemize{
#' \item "raysearch.tag" dictionary is provided by \emph{RaySearch laboratories} \strong{\[2\]}
#' }

#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{nema}{espadon}
#' @references \strong{\[2\]} \insertRef{raysearch}{espadon}
#' @examples
#' str (dicom.tag.dictionary ())
#' str (dicom.tag.dictionary (NULL))
#' @export
dicom.tag.dictionary <- function (add.dict = c("raysearch.tag")) {
  db <- dicomtag
  if ("raysearch.tag" %in% add.dict) db <- rbind(db, raystationdt)
  db <- db[order(db$tag), ]
  row.names(db) <- NULL
  return(db)
}
