######################################################################################
#' Regions of Interest (RoI) indices
#' @description The select.names function allows you to select words from a
#' vector of words, according to several criteria, eliminating spaces and case.
#' @param names Words vector
#' @param roi.name Vector of words to compare to \code{names}. By default 
#' \code{roi.name = NULL}. See Details
#' @param roi.sname Vector of words or parts of words to compare. By default 
#' \code{roi.sname = NULL}. See Details
#' @param roi.idx Index vector. By default \code{roi.idx = NULL}. See Details.
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all 
#' \code{NULL}, then all RoI are selected.
# @note The returned index vector is sorted from smallest to largest.
#' @return Returns the indices of the elements of the word vector \code{names} 
#' satisfying one or more of the following conditions:
#' \itemize{
#' \item ASCII // TRANSLIT transcriptions, without spaces, of \code{names} and 
#' \code{roi.name}, are identical.
#' \item ASCII // TRANSLIT transcriptions, without spaces of \code{roi.sname}
#' are identical to part of ASCII // TRANSLIT transcriptions, without spaces of
#' \code{names}.
#' \item \code{names} indices belong to the index vector \code{roi.idx}.
#' }
#' @examples
#' # loading patient objects
#' names <- c ("Eye left", "EyeR", "OPTICAL nerve L", "opical nervR", "chiasma")
#'
#' # RoI selection.
#' select.names (names = names, roi.name = c("eye left", "eye right"))
#' select.names (names = names, roi.sname = c("eye", "ner"))
#' select.names (names = names, roi.idx = 4:9)

#' @export
select.names <- function (names, roi.name = NULL, roi.sname = NULL,
                          roi.idx = NULL) {
  if (is.null(names)) return(NULL)
  if (all (is.null (roi.name)) & all (is.null (roi.sname)) & all (is.null (roi.idx))) return (1:length (names))
  s.idx <- c()
  
  if (!all (is.null (roi.idx))) s.idx <- c(s.idx, roi.idx[roi.idx %in% 1:length (names)])
  
  a <- iconv (names, to="ASCII//TRANSLIT")
  a <- toupper (gsub("[[:space:],_]", "", a))
  b <- iconv (roi.name, to="ASCII//TRANSLIT")
  b <- toupper (gsub("[[:space:],_]", "", b))
  # if (!all (is.null (roi.name))) s.idx <- c(s.idx, which (a %in% b))
  if (!all (is.null (roi.name))) s.idx <- c(s.idx, match(b,a))
  s.idx <- s.idx[!is.na(s.idx)]
  
  b <- iconv (roi.sname, to="ASCII//TRANSLIT")
  b <- toupper (gsub("[[:space:],_]", "", b))
  if (!all (is.null (roi.sname))) {
    # result <- sapply (b, function (test) return (grepl (pattern=test, a)))
    # if (any(result)) s.idx <- c(s.idx, which (apply (result, 1, function (v) sum (v) > 0)))
    result <- unlist(lapply (b, function (test) return (grep (pattern=test, a))))
    if (length(result) > 0) s.idx <- c(s.idx, result)
  }
  # s.idx <- sort (unique (s.idx))
  s.idx <-  (unique (s.idx))
  if (length(s.idx)==0) return (NULL)
  return (s.idx)
}
