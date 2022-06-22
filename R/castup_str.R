################################################################################
#' Cast of a character string
#' @description The \code{castup.str} function converts a word to upper case,
#' without accents and spaces.
#' @param st character string
#' @return Returns the ASCII//TRANSLIT transcription of the word \code{st},
#' without accents, spaces and in capitals.
#' @seealso \link[espadon]{castlow.str}.
#' @examples
#' castup.str (st = c("Right eye", "Left_Lung", "Right-Lung"))
#' @export
castup.str <- function (st) {
  dum <- iconv (st, to="ASCII//TRANSLIT")
  return (toupper (gsub("[[:space:],_]", "", dum)))
}
