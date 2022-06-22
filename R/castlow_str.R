################################################################################
#' Cast of a character string
#' @description The \code{castlow.str} function converts a word to lowercase,
#' without accents and spaces.
#' @param st character string
#' @return Returns the ASCII//TRANSLIT transcription of the word \code{st},
#' without accents, spaces and in lowercase letters.
#' @seealso \link[espadon]{castup.str}.
#' @examples
#' castlow.str (st = c("Right eye", "Left_Lung", "Right-Lung"))

#' @export
castlow.str <- function (st) {
  dum <- iconv (st, to="ASCII//TRANSLIT")
  return (tolower (gsub("[[:space:],_]", "", dum)))
}