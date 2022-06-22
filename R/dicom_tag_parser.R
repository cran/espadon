#' DICOM TAG parser
#' @description the \code{dicom.tag.parser} function decodes the content between 
#' two DICOM raw data addresses.
#' @param start Positive integer. Index of the first raw data to parse in the 
#' \code{dicom.raw.data}.
#' @param stop Positive integer. Index of the last raw data to parse in the 
#' \code{dicom.raw.data}.
#' @param VR Character string, representing the value representation of DICOM 
#' data. See DICOM standard.
#' @param endian Character string, equal to "little" or "big". 
#' @param dicom.raw.data Raw vector, representing the binary extraction of the 
#' DICOM file.
#' @param try.parse Boolean. If \code{TRUE}, the value, with an undocumented VR,
#' is considered, as far as possible, as a string.
#'
#' @return Returns the \code{dicom.raw.data} content between the addresses 
#' \code{start} and \code{stop}. Depending on the representation of the value 
#' (\code{VR}), it can be a character string or a numerical vector.
#' @export
#'
#' @examples
#' # creation of the toy.dicom.raw () addresses dataframe:
#' df <- dicom.browser (toy.dicom.raw ())
#' 
#' # search for modality of toy.dicom.raw ()
#' idx <- grep ("^[(]0008,0060[)]$", df$tag)
#' modality  <- dicom.tag.parser (df$start[idx], df$stop[idx], df$VR[idx], 
#'                                df$endian[idx], toy.dicom.raw ())
#' modality


dicom.tag.parser <- function (start, stop, VR, endian, dicom.raw.data, try.parse = FALSE) {
  # .dicom.parse.tag <- function (idx, dicom.df, dicom.raw.data) {
  s <- start[1] #dicom.df$start[idx[1]]
  e <- stop[1] #dicom.df$stop[idx[1]]
  
  if (is.na(s) | is.na(e)) return (NA)
  m <- dicom.raw.data[s:e]
  if (nchar (VR[1])>2)  return (m)
  
  switch(VR[1],
         "UL" = {
           v =readBin (m, what="int", n= (e-s+1)/4, size = 4, endian = endian[1]) 
           n <- v[v<0]
           v[v<0]= 4294967296 + n
           return (v)},
         "US" = {readBin (m, what="int", n= (e-s+1)/2, size = 2, signed =FALSE, endian = endian[1]) },
         "OW" = {readBin (m, what="int", n= (e-s+1)/2, size = 2, signed =FALSE, endian = endian[1]) },
         "OB" = {readBin (m, what="integer", n= (e-s+1), size = 1, signed =FALSE) },
         # "UN" = {readBin (m, what="int", n= (e-s+1), size = 1, signed =FALSE, endian = endian[1]) },
         "SL" = {readBin (m, what="int", n= (e-s+1)/4, size = 4, endian = endian[1])},
         "SS" = {readBin (m, what="int", n= (e-s+1)/2, size = 2, signed =TRUE, endian = endian[1])},
         "FD" = {readBin (m, what="double", n= (e-s+1)/8, size = 8, endian = endian[1])},
         "OD" = {readBin (m, what="double", n= (e-s+1)/8, size = 8, endian = endian[1])},
         "FL" = {readBin (m, what="double", n= (e-s+1)/4, size = 4, endian = endian[1])},
         "OF" = {readBin (m, what="double", n= (e-s+1)/4, size = 4, endian = endian[1])},
         "AT" = {if (endian[1]=="little") m <- m[c(2,1,4,3)]
         return (toupper(paste("(", m[1], m[2],",", m[3],m[4],")", sep="")))},
         "SQ" = {return(NULL)},
         "00" = {return(NULL)},
         "UN" = {
           if (try.parse){
             id0 <- which(m==as.raw(0))[1]
             if (is.na(id0)) return (rawToChar (m))
             if (id0 == 1) return("")
             return (rawToChar (m[1:(id0-1)]))
           } else {
             return (m)
           }
         },
         {
           id0 <- which(m==as.raw(0))[1]
           if (is.na(id0)) return (rawToChar (m))
           if (id0 == 1) return("")
           return (rawToChar (m[1:(id0-1)])) 
         }
  )
  
  
}