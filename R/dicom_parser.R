#' Conversion of DICOM raw data into a dataframe or a list of DICOM TAG information
#' @description The \code{dicom.parser} function creates a dataframe or a list from 
#' DICOM raw data. The created dataframe or list provides information about the 
#' content of the DICOM TAGs included in the raw data.
#'
#' @param dicom.raw.data Raw vector, representing the binary extraction of the DICOM file.
#' @param as.txt Boolean. If \code{as.txt = TRUE}, the function returns a 
#' dataframe, a list otherwise.
#' @param nested.list Boolean. Only used if \code{as.txt = FALSE}. If 
#' \code{nested.list = FALSE}, the returned list consists  of nested lists.
#' @param try.parse Boolean. If \code{TRUE}, the tag with unknown DICOM VR 
#' (value representation) is converted into string if possible.
# @param txt.sep String. Used if \code{as.txt = TRUE}. See Note.
# @param txt.length Positive integer. Used if \code{as.txt = TRUE}. See Note.
#' @param txt.sep String. Used if \code{as.txt = TRUE}. Separator of the tag value elements.
#' @param txt.length Positive integer. Used if \code{as.txt = TRUE}. Maximum number 
#' of letters in the representation of the TAG value.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.

# @note If \code{as.txt = TRUE}, and if the TAG contains a non-ASCII value, 
# then it will be represented by the concatenation of \code{txt.length} at most 
# values, separated by \code{txt.sep}.

#' @return Returns a list of elements or a dataframe, depending on  \code{as.list}. 
#' @return If it returns a dataframe, the columns are names TAG, VR (value representation), 
#' VM (value multiplicity), loadsize and Value. The field \code{$Value} is a string 
#' representation of the true value.
#' @return If it returns a list, each of its elements, named by a TAG, is either 
#' a vector or a string, depending of the TAG included in \code{dicom.raw.data}.

#' @seealso \link[espadon]{dicom.raw.data.loader}, \link[espadon]{dicom.tag.parser} 


#' @examples
#' # content of the dummy raw data toy.dicom.raw (), as a list.
#' L <- dicom.parser (toy.dicom.raw (), as.txt = FALSE)
#' str(L[40:57])
#' 
#' L <- dicom.parser (toy.dicom.raw (), as.txt = FALSE, nested.list = TRUE)
#' str(L[40:45])
#' 
#' # content of the dummy raw data toy.dicom.raw (), as a dataframe.
#' L <- dicom.parser (toy.dicom.raw (), as.txt = TRUE)
#' str (L)

#' @export
dicom.parser <- function (dicom.raw.data, as.txt=TRUE, nested.list = FALSE, try.parse = FALSE, txt.sep = "\\", 
                          txt.length = 100, tag.dictionary = dicom.tag.dictionary ()){
  if (!is.raw(dicom.raw.data)) {
    warning ("dicom.raw.data should be raw data :).")
    return (NULL)
  }
  dicom.df <- dicom.browser(dicom.raw.data, tag.dictionary=tag.dictionary)
  if (is.null(dicom.df)) return (NULL)
  L <- lapply(1:nrow(dicom.df), function (idx) dicom.tag.parser (dicom.df$start[idx],dicom.df$stop[idx],
                                                         dicom.df$VR[idx], dicom.df$endian[idx],
                                                         dicom.raw.data,try.parse=try.parse))
  names(L) <- dicom.df$tag
  
  if (!as.txt & !nested.list) return(L)
  if (nested.list) {
    n <- names(L)
    encaps <- as.numeric(sapply(n,function(st) length(unlist(strsplit(st," ")))))
    parent.encaps.idx <- sapply( 1:length(encaps), function(i)
      rev(which((encaps[1:i] == encaps[i]-1)))[1])
    seq.idx <- !is.na(match(1:length(encaps),unique (parent.encaps.idx[!is.na(parent.encaps.idx)])))
    
    l <- lapply(1: length(L),function(i) NULL)
    names(l)  <- sapply(names(L), function(str) rev(unlist(strsplit(str," ")))[1])
    l[!seq.idx] <- L[!seq.idx]
    
    for (i in max(encaps):2){
      ei <- which(encaps==i)
      li <- unique(parent.encaps.idx[ei])
      for (ej in li)
        l[[ej]] <-   l[ei[ej==parent.encaps.idx[ei]]]
    }
    l[encaps>1] <- NULL
    return (l)
  }
  
  db <- dicom.df[, c(1:2)]
  colnames( db) <- c ("TAG", "VR")
  dbtag <- sapply(db$TAG, function(t) rev(unlist(strsplit(t,"[ ]")))[1])
  db$VM <- tag.dictionary[match(dbtag,tag.dictionary$tag),"name"]
  
  db$loadsize <- sapply(L, function (l) length(l))
  
  encoding <- L[["(0008,0005)"]]
  if(!is.null(encoding)) {
    conv_idx <- grep(paste0("^",tolower (gsub("[[:space:],_,-]", "", encoding)),"$"),
                     tolower (gsub("[[:space:],_,-]", "", iconvlist())))
    if (length(conv_idx)>0) {encoding <- iconvlist()[conv_idx[1]]
    } else encoding <- NULL
  }
  
  db$Value <- sapply(L, function (l) {
    if (length(l)==0) return("")
    if (is.na(l[1])) return("")
    if (is.character(l) & !is.null(encoding)){ lc <- iconv(l, encoding)
    } else {lc <- as.character(l)}
    llc <- length(lc)
    cum <- cumsum (nchar(lc) + nchar(txt.sep))
    if (cum[llc] <= txt.length + nchar(txt.sep))  return(paste(lc[1:llc], collapse=txt.sep))
    idx <- which(cum<=txt.length-3)
    if (length(idx)==0) return (paste0(substr(lc,1,txt.length-3),"..."))
    n <- rev (idx)[1]
    return(paste(paste(lc[1:n], collapse=txt.sep),"...",sep=txt.sep))
    # if (txt.length < length(l)) return(paste(paste(l[1:txt.length], collapse=txt.sep),"...",sep=txt.sep))
    # return (paste(l, collapse=txt.sep))
    })
  
  db[is.na(db)] <- ""
  db[db=="NA"] <- ""
  db[,c(1:3,5)] <- lapply(db[c(1:3,5)],as.character)
  db[,4]  <- as.integer(db[,4])
  # Encoding(db$Value) <- "UTF-8"
  
  return(db)
}
