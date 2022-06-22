############################################################################
#' Converting .Rdcm files to .xlsx files
#' @description A *.Rdcm file contains the list of contents, in dataframe form,
#' of the DICOM files of the same object. 
#' The \code{xlsx.from.Rdcm} function creates, from a *.Rdcm file, an Excel file, 
#' in which each page contains the dataframe representation of a DICOM file of
#' the same object.
#' @param Rdcm.filenames String vector, representing the *.Rdcm filenames to be 
#' converted.
#' @param dest.dirname  String vector of the same length as \code{Rdcm.filenames}, 
#' indicating the directory where the *.xlsx files will be created. 
#' @param txt.sep String. Used if \code{as.txt = TRUE}. Separator of the tag value elements.
#' @param txt.length Positive integer. Used if \code{as.txt = TRUE}. Maximum number 
#' of letters in the representation of the TAG value.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @return Returns a boolean vector, establishing the existence of the
#' created Excel files which have the same basenames as the *.Rdcm files.
#' @examples
#' # First, create a Rdcm file from toy.dicom.raw () to a temporary file for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "PM_rtplan", tmpdir = pat.dir, fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' dicom.to.Rdcm.converter (dcm.filename, pat.dir, update = TRUE)
#' file.remove (dcm.filename)
#' list.files (pat.dir)
#' 
#' # Creating an Excel file
#' Rdcm.filenames <- list.files (pat.dir, pattern = "[.]Rdcm$",
#'                               recursive = TRUE, full.names = TRUE)
#' xlsx.from.Rdcm (Rdcm.filenames)
#' list.files (pat.dir)
#' 
#' # Cleaning temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @import openxlsx
#' @export
xlsx.from.Rdcm <- function(Rdcm.filenames, dest.dirname = dirname(Rdcm.filenames),
                           txt.sep = "\\", txt.length = 100, 
                           tag.dictionary = dicom.tag.dictionary ()) {
  if (length(Rdcm.filenames) != length (dest.dirname)){
    warning ("Rdcm.filename and dest.dirname must have the same length.")
    return(rep(FALSE,Rdcm.filenames))
  }
  
  Rdcm.flag <- grepl("[.]Rdcm$" ,basename(Rdcm.filenames))
  Rdcm.f <- Rdcm.filenames[Rdcm.flag]
  dest.dirname <- dest.dirname[Rdcm.flag]
  if (length(Rdcm.f)==0) {
    warning ("Rdcm.filename are not Rdcm files.")
    return(Rdcm.flag)
  }
  xls.filenames<-rep("",length(Rdcm.f))
  for(f.idx in 1:length(Rdcm.f)){
    f <- Rdcm.f[f.idx]
    xls.filenames[f.idx] <- file.path (dest.dirname[f.idx], 
                                       paste(unlist(strsplit(basename(f),"[.]Rdcm$")),"xlsx", sep="."))
    
    if (file.exists(Rdcm.f[f.idx])){
      l <- load.Rdcm.raw.data (Rdcm.f[f.idx], address = TRUE, data=TRUE)
      if (l$from.dcm){
        wb <- createWorkbook ()
        for(sheet.idx in 1:length(l$data)) {
          
          db <- data.frame(TAG=names(l$data[[sheet.idx]]))
          dbtag <- sapply(db$TAG, function(t) rev(unlist(strsplit(t,"[ ]")))[1])
          
          db$VR <- l$address[[sheet.idx]]$VR
          db$VM <- tag.dictionary[match(dbtag, tag.dictionary$tag),"name"]
          db$loadsize <- sapply(l$data[[sheet.idx]], function (l_) length(l_))
          
          db$Value <- sapply(l$data[[sheet.idx]], function (l_) {
            if (length(l_)==0) return("")
            if (is.na(l_[1])) return("")
            lc <- as.character(l_)
            llc <- length(lc)
            cum <- cumsum (nchar(lc) + nchar(txt.sep))
            if (cum[llc] <= txt.length + nchar(txt.sep))  return(paste(lc[1:llc], collapse=txt.sep))
            idx <- which(cum<=txt.length-3)
            if (length(idx)==0) return (paste0(substr(lc,1,txt.length-3),"..."))
            n <- rev (idx)[1]
            return(paste(paste(lc[1:n], collapse=txt.sep),"...",sep=txt.sep))
          })
          
          # db$Value <- sapply(l$data[[sheet.idx]],function (l_) paste(l_[1:min(80,length(l_))], collapse="\\"))
          db[is.na(db)] <- ""
          db[db=="NA"] <- ""
          db[,c(1:3,5)] <- lapply(db[c(1:3,5)],as.character)
          db[,4]  <- as.integer(db[,4])
          
          if ("(0008,0005)" %in% db$TAG)Encoding(db$Value) <- db$Value[db$TAG==("(0008,0005)")]#"UTF-8"
          openxlsx::addWorksheet (wb, as.character (sheet.idx))
          openxlsx::setColWidths (wb, as.character (sheet.idx), cols = 1, widths = 30)
          openxlsx::setColWidths (wb, sheet = as.character (sheet.idx), cols=2:ncol(db), widths = "auto")
          openxlsx::writeDataTable (wb, as.character(sheet.idx), x=db)
        }
        
        openxlsx::saveWorkbook(wb, xls.filenames[f.idx], overwrite = TRUE)
      }
    }
  }
  Rdcm.flag[Rdcm.flag] <- file.exists(xls.filenames)
  return(Rdcm.flag)
}