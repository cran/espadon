#' DICOM anonymizer
#' @description the \code{dicom.raw.data.anonymizer} function anonymizes 
#' \code{dicom.raw.data}.
#' @param dicom.raw.data Raw vector, representing the binary extraction of the 
#' DICOM file.
#' @param offset Integer, default to 0. Each date of the DICOM will be shifted 
#' by this offset expressed in days.
#' @param new.PIN Character string, representing the PIN remplacing the old one.
#' @param reset.private.tag Boolean, if \code{TRUE}, the value of tags that are 
#' not in the \code{tag.dictionary} is removed.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, whose structure it must keep. This 
#' dataframe is used to parse DICOM files.
#' @return Returns an anonymyzed raw vector. See Note.
#' @note The raw data is anonymized as follows:  
#' \itemize{
#' \item Each date of the DICOM file will be shifted by \code{offset} expressed in days.
#' \item Each patient's name, and patient'ID are remplaced by \code{new.PIN}
#' \item All other patient data are deleted, except age, weight, height, gender 
#' and shifted birthday.
#' \item All address, phone, physician, operator, author, reviewer, service.
#' \item If \code{reset.private.tag = TRUE}, the values of the tags not contained in the \code{tag.dictionary} are 
#' deleted.
#' }
#' @export
#' @examples
#' # pseudomization of the dummy raw data toy.dicom.raw ()
#' an.raw.data <- dicom.raw.data.anonymizer (toy.dicom.raw (), offset = -2)
#' data <- dicom.parser (toy.dicom.raw ())
#' an.data <- dicom.parser (an.raw.data) 
#'
#' # Checking for differences
#' flag.dif <- data$Value !=  an.data$Value  
#' df <- cbind (data[flag.dif, c ("VM","Value")], an.data[flag.dif, "Value"])      
#' colnames (df) <-   c ("VM","old Value","new Value")    
#' df                                     
#'
#' # save data in a the new file
#' #############################
#' # new.file.name <- "an.dcm"
#' # zz <- file (new.file.name, "wb")
#' # writeBin (an.raw.data, zz, size = 1)
#' # close (zz)
dicom.raw.data.anonymizer <- function( dicom.raw.data, offset = 0 , new.PIN = "Anonymous ",
                                       reset.private.tag = FALSE,
                                       tag.dictionary = dicom.tag.dictionary ()){
  if (nchar(new.PIN)%%2 != 0) new.PIN <- paste(new.PIN," ",sep = "")
  
  dicom.df <- dicom.browser (dicom.raw.data, full.info = TRUE, tag.dictionary = tag.dictionary)
  if (is.null(dicom.df)) {
    warning ("not dicom compliant.")
    return(FALSE)
  }

  last.tag <-sapply(dicom.df$tag, function(t) rev(unlist(strsplit(t, "[ ]")))[1])
  
  m <- match(last.tag, tag.dictionary$tag)
  name<- tag.dictionary$name[m]
  not.item <- grepl("[)]$",last.tag)
  has.load <- !is.na(dicom.df$start) & !is.na(dicom.df$stop)
  if (reset.private.tag){ 
    rpt <- which(is.na(name) & not.item & !grepl ('0000[)]$',last.tag) & has.load)
  }else {
    rpt <- integer(0)
    }
  
  res.idx <- which (grepl ('physician|operator|author[ ]|author$|reviewer',tolower (name)) & has.load)
  pat.idx <- which (grepl ('^[(]0010,',tolower (last.tag)) & has.load)
  reset.idx <- unique (sort (c (res.idx,#[grepl('requesting|name|address|phone',tolower (name[res.idx]))],
                                # which (grepl ('undocumented',tolower (name)) & has.load), 
                                # which (dicom.df$VR=="UN" & has.load),
                                rpt,
                                pat.idx[-grep ('0000[)]$|0010[)]$|0020[)]$|0030[)]$|0040[)]$|1010[)]$|1020[)]$|1030[)]$',last.tag[pat.idx])],
                                which (last.tag=="(4008,0119)" & has.load),
                                which (last.tag=="(4008,011A)"& has.load))))
  
  pat.idx <- which(grepl("[(]0010,0020[)]$|[(]0010,0010[)]$", dicom.df$tag) & has.load)
  
  identity.removed <- which(last.tag=="[(]0012,0062[)]")
  
  DA.idx <- which (dicom.df$VR=="DA" & has.load)
  DT.idx <- which (dicom.df$VR=="DT" & has.load)
  
  modify.idx <- sort(c(pat.idx,reset.idx,identity.removed), decreasing = TRUE)
  
  
  # anonymise date
  if (length(DA.idx)>0) {
    new.DA <- sapply(DA.idx, function(idx) {
      as.character(format(as.Date(dicom.tag.parser (dicom.df$start[idx],dicom.df$stop[idx],
                                                    dicom.df$VR[idx], dicom.df$endian[idx],
                                                    dicom.raw.data), format="%Y%m%d") + offset, format="%Y%m%d") )})
    for (idx in 1:length(DA.idx)) dicom.raw.data[dicom.df$start[DA.idx[idx]]:dicom.df$stop[DA.idx[idx]]] <- charToRaw(new.DA[idx])
  }
  
  if (length(DT.idx)>0) {
    new.DT <- sapply(DT.idx, function(idx) {
      DT <- dicom.tag.parser (dicom.df$start[idx],dicom.df$stop[idx],
                              dicom.df$VR[idx], dicom.df$endian[idx],
                              dicom.raw.data)
      DA <-substr(DT,1,8) 
      if (nchar(DT)>8) {reste <- substr(DT,9,nchar(DT)) 
      } else {reste <- ""}
      paste(as.character(format(as.Date(DA, format="%Y%m%d") + offset, format="%Y%m%d")),reste,sep="")
    })
    for (idx in 1:length(DT.idx)) dicom.raw.data[dicom.df$start[DT.idx[idx]]:dicom.df$stop[DT.idx[idx]]] <- charToRaw(new.DT[idx])
  }
  
  if (length(modify.idx)>0) {
    
    for (idx in modify.idx) {
      new.raw <- NULL
      part1 <- NULL
      part3 <- NULL
      if (idx %in% pat.idx) {
        new.raw <- charToRaw(new.PIN)
      } else if (idx %in% identity.removed) {new.raw <- charToRaw("YES ")}
      if (dicom.df$start[idx]> 1) part1 <- dicom.raw.data[1:(dicom.df$start[idx]-1)]
      if (length(dicom.raw.data)> dicom.df$stop[idx]) part3 <- dicom.raw.data[(dicom.df$stop[idx] + 1): length(dicom.raw.data)]
      #transformer la partie du dicom.raw.data
      dicom.raw.data <- c(part1, new.raw, part3)
      
      new.raw.load <- packBits(intToBits(as.raw(length (new.raw))), type="raw")
      dicom.raw.data.loadsize <- dicom.df$load.stop[idx]-dicom.df$load.start[idx]+1
      # if (dicom.raw.data.loadsize==4){
      #   if (dicom.df$endian[idx]=="little") {
      #     dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <-new.raw.load
      #   } else { dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <-rev(new.raw.load)}
      # } else {
      if (dicom.df$endian[idx]=="little") {
        dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <-new.raw.load[1:dicom.raw.data.loadsize]
      } else {
        dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <-rev(new.raw.load[1:dicom.raw.data.loadsize]) 
      }
      # }
      
      
      if (nchar(dicom.df$encaps.load[idx])>0){
        to.substract <- dicom.df$stop[idx]-dicom.df$start[idx]+1 - length(new.raw)
        level.to.change <- suppressWarnings (as.numeric(unlist(strsplit(dicom.df$encaps.load[idx],'[ ]'))))
        if (length(level.to.change)>0) {
          tag.to.change <- unlist(strsplit(dicom.df$tag[idx],'[ ]'))
          tag.to.change <- sapply (1:length(level.to.change), function (i) paste(tag.to.change[1:i], collapse=" "))
          tag.to.change.idx <- match(tag.to.change,dicom.df$tag)
          tag.to.change.idx  <- tag.to.change.idx[!is.na(level.to.change)]
          #--> reset idx
          for (i in tag.to.change.idx){
            size <- dicom.df$load.stop[i]-dicom.df$load.start[i]+1
            raw.load <- packBits(intToBits(readBin (dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]], what="int", n= 1,
                                                    size = size,
                                                    endian = dicom.df$endian[i]) - to.substract), type="raw")
            
            if (size==4){
              if (dicom.df$endian[i]=="little") {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-raw.load
              } else { dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-rev(raw.load)}
            } else {
              if (dicom.df$endian[i]=="little") {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-raw.load[1:2]
              } else {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-rev(raw.load[1:2]) }
            }
            
          }
        }
      }
    }
  }
  
  

  idx.group  <- grep("0000[)]$", dicom.df$tag)
  tag.group <- sapply(dicom.df$tag[idx.group], function (t) substr(t,1,nchar(t)-5))
  tag.reg <-  paste0("^", gsub(")","[)]", gsub("(","[(]",tag.group, fixed = TRUE), fixed = TRUE))

  impacted <- sapply(tag.reg,function(t) any(grepl(t, dicom.df$tag[modify.idx])))
  if (any(impacted)){
    dicom.df <- dicom.browser(dicom.raw.data,full.info = TRUE)

    length.group <- sapply(tag.reg[impacted], function(t){
      idx <- range(grep(t,dicom.df$tag))+1
      e <- ifelse(idx[2]<= nrow(dicom.df),dicom.df$tag.start[idx[2]], length(dicom.raw.data)+1)
      e-dicom.df$tag.start[idx[1]]})

    idx.group_ <-idx.group[impacted]
    for (i in 1:length(idx.group_)){
      size <- dicom.df$stop[idx.group_[i]]-dicom.df$start[idx.group_[i]]+1
      raw.load_ <- packBits(intToBits(length.group[i]), type="raw")

      if(!is.na(size)){
        if (size==4){
          if (dicom.df$endian[i]=="little") {raw.load <- raw.load_
          } else {raw.load <-rev(raw.load_)}
        } else {
          if (dicom.df$endian[i]=="little") {raw.load<- raw.load_[1:2]
          } else { raw.load<-rev(raw.load_[1:2])}
        }

        dicom.raw.data[dicom.df$start[idx.group_[i]]:dicom.df$stop[idx.group_[i]]] <-raw.load
        
      }
    }
  }

  return(dicom.raw.data)
}

