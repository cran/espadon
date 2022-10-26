#' Change TAG value in DICOM raw data
#' @description The \code{dicom.set.tag.value} function changes, in the DICOM 
#' raw data, the values of the TAG whose VR is a string of characters.
#' @param dicom.raw.data Raw vector, representing the binary extraction of the DICOM file.
#' @param tag String vector, representing the list of tags whose value is to be 
#' changed. See note 1.
#' @param tag.value String vector,representing the list of new tag values.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @note 1- The list of tags included in the DICOM file are given by the first columns 
#' of the dataframe provided by the functions \link[espadon]{dicom.browser} and 
#' \link[espadon]{dicom.parser}.
#' @note 2-  The \code{dicom.set.tag.value} function may take some processing time. 
#' To minimize this time, it is recommended to prepare in advance all the tags to 
#' be modified, and use the \code{dicom.set.tag.value} function only once, as shown in 
#' the example.
#' @return Returns a raw vector, with new tag values.
#' @export
#' @examples

#' # change the value of tags "(0010,0010)" and "(0010,0020)" in the
#' # dummy raw data toy.dicom.raw ()
#' new.raw.data <- dicom.set.tag.value (toy.dicom.raw (), 
#'                                      tag =  c ("(0010,0010)", "(0010,0020)"),
#'                                      tag.value = c ("unknown", "000001"))
#' # change control 
#' data <- dicom.parser (new.raw.data) 
#' data[data$TAG %in% c ("(0010,0010)", "(0010,0020)"), ]
#' 
#' # save data in a the new file
#' #############################
#' # new.file.name <- "new.dcm"
#' # zz <- file (new.file.name, "wb")
#' # writeBin (new.raw.data  , zz, size = 1)
#' # close (zz)

dicom.set.tag.value <- function (dicom.raw.data, tag, tag.value, 
                                 tag.dictionary = dicom.tag.dictionary ()) {
  # now <- Sys.time()
  dicom.df <- dicom.browser (dicom.raw.data,   full.info = TRUE, tag.dictionary = tag.dictionary)
  # Sys.time()-now
  if (is.null(dicom.df)) {
    warning ("not dicom compliant.")
    return (dicom.raw.data)
  }
  if (length(tag.value)==1) tag.value <- rep (tag.value, length(tag))
  if (length(tag)!=length(tag.value)) {
    warning ("tag and tag.value must have same length or tag.value length must be 1.")
    return (dicom.raw.data)
  }
  modify.idx  <- as.numeric(sapply (tag, function (t) {
    idx <- match(t,dicom.df$tag)
    if (length(idx)==0) return(NA)
    return(idx)
  }))
  
  modify.idx.flag <- !is.na(modify.idx)
  if (any(!modify.idx.flag)) {
    for(w.idx in which(!modify.idx.flag)) message (paste ("tag", tag[w.idx],"does not exist."))
  }
  
  tag_ <- tag[modify.idx.flag]
  tag.value_ <- tag.value[modify.idx.flag]
  modify.idx <- modify.idx[modify.idx.flag]
  
  if (length(modify.idx)==0) return (dicom.raw.data)
  
  flag.ascii <- !is.na(match (dicom.df$VR[modify.idx], c("AE", "AS", "CS", "DA", 
                                                         "DS", "DT", "IS", "LO", 
                                                         "LT", "PN", "SH", "ST",
                                                         "TM", "UI", "UN", "UT")))
  if (any(!flag.ascii)) {
    for(w.idx in which(!flag.ascii)) message (paste ("tag",  tag[w.idx],
                                                     "does not have an ASCII VR and is not modified."))
  } 
  tag_ <- tag_[flag.ascii]
  tag.value_  <- tag.value_[flag.ascii]
  modify.idx <- modify.idx[flag.ascii]
  
  if (length(modify.idx)==0) return (dicom.raw.data)
    
  conformity.flag<- rep(TRUE, length(modify.idx))
  for (conf.idx in 1: length (conformity.flag)) {
    switch(dicom.df$VR[modify.idx[conf.idx]],
    "AE" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- !is.na(match(new.raw, as.raw(c(10, 12, 13, 27, 92))))

      if (length(new.raw)>16 | any(m)) {
        message (paste ("tag", tag_[conf.idx], "is not AE compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      } 
    },
    "AS" = {
      new.raw <- charToRaw (tag.value_ [conf.idx])
      m <- match(new.raw[1:3], as.raw(48:57))

      if (length(new.raw)!=4 | any(is.na(m)) | !(new.raw[4] %in% charToRaw ("DWMY"))){
        message (paste ("tag", tag_[conf.idx], "is not AS compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
        }
    },
    "CS" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- match(new.raw, as.raw(c(65:90,48:57,95,32)))
      
      if (length(new.raw)>16 | any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not CS compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
      tag.value_ [conf.idx] <- rawToChar(new.raw)
    },
    "DA" = {
      new.raw <- charToRaw (tag.value_ [conf.idx])
      m <- match(new.raw[1:8], as.raw(48:57))
      if (length(new.raw)!=8 | any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not DA compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "DS" = {
      new.raw <- charToRaw (tag.value_ [conf.idx])
      m <- match(new.raw,as.raw(32))
      new.raw <- new.raw[is.na(m)]
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- match(new.raw, charToRaw("0123456789+-Ee. \\"))
      
      if (any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not DS compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "DT" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- match(new.raw, charToRaw("0123456789+-. "))
      if (length(new.raw)>26 | any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not DT compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "IS" = {
      mem <- suppressWarnings(as.numeric(tag.value_ [conf.idx]))
      if (is.na(mem)){
        message (paste ("tag", tag_[conf.idx], "is not IS compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      } else if ((mem < -2^31) | (mem > +2^31 - 1 )){
        message (paste ("tag", tag_[conf.idx], "must be betwenn -2^31 and (2^31-1). This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      } else {
        new.raw <- charToRaw (tag.value_ [conf.idx])
        m <- match(new.raw,as.raw(32))
        new.raw <- new.raw[is.na(m)]
        if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
        tag.value_ [conf.idx] <- rawToChar(new.raw)
        
        m <- match(new.raw, charToRaw("0123456789+- "))
        if (length(new.raw)>12 | any(is.na(m))){
          message (paste ("tag", tag_[conf.idx], "is not IS compliant. This TAG is not modified."))
          conformity.flag[conf.idx] <- FALSE
        }
      }
    },
    "LO" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- !is.na(match(new.raw,as.raw(c(0:26,28:31,92))))
      if (length(new.raw)>64 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not LO compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "LT" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
      
      if (length(new.raw)>64 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not LT compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    }, 
    "PN" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      le <- 0
      if (nchar(tag.value_ [conf.idx])>0) le <- max(nchar(unlist(strsplit(tag.value_ [conf.idx],"^", fixed=TRUE))))
      m <- !is.na(match(new.raw, as.raw(c(10, 12, 13, 92))))
      if (le>64 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not PN compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "SH" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      m <- !is.na(match(new.raw, as.raw(c(0:26,28:31,92))))
      if (length(new.raw)>16 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not SH compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "ST" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
      if (nchar(tag.value_ [conf.idx])>1024 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not SH compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    },
    "TM" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- match(new.raw, charToRaw("0123456789. "))
      if (length(new.raw)>16 | any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not TM compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    }, 
    "UI" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(0))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- match(new.raw, as.raw(c(0,46,48:57)))
      if (length(new.raw)>64 | any(is.na(m))){
        message (paste ("tag", tag_[conf.idx], "is not UI compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    }, 
    "UN" = {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
    }, 
    "UT"= {
      new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
      if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
      tag.value_ [conf.idx] <- rawToChar(new.raw)
      
      m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
      if (length(new.raw)>2^32-2 | any(m)){
        message (paste ("tag", tag_[conf.idx], "is not UT compliant. This TAG is not modified."))
        conformity.flag[conf.idx] <- FALSE
      }
    })
  }
  
  
  tag_ <- tag[conformity.flag]
  tag.value_  <- tag.value_[conformity.flag]
  modify.idx <- modify.idx[conformity.flag]
  
  if (length(modify.idx)==0) return (dicom.raw.data)
  
  order.idx <- order(modify.idx, decreasing = TRUE)
  modify.idx <- modify.idx[order.idx]
  tag_ <- tag_[order.idx]
  tag.value_ <- tag.value_[order.idx] 
  
  for (m.idx in 1:length(modify.idx)) {
    idx <- modify.idx[m.idx]
    new.raw <- charToRaw(tag.value_[m.idx])
    part1 <- NULL
    part3 <- NULL
    # if (dicom.df$start[idx]> 1) part1 <- dicom.raw.data[1:(dicom.df$start[idx]-1)]
    if (dicom.df$load.stop[idx]> 0) part1 <- dicom.raw.data[1:(dicom.df$load.stop[idx])]
    if (nrow(dicom.df)>idx) part3 <- dicom.raw.data[dicom.df$tag.start[idx+1]: length(dicom.raw.data)]
    # if (length(dicom.raw.data)> dicom.df$stop[idx]) part3 <- dicom.raw.data[(dicom.df$stop[idx] + 1): length(dicom.raw.data)]
    #transformer la partie du dicom.raw.data
    dicom.raw.data <- c(part1, new.raw, part3)
    
    new.raw.load <- packBits(intToBits(length (new.raw)), type="raw")
    dicom.raw.data.loadsize <- dicom.df$load.stop[idx]-dicom.df$load.start[idx]+1
    
    if (dicom.df$endian[idx]=="little") {
      dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <- new.raw.load[1:dicom.raw.data.loadsize]
    } else {
      dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <- rev(new.raw.load[1:dicom.raw.data.loadsize]) 
    }
    # }
    
    
    if (nchar(dicom.df$encaps.load[idx])>0){
      if (!is.na(dicom.df$start[idx])){
        to.substract <- dicom.df$stop[idx]-dicom.df$start[idx]+1 - length(new.raw)
      } else {
        to.substract <- - length(new.raw)
      }
      level.to.change <- suppressWarnings (as.numeric(unlist(strsplit(dicom.df$encaps.load[idx],'[ ]'))))
      if (length(level.to.change)>0) {
        tag.to.change <- unlist(strsplit(dicom.df$tag[idx],'[ ]'))
        tag.to.change <- sapply (1:length(level.to.change), function (i) paste(tag.to.change[1:i], collapse=" "))
        tag.to.change.idx <- match(tag.to.change,dicom.df$tag)
        tag.to.change.idx  <- tag.to.change.idx[!is.na(level.to.change)]
        
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
  
  return (dicom.raw.data)
}