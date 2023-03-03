#' Change TAG value in DICOM raw data
#' @description The \code{dicom.set.tag.value} function changes, in the DICOM 
#' raw data, the values of the TAG whose VR is a string of characters.
#' @param dicom.raw.data Raw vector, representing the binary extraction of the DICOM file.
#' @param tag String vector, representing the list of tags whose value is to be 
#' changed. See note 1.
#' @param tag.value String vector,representing the list of new tag values.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @param ... Additional arguments \code{dicom.browser} when previously calculated by 
#' \link[espadon]{dicom.browser}  with argument \code{full.info = TRUE}.
#' @note 1- The list of tags included in the DICOM file are given by the first columns 
#' of the dataframe provided by the functions \link[espadon]{dicom.browser} and 
#' \link[espadon]{dicom.parser}.
#' @note 2-  The \code{dicom.set.tag.value} function may take some processing time. 
#' To minimize this time, it is recommended to prepare in advance all the tags to 
#' be modified, and use the \code{dicom.set.tag.value} function only once, as shown in 
#' the example.
#' @return Returns a raw vector, with new tag values.

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
#' @export
dicom.set.tag.value <- function (dicom.raw.data, tag, tag.value, 
                                 tag.dictionary = dicom.tag.dictionary (), ...) {
  
  
  if (length(tag.value)==1) tag.value <- rep (tag.value, length(tag))
  if (length(tag)!=length(tag.value)) {
    warning ("tag and tag.value must have same length or tag.value length must be 1.")
    return (dicom.raw.data)
  }
  
  dicom.df <- NULL
  args <- list(...)
  if (!is.null(args[['dicom.browser']])) {
    dicom.df <- args[['dicom.browser']]
    if (ncol(dicom.df)!= 9) {
      dicom.df <- NULL
    }else {
      nb <- ifelse(is.na(dicom.df$stop[nrow(dicom.df)]), 
                   dicom.df$load.stop[nrow(dicom.df)],dicom.df$stop[nrow(dicom.df)])
      if (is.na(nb)) nb <- length(dicom.raw.data) #pas de vérif dans ce cas
      if (nb != length(dicom.raw.data)) dicom.df <- NULL
    }
  } 
  if (is.null(dicom.df)) dicom.df <- dicom.browser(dicom.raw.data, full.info = TRUE,
                                                   tag.dictionary=tag.dictionary)
  
  
  if (is.null(dicom.df)) {
    warning ("not dicom compliant.")
    return (dicom.raw.data)
  }

  modify.idx <- match(tag, dicom.df$tag)
  
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
  tag.raw.l <- lapply(tag.value_, function(i) raw(0))
  for (conf.idx in 1: length (conformity.flag)) {
    switch(dicom.df$VR[modify.idx[conf.idx]],
           "AE" = {
             new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
             if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
             tag.value_ [conf.idx] <- rawToChar(new.raw)
             tag.raw.l[[conf.idx]] <- new.raw
             m <- !is.na(match(new.raw, as.raw(c(10, 12, 13, 27, 92))))
             
             if (length(new.raw)>16 | any(m)) {
               message (paste ("tag", tag_[conf.idx], "is not AE compliant. This TAG is not modified."))
               conformity.flag[conf.idx] <- FALSE
             } 
           },
           "AS" = {
             new.raw <- charToRaw (tag.value_ [conf.idx])
             m <- match(new.raw[1:3], as.raw(48:57))
             tag.raw.l[[conf.idx]] <- new.raw
             if (length(new.raw)!=4 | any(is.na(m)) | !(new.raw[4] %in% charToRaw ("DWMY"))){
               message (paste ("tag", tag_[conf.idx], "is not AS compliant. This TAG is not modified."))
               conformity.flag[conf.idx] <- FALSE
             }
           },
           "CS" = {
             new.raw <- charToRaw (trimws(tag.value_ [conf.idx]))
             if (length(new.raw)%%2==1) new.raw <- c(new.raw, as.raw(32))
             tag.value_ [conf.idx] <- rawToChar(new.raw)
             tag.raw.l[[conf.idx]] <- new.raw
             m <- match(new.raw, as.raw(c(65:90,48:57,95,32)))
             
             if (length(new.raw)>16 | any(is.na(m))){
               message (paste ("tag", tag_[conf.idx], "is not CS compliant. This TAG is not modified."))
               conformity.flag[conf.idx] <- FALSE
             }
           },
           "DA" = {
             new.raw <- charToRaw (tag.value_ [conf.idx])
             
             m <- match(new.raw[1:8], as.raw(48:57))
             if (length(new.raw)!=8 | any(is.na(m))){
               message (paste ("tag", tag_[conf.idx], "is not DA compliant. This TAG is not modified."))
               conformity.flag[conf.idx] <- FALSE
             }
             tag.raw.l[[conf.idx]] <- new.raw
           },
           "DS" = {
             new.raw <- charToRaw (tag.value_ [conf.idx])
             
             m <- match(new.raw,as.raw(32))
             new.raw <- new.raw[is.na(m)]
             if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
             tag.value_ [conf.idx] <- rawToChar(new.raw)
             tag.raw.l[[conf.idx]] <- new.raw
             
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
             tag.raw.l[[conf.idx]] <- new.raw
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
               tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
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
             tag.raw.l[[conf.idx]] <- new.raw
           }, 
           "UT"= {
             new.raw <- charToRaw (trimws(tag.value_ [conf.idx], which="right"))
             if (length(new.raw)%%2==1) new.raw <- c (new.raw,as.raw(32))
             tag.value_ [conf.idx] <- rawToChar(new.raw)
             tag.raw.l[[conf.idx]] <- new.raw
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
  tag.raw.l<- tag.raw.l[conformity.flag]
  if (length(modify.idx)==0) return (dicom.raw.data)
  
  ###new section (replaces old section)
  
  order.idx <- order(modify.idx, decreasing = FALSE)
  modify.idx <- modify.idx[order.idx]
  tag_ <- tag_[order.idx]
  tag.value_ <- tag.value_[order.idx] 
  tag.raw.l<- tag.raw.l[order.idx]
  
  
  impacted.group <- FALSE
  idx.group  <- grep("0000[)]$", dicom.df$tag)
  if (length(idx.group)>0){
    tag.group <- sapply(dicom.df$tag[idx.group], function (t) substr(t,1,nchar(t)-5))
    tag.reg <-  paste0("^", gsub(")","[)]", gsub("(","[(]",tag.group, fixed = TRUE), fixed = TRUE))
    
    impacted.group <- sapply(tag.reg,function(t) any(grepl(t, dicom.df$tag[modify.idx])))
  }
  

  tag_value_w <- as.numeric(sapply(tag.raw.l, length))
  
  #on vérifie si cela à modifier les charges 
  old.load.width <- dicom.df$stop[modify.idx] + 1 - dicom.df$start[modify.idx]
  to.add <- tag_value_w - old.load.width
  nato.add <- is.na(to.add) # là oùil n'y avait rien avant
  to.add[nato.add] <- tag_value_w[nato.add]
  
  # modify.full.tag <- strsplit(dicom.df$tag[modify.idx[le.ok]],'[ ]')
  modify.full.tag <- strsplit(dicom.df$tag[modify.idx[]],'[ ]')
  modify.full.tag.le <- sapply(modify.full.tag,length)
  
  modify.full.tag <- lapply(1:length(modify.full.tag), function(i){
    t<- data.frame(tag =sapply (1:length(modify.full.tag[[i]]), 
                                function (j) 
                                  paste(modify.full.tag[[i]][1:j], collapse=" ")))
    t$add =to.add[i]
    t
  })
  
  modify.load.tab <- do.call(rbind ,modify.full.tag)
  modify.load.tab$add <- as.numeric(modify.load.tab$add)
  
  # modify.load.tab <- aggregate(modify.load.tab,by =list( as.factor(modify.load.tab$tag)),
  #                              FUN= function(x){ sum(suppressWarnings(as.numeric(x)),na.rm = T)})
  
  L <- by (modify.load.tab, modify.load.tab$tag, FUN= function(x){ sum(x$add,na.rm = T)})
  modify.load.tab <- data.frame(tag = names(L),add=as.numeric(L))
  
  # modify.load.tab$tag <- NULL
  # colnames(modify.load.tab) <- c("tag","add")
  
  modify.load.tab$line <- match(modify.load.tab$tag, dicom.df$tag)
  modify.load.tab$load.start <- dicom.df$load.start[modify.load.tab$line]
  modify.load.tab$load.stop <- dicom.df$load.stop[modify.load.tab$line]  
  modify.load.tab$endian <- dicom.df$endian[modify.load.tab$line] =="little"
  l.load <- lapply(1:nrow(modify.load.tab), function(t) raw())  
  nna.idx <- which(!is.na(modify.load.tab$load.start) & !is.na(modify.load.tab$load.stop))
  l.load[nna.idx] <- lapply(nna.idx, function(i) 
    dicom.raw.data[modify.load.tab$load.start[i]:modify.load.tab$load.stop[i]])
  
  modify.load.tab$old.load <- NA
  dum <- .raw.to.value(l.load[nna.idx], modify.load.tab$endian[nna.idx])
  nm1.f <- dum!=-1
  modify.load.tab$old.load[nna.idx][nm1.f] <- dum[nm1.f]
  
  l.load[nna.idx[nm1.f]] <- .value.to.raw(
    value =modify.load.tab$old.load[nna.idx[nm1.f]] + modify.load.tab$add[nna.idx[nm1.f]],
    raw.length = modify.load.tab$load.stop[nna.idx[nm1.f]]-modify.load.tab$load.start[nna.idx[nm1.f]] + 1,
    little.endian = modify.load.tab$endian[nna.idx[nm1.f]])
  
  raw.idx <- unlist(lapply(nna.idx[nm1.f], function(i) 
    modify.load.tab$load.start[i]:modify.load.tab$load.stop[i]))
  dicom.raw.data[raw.idx] <- unlist(l.load[nna.idx[nm1.f]])
  
  # tableau tag à modifier
  modify.tab <- modify.load.tab[match(dicom.df$tag[modify.idx], modify.load.tab$tag), ]
  modify.tab$start <- dicom.df$start[modify.tab$line]
  modify.tab$stop <- dicom.df$stop[modify.tab$line]
  
  #on s'occupe des groupes
  if (any (impacted.group)){
    idx.group_ <-idx.group[impacted.group]
    group.to.add <-sapply(tag.reg[impacted.group], function (tr) 
      sum(modify.tab$add[grepl(tr, modify.tab$tag)]))
    l.value <- lapply(idx.group_, function(i) dicom.raw.data[dicom.df$start[i]:dicom.df$stop[i]])
    
    l.value <- .value.to.raw (value = .raw.to.value (l.value, dicom.df$endian[idx.group_]=="little") + group.to.add,
                              raw.length= dicom.df$stop[idx.group_]-dicom.df$start[idx.group_] + 1, 
                              dicom.df$endian[idx.group_]=="little")
    
    raw.idx <- unlist(lapply(idx.group_, function(i) dicom.df$start[i]:dicom.df$stop[i]))
    dicom.raw.data[raw.idx] <- unlist(l.value)
  }
  #on s'occupe des tags à modifier
  #--> on découpe d'abord les raw data en liste
  sta1 <- modify.tab$start
  f <- is.na(sta1)
  sta1[f] <- modify.tab$load.stop[f] + 0.5
  sto1 <- sta1 + modify.tab$old.load -1
  sta2 <- modify.tab$stop + 1
  sta2[f] <-  modify.tab$load.stop[f] + 1
  M.cut <- matrix(sort(c(1,modify.tab$load.stop,sta1,sto1,sta2,length(dicom.raw.data))), ncol=2, byrow = TRUE)
  depass <- which(M.cut>length(dicom.raw.data), arr.ind = TRUE)[1]
  if(!is.na(depass)) M.cut <- M.cut[-depass, ]
  rd.l <- lapply(1:nrow(M.cut), function(i) {
    if (M.cut[i,1] %% 1 == 0) return(dicom.raw.data[M.cut[i,1]:M.cut[i,2]])
    return(raw(0))
  })
  rdl.idx <- seq(2,nrow(M.cut),2)
  rd.l[rdl.idx] <- tag.raw.l
  return(unlist(rd.l))
  
}

###old section
  # order.idx <- order(modify.idx, decreasing = TRUE)
  # modify.idx <- modify.idx[order.idx]
  # tag_ <- tag_[order.idx]
  # tag.value_ <- tag.value_[order.idx] 
  # 
  # 
  # for (m.idx in 1:length(modify.idx)) {
  #   idx <- modify.idx[m.idx]
  #   new.raw <- charToRaw(tag.value_[m.idx])
  #   part1 <- NULL
  #   part3 <- NULL
  #   follow.flag <- rep(FALSE, nrow(dicom.df))
  #   
  #   if (dicom.df$load.stop[idx]> 0) part1 <- dicom.raw.data[1:(dicom.df$load.stop[idx])]
  #   if (nrow(dicom.df)>idx) {
  #     part3 <- dicom.raw.data[dicom.df$tag.start[idx+1]: length(dicom.raw.data)]
  #     follow.flag[(idx+1):nrow(dicom.df)] <- TRUE
  #   }
  #   #transformer la partie du dicom.raw.data
  #   dicom.raw.data <- c(part1, new.raw, part3)
  #   
  #   new.raw.load <- packBits(intToBits(length (new.raw)), type="raw")
  #   dicom.raw.data.loadsize <- dicom.df$load.stop[idx]-dicom.df$load.start[idx]+1
  #   
  #   if (dicom.df$endian[idx]=="little") {
  #     dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <- new.raw.load[1:dicom.raw.data.loadsize]
  #   } else {
  #     dicom.raw.data[dicom.df$load.start[idx]:dicom.df$load.stop[idx]] <- rev(new.raw.load[1:dicom.raw.data.loadsize]) 
  #   }
  # 
  #   # MAJ dicom.df
  #   if (!is.na(dicom.df$start[idx])){ # le tag ne contenait rien
  #     to.add <- length(new.raw) - dicom.df$stop[idx] + dicom.df$start[idx] - 1 
  #   } else {
  #     to.add <-  length(new.raw)
  #   }
  #   dicom.df[follow.flag, c("start","stop","load.start","load.stop","tag.start")]  <- 
  #     dicom.df[follow.flag, c("start","stop","load.start","load.stop","tag.start")] + to.add
  #   
  #   dicom.df$stop[idx] <- dicom.df$stop[idx]+ to.add
  #   
  #   # dicom.tag.parser (dicom.df$start[idx], 
  #   #                   dicom.df$stop[idx], dicom.df$VR[idx], 
  #   #                   dicom.df$endian[idx], dicom.raw.data)
  #   
  #   
  #   
  #   if (nchar(dicom.df$encaps.load[idx])>0){
  #   
  #     level.to.change <- suppressWarnings (as.numeric(unlist(strsplit(dicom.df$encaps.load[idx],'[ ]'))))
  #     if (length(level.to.change)>0) {
  #       tag.to.change <- unlist(strsplit(dicom.df$tag[idx],'[ ]'))
  #       encaps <- length(tag.to.change)
  #       tag.to.change <- sapply (1:length(level.to.change), function (i) paste(tag.to.change[1:i], collapse=" "))
  #       tag.to.change.idx <- match(tag.to.change,dicom.df$tag)
  #       tag.to.change.idx  <- tag.to.change.idx[!is.na(level.to.change)]
  #       
  #       for (i in tag.to.change.idx){
  #         size <- dicom.df$load.stop[i]-dicom.df$load.start[i]+1
  #         raw.load <- packBits(intToBits(readBin (dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]], what="int", n= 1,
  #                                                 size = size,
  #                                                 endian = dicom.df$endian[i]) + to.add), type="raw")
  #         
  #         if (size==4){
  #           if (dicom.df$endian[i]=="little") {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-raw.load
  #           } else { dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-rev(raw.load)}
  #         } else {
  #           if (dicom.df$endian[i]=="little") {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-raw.load[1:2]
  #           } else {dicom.raw.data[dicom.df$load.start[i]:dicom.df$load.stop[i]] <-rev(raw.load[1:2]) }
  #         }
  #       }
  #       
  #       # regular.tag.to.change <- rev(paste0("^",gsub("(","[(]",gsub(")","[)]",tag.to.change,fixed = TRUE),fixed = TRUE)))
  #       # first <- idx
  #       # for(encaps.tag in regular.tag.to.change){
  #       #   load.to.change.idx <- grep(encaps.tag, dicom.df$tag)
  #       #   load.to.change.idx <- load.to.change.idx[load.to.change.idx<=first]
  #       #   
  #       #   dicom.df$encaps.load[load.to.change.idx] <- sapply(load.to.change.idx, function (i){
  #       #     u <- as.numeric(unlist(strsplit(dicom.df$encaps.load[i]," ")))
  #       #     if (length(u)==0) return("")
  #       #     f <- c(1:length(u))<encaps
  #       #     paste0(paste(  c(u[f]+ to.add, u[!f]), collapse = " ")," ")
  #       #   })
  #       #   first <- load.to.change.idx[1]-1
  #       #   encaps <- encaps - 1
  #       # }
  #     }
  #   }
  # }
  # 
  # 
  # 
  # idx.group  <- grep("0000[)]$", dicom.df$tag)
  # tag.group <- sapply(dicom.df$tag[idx.group], function (t) substr(t,1,nchar(t)-5))
  # tag.reg <-  paste0("^", gsub(")","[)]", gsub("(","[(]",tag.group, fixed = TRUE), fixed = TRUE))
  # 
  # impacted <- sapply(tag.reg,function(t) any(grepl(t, dicom.df$tag[modify.idx])))
  # if (any(impacted)){
  #   # dicom.df <- dicom.browser(dicom.raw.data,full.info = TRUE, tag.dictionary =tag.dictionary)
  #   
  #   length.group <- sapply(tag.reg[impacted], function(t){
  #     idx <- range(grep(t,dicom.df$tag))+1
  #     e <- ifelse(idx[2]<= nrow(dicom.df),dicom.df$tag.start[idx[2]], length(dicom.raw.data)+1)
  #     e-dicom.df$tag.start[idx[1]]})
  #   
  #   idx.group_ <-idx.group[impacted]
  #   for (i in 1:length(idx.group_)){
  #     size <- dicom.df$stop[idx.group_[i]]-dicom.df$start[idx.group_[i]]+1
  #     raw.load_ <- packBits(intToBits(length.group[i]), type="raw")
  #     
  #     if(!is.na(size)){
  #       if (size==4){
  #         if (dicom.df$endian[i]=="little") {raw.load <- raw.load_
  #         } else {raw.load <-rev(raw.load_)}
  #       } else {
  #         if (dicom.df$endian[i]=="little") {raw.load<- raw.load_[1:2]
  #         } else { raw.load<-rev(raw.load_[1:2])}
  #       }
  #       
  #       dicom.raw.data[dicom.df$start[idx.group_[i]]:dicom.df$stop[idx.group_[i]]] <-raw.load
  #       
  #     }
  #   }
  # }
  # 
  # return (dicom.raw.data)
# }


.raw.to.value <- function(raw.list, little.endian ){
  raw.list.le <- as.numeric(sapply(raw.list,length))
  le4.f <- raw.list.le==4
  le2.f <- raw.list.le==2
  # little.endian <- db$endian=="little"
  
  load.v <- rep(0, length(raw.list))
  f <- le2.f & little.endian
  if (any(f)) load.v [f] <- readBin(as.raw(unlist (raw.list[f])), what = "int", n = sum(f),
                                    size = 2, signed = FALSE, endian = "little")
  f <- le4.f & little.endian
  if (any(f)) load.v [f] <- readBin(as.raw(unlist (raw.list[f])), what = "int", n = sum(f),
                                    size = 4,  endian = "little")
  f <- le2.f & !little.endian
  if (any(f)) load.v [f] <- readBin(as.raw(unlist (raw.list[f])), what = "int", n = sum(f),
                                    size = 2, signed = FALSE, endian = "big")
  f <- le4.f & !little.endian
  if (any(f)) load.v [f] <- readBin(as.raw(unlist (raw.list[f])), what = "int", n = sum(f),
                                    size = 4,  endian = "big")
  return (load.v)
}

.value.to.raw <- function(value, raw.length, little.endian){
  le4.f <- raw.length==4
  le2.f <- raw.length==2
  
  l <- lapply(1:length(value), function(i) raw(0))
  f <- le2.f & little.endian
  if (any(f)) l [f]  <- lapply(value[f],function(v) packBits(intToBits(v), type="raw")[1:2])
  
  f <- le4.f & little.endian
  if (any(f)) l [f]  <- lapply(value[f],function(v) packBits(intToBits(v), type="raw"))
  
  f <- le2.f & !little.endian
  if (any(f)) l [f]  <- lapply(value[f],function(v) packBits(intToBits(v), type="raw")[2:1])
  
  f <- le4.f & !little.endian
  if (any(f)) l [f]  <- lapply(value[f],function(v) packBits(intToBits(v), type="raw")[4:1])
  
  return(l)
  
}