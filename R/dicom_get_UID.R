
.dicom.get.UID <- function (dicom.raw.data, white.list = c("instance","reference"), 
                            black.list = c("frame of reference","class"),
                            tag.dictionary = dicom.tag.dictionary (),...){
  
  
  dicom.df <- NULL
  args <- list(...)
  if (!is.null(args[['dicom.browser']])) {
    dicom.df <- args[['dicom.browser']]
    nb <- ifelse(is.na(dicom.df$stop[nrow(dicom.df)]), 
                 dicom.df$load.stop[nrow(dicom.df)],dicom.df$stop[nrow(dicom.df)])
    if (is.na(nb)) nb <- length(dicom.raw.data) #pas de vérif dans ce cas
    if (nb != length(dicom.raw.data)) dicom.df <- NULL
  } 
  if (is.null(dicom.df)) dicom.df <- dicom.browser(dicom.raw.data, full.info = TRUE,
                                                   tag.dictionary=tag.dictionary)
  
  last.tag <-sapply(dicom.df$tag, function(t) rev(unlist(strsplit(t, "[ ]")))[1])
  m <- match(last.tag, tag.dictionary$tag)
  name<- tag.dictionary$name[m]
  vr <- tag.dictionary$VR[m]
  
  white.list = tolower(white.list)
  white.list <- gsub("[ ]","[ ]", tolower(white.list))
  white.list <- gsub("[[ ]]","[ ]", white.list, fixed = TRUE)
  
  black.list = tolower(black.list)
  black.list <- gsub("[ ]","[ ]", tolower(black.list))
  black.list <- gsub("[[ ]]","[ ]", black.list, fixed = TRUE)
  
  vrUI <- vr=="UI"
  vrUI[is.na(vrUI)] <- FALSE
  
  has.load <- !is.na(dicom.df$start) & !is.na(dicom.df$stop)
  keep <- apply(do.call(rbind,lapply(white.list,function(l) grepl(l,tolower(name)))),
                2,sum)>0
  ui.idx <- which(apply(do.call(rbind,c(lapply(black.list,function(l) !grepl(l,tolower(name))), 
                                        list(keep,vrUI,has.load))),2,prod)==1)
  
  # ui.idx <- which(vr=="UI" & !grepl("[ ]class", tolower(name)) &  
  #                   !grepl("frame of reference", tolower(name)))
  
  Value <- sapply(ui.idx,
                  function(i) dicom.tag.parser (dicom.df$start[i], 
                                                dicom.df$stop[i], dicom.df$VR[i], 
                                                dicom.df$endian[i], dicom.raw.data))
  ui.idx.f <- !grepl("^1[.]2[.]840[.]10008",Value)
  ui.idx <- ui.idx[ui.idx.f]
  Value <- Value[ui.idx.f]
  
  tab<- data.frame(line = ui.idx, tag= dicom.df$tag[ui.idx], name =name[ui.idx], 
                   VR =vr[ui.idx], Value = Value)
  rownames(tab) <- NULL
  
  
  return(tab)
}


.Rdcm.get.UID <- function (Rdcm.files, white.list = c("SOP","reference"), 
                           black.list = c("frame of reference","class"),
                           tag.dictionary = dicom.tag.dictionary ()){
  if (length(Rdcm.files)==0) {
    warning ("no files to analyze.")
    return (NULL)
  } 
  flag <- dir.exists(Rdcm.files)
  Rdcm.dir <- Rdcm.files[flag]
  Rdcm.filenames1 <-  list.files(Rdcm.dir,recursive = TRUE,full.names = TRUE)
  Rdcm.filenames2 <- Rdcm.files[!flag]
  Rdcm.filenames2 <- Rdcm.filenames2[file.exists(Rdcm.filenames2)]
  Rdcm.filenames <- c(Rdcm.filenames1,Rdcm.filenames2)
  if (length(Rdcm.filenames)==0) {
    warning ("no files to analyze.")
    return (NULL)
  } 
  
  Rdcm.filenames <-Rdcm.filenames[grepl("[.]Rdcm$",Rdcm.filenames)]
  if (length(Rdcm.filenames)==0) {
    warning ("no files to analyze.")
    return (NULL)
  } 
  
  tab <- do.call (rbind, lapply(Rdcm.filenames,function(object.name) {
    obj <- load.Rdcm.raw.data(object.name)
    if (length(obj$address) == 0) return(NULL)
    dicom.fs <- obj$header$object.info$dicom.file
    if (length(dicom.fs) == 0) dicom.fs <- rep("", length(obj$address))
    tab <- do.call(rbind, lapply (1:length(obj$address), function(i.df) {
      dicom.df <- obj$address[[i.df]]
      last.tag <-sapply(dicom.df$tag, function(t) rev(unlist(strsplit(t, "[ ]")))[1])
      m <- match(last.tag, tag.dictionary$tag)
      name<- tag.dictionary$name[m]
      vr <- tag.dictionary$VR[m]
      
      white.list = tolower(white.list)
      white.list <- gsub("[ ]","[ ]", tolower(white.list))
      white.list <- gsub("[[ ]]","[ ]", white.list, fixed = TRUE)
      
      black.list = tolower(black.list)
      black.list <- gsub("[ ]","[ ]", tolower(black.list))
      black.list <- gsub("[[ ]]","[ ]", black.list, fixed = TRUE)
      
      vrUI <- vr=="UI"
      vrUI[is.na(vrUI)] <- FALSE
      
      keep <- apply(do.call(rbind,lapply(white.list,function(l) grepl(l,tolower(name)))),
                    2,sum)>0
      ui.idx <- which(apply(do.call(rbind,c(lapply(black.list,function(l) !grepl(l,tolower(name))), 
                                            list(keep,vrUI))),2,prod)==1)
      
      
      Value  <- as.character(unlist(obj$data[[i.df]][ui.idx]))
      ui.idx.f <- !grepl("^1[.]2[.]840[.]10008",Value)
      ui.idx <- ui.idx[ui.idx.f]
      Value <- Value[ui.idx.f]
      df <- data.frame(line = ui.idx, tag= dicom.df$tag[ui.idx], name =name[ui.idx], 
                       VR =vr[ui.idx], Value = Value)
      df$dicom.file <- dicom.fs[i.df]
      df
    }))}))
  rownames(tab) <- NULL
  return(tab)
}
