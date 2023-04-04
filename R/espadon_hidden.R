.onUnload <- function (libpath) {
  library.dynam.unload("espadon", libpath)
}

# .onAttach <- function(libname,pkgname) {
#   packageStartupMessage(c('Make sure you always have the latest version of espadon.\n',
#                           'For checking, copy-paste :\n',
#                           '   db <- available.packages(filters = c("duplicates"))\n',
#                           '   if (packageVersion ("espadon") < as.character (db["espadon","Version"])){\n',
#                           '       message ("you need to update espadon")\n',
#                           '   } else {message ("espadon is up to date")}\n '))
# }

#####################################################################################
.espadon.version <- function(){
  return ("1.3.2")
}

#####################################################################################

.espadon.UID <- function(){
  return ("1.2.826.0.1.3680043.10.1002")
}
######################################################
#' @import progress
#' @import qs

.create.Rdcm.object <- function (dcm.filenames,existing.obj = NULL, verbose = TRUE,
                                 update = FALSE, save.flag = TRUE, save.dir =NULL, 
                                 only.header = FALSE, ignore.duplicates = FALSE,
                                 tag.dictionary = dicom.tag.dictionary ()){
  
  tempqs.dir <- file.path (tempdir(), "tempqs.dir") 
  if (!dir.exists(tempqs.dir)) dir.create (tempqs.dir, recursive = TRUE) 
  on.exit(
    expr = {
      if (dir.exists(tempqs.dir)) unlink (tempqs.dir, recursive = TRUE)
    })

  colnames.df <- c("reference", 
                   "acquisition.date", "study.date","study.time","creation.date", 
                   "modality", 
                   "SOP.ID","SOP.UID","SOP.type",
                   "transfer.syntax.UID", "implementation.ID",
                   "scanning.sequence", "study.description", "serie.description", 
                   "study.ID", "study.UID","serie.UID", 
                   "PID", "birthday","sex",
                   "slice.nb",
                   "file.idx", 
                   "ref.label", ##
                   "SOP.label.nb","error")
  tag <- c("[(]0020,0052[)]$", 
           "^[(]0008,0021[)]$|^[(]0008,0022[)]$|^[(]0008,0023[)]$|^[(]3006,0008[)]$|^[(]300A,0006[)]$", 
           "^[(]0008,0020[)]$", 
           "^[(]0008,0030[)]$", "^[(]0008,0012[)]$", 
           "^[(]0008,0060[)]$",
           "^[(]0008,0016[)]$","^[(]0008,0018[)]$","^[(]0008,0008[)]$",
           "^[(]0002,0010[)]$", "^[(]0002,0012[)]$", 
           "^[(]0018,0020[)]$", "^[(]0008,1030[)]$", "^[(]0008,103E[)]$", 
           "^[(]0020,0010[)]$", "^[(]0020,000D[)]$", "^[(]0020,000E[)]$",
           "^[(]0010,0020[)]$", "^[(]0010,0030[)]$","^[(]0010,0040[)]$", 
           "^[(]0020,0013[)]$")
  
  if (!is.null(existing.obj)){
    existing.obj[is.na(existing.obj)]<- ""
    old.object.df_ <-   castlow.str (apply (existing.obj[,c("study.date","study.time","creation.date", 
                                                            "modality", 
                                                            "SOP.ID", "SOP.type",
                                                            "transfer.syntax.UID", "implementation.ID",
                                                            "scanning.sequence", "study.description", "serie.description", 
                                                            "study.ID", "study.UID","serie.UID")],1,paste, collapse=""))
  } else old.object.df_<- NULL
  
  
  df <- data.frame (matrix(data=rep("", length(dcm.filenames)*length(colnames.df)), 
                           ncol = length(colnames.df)))
  colnames(df)<-colnames.df
  df$error <- FALSE
  out.idx <- 1
  
  # lapply (1:length(in.file.names), function (f.idx){
  #   m <- .dicom.load.raw.data(in.file.names[f.idx])
  #   .get.dcm.info(m, tag=tag, tag_stop="(0020,0052)")
  # })
  if (verbose) pb <- progress_bar$new(format = " downloading [:bar] :percent",
                                      total = length(dcm.filenames), width= 60)
  
  
 
  pb.idx <- 0
  for (f.idx in 1:length(dcm.filenames)){
    m <- dicom.raw.data.loader(dcm.filenames[f.idx])
    dicom.df <- dicom.browser (m,tag.dictionary = tag.dictionary)
    
    if (is.null(dicom.df)) {
      L <- NULL
    } else {
      # L <- list (data = lapply(1:nrow(dicom.df), function (idx) 
      #   dicom.tag.parser (dicom.df$start[idx], dicom.df$stop[idx], 
      #                     dicom.df$VR[idx], dicom.df$endian[idx],
      #                     m, try.parse= FALSE)),
      #   address=dicom.df, filename = dcm.filenames[f.idx])
      L <- list (data = lapply(1:nrow(dicom.df), function (idx){
        if ((is.na(dicom.df$start[idx]) | is.na(dicom.df$stop[idx])) & dicom.df[idx,"tag"]!="(7FE0,0010)") return (NA)
        if (dicom.df[idx,"tag"]=="(7FE0,0010)"){
          if (is.na(dicom.df$start[idx]) | is.na(dicom.df$stop[idx])) return(raw(0))
              return(m [dicom.df$start[idx]:dicom.df$stop[idx]])
        } 
        dicom.tag.parser (dicom.df$start[idx], dicom.df$stop[idx], 
                          dicom.df$VR[idx], dicom.df$endian[idx],
                          m, try.parse= FALSE)}),
        address=dicom.df, filename = dcm.filenames[f.idx])
      names(L$data) <- dicom.df$tag
    }
    
    if (!is.null(L)) {
      n.L <- names (L$data)
      df[f.idx,1:length(tag)] <- sapply(tag, function (t) {
        value<-sapply(grep(t,n.L),function(i) L$data[[i]])
        if (length(value)==0) return ("")
        value <- value [!is.na(value)]
        if (length(value)==0) return ("")
        value <- value [ trimws(value)!=""]
        if (length(value)==0) return ("")
        return (sort(value)[1])
      })
      
      if (!is.null(L$data[["(7FE0,0010)"]]) & length(L$data[["(7FE0,0010)"]])==0) df[f.idx,c("error")] <- TRUE
      
      if (!is.null(old.object.df_)) {
        ind <- which( old.object.df_==
                        castlow.str (apply (df[f.idx,c("study.date","study.time","creation.date", 
                                                       "modality", 
                                                       "SOP.ID", "SOP.type",
                                                       "transfer.syntax.UID", "implementation.ID",
                                                       "scanning.sequence", "study.description", "serie.description", 
                                                       "study.ID", "study.UID","serie.UID")],1,paste, collapse="")) )
        if (length(ind)!=0) df[f.idx,c("ref.label", "SOP.label.nb")] <- existing.obj[ind,c("ref.label",  "SOP.label.nb")]
      }
      
      if ((update) | (df$SOP.label.nb[f.idx]=="" & !update)){
        tf <- tempfile(paste("qs",out.idx,"_",sep=""), tmpdir = tempqs.dir)
        df[f.idx, ]$file.idx <- basename(tf)
        qsave(L, file = tf)
        out.idx <- out.idx+1
      }
    }
    if (verbose) pb$tick()
    
  }
  
  df$cd <-  format (file.info(dcm.filenames)$mtime,format = "%Y%m%d")
  
  df <- df[df$SOP.ID !="",]
  if (nrow(df)==0) return (NULL)
  df <- df[order(df$study.date,df$acquisition.date,df$PID,df$modality, df$SOP.ID, 
                 df$implementation.ID, df$SOP.type, 
                 df$study.ID, df$study.UID, as.numeric(df$slice.nb)), ]
  
  # object.df <- unique (df[ ,c(1:17)])
  object.df <- unique (df[ ,c("reference", 
                              "acquisition.date", "study.date","study.time","creation.date", 
                              "modality", 
                              "SOP.ID","SOP.type","transfer.syntax.UID", "implementation.ID",
                              "scanning.sequence", "study.description", "serie.description", 
                              "study.ID", "study.UID","serie.UID",
                              "ref.label","SOP.label.nb")])
  object.df$PID <-  df$PID[match(row.names(object.df), row.names(df))]
  object.df$outfilename= ""
  object.df$input_file_nb = 1
  object.df$cd <- df$cd[match(row.names(object.df), row.names(df))]
  row.names(object.df) <- NULL
  
  
  if (any(object.df$ref.label=="")){ #| any(object.df$SOP.label.nb=="")){
    #on cherche le label des referentiel, en utilisant ceux qui existent deja dans le old
    ref.tab <- unique(object.df[,c("ref.label","reference")])
    ref.tab <- ref.tab[order(ref.tab$ref.label, decreasing=TRUE),]
    ref.tab <- ref.tab[!duplicated(ref.tab$reference),]
    ref.tab <- ref.tab[order(ref.tab$ref.label),]
    ref.possible <- "1": as.character(nrow (ref.tab))
    ref.possible <- ref.possible[is.na(match(ref.possible,ref.tab$ref.label [ ref.tab$ref.label!=""]))]
    
    ref.name <- table (object.df$reference)
    ref.name <- names (ref.name) [order(ref.name,decreasing =TRUE)]
    ref.tab <- ref.tab[match(ref.tab$reference,ref.name), ]
    vide <- which(ref.tab$ref.label=="")
    if (length(vide)>0){
      ref.tab[vide,"ref.label"]<-ref.possible[1:length(vide)]
    }
    object.df$ref.label <- ref.tab$ref.label [match (object.df$reference, ref.tab$reference)]
  }
  object.df$temp <- ""
  object.df$SOP.UID <- ""
  object.df$OK <- TRUE
  paste.df <- castlow.str(apply(df[,c("acquisition.date", "study.date","study.time","creation.date", 
                                      "modality","SOP.ID","SOP.type", 
                                      "transfer.syntax.UID", "implementation.ID",
                                      "scanning.sequence", 
                                      "study.ID", "study.UID","serie.UID")], 1,paste,collapse=""))
  paste.object.df <- castlow.str(apply(object.df[,c("acquisition.date", "study.date","study.time","creation.date", 
                                                    "modality", "SOP.ID", "SOP.type",
                                                    "transfer.syntax.UID","implementation.ID",
                                                    "scanning.sequence", 
                                                    "study.ID", "study.UID","serie.UID")], 1,paste,collapse=""))
  #on met à jour le nombre de files, et le début du nom de sauvegarde de l'objet
  for (obj.idx in 1:nrow(object.df)){
    fl <- paste.df==paste.object.df[obj.idx]
    obj.info <- df[fl,]
    obj.info$slice.nb <- as.numeric(obj.info$slice.nb)
    obj.info <- obj.info[order (obj.info$slice.nb),]
    object.df$temp[obj.idx] <- paste (unique(obj.info$file.idx [!is.na(obj.info$file.idx)]), collapse=";")
    object.df$SOP.UID[obj.idx] <- paste (sort(unique(obj.info$SOP.UID)), collapse=";")
    object.df$OK[obj.idx] <- !any(obj.info$error)
    study.date <- unique(obj.info$study.date)
    study.date <- study.date[study.date!=""]
    study.date <- ifelse(length(study.date)==0, "", min(study.date))
    
    creation.date  <- unique(obj.info$creation.date)
    creation.date <- creation.date[creation.date!=""]
    creation.date <- ifelse(length(creation.date)==0, "", min(creation.date))
    
    cd <-  min(unique(object.df$cd[obj.idx]))
    date <- ifelse(study.date=="",ifelse (creation.date=="",cd,creation.date) ,study.date)
    object.df$input_file_nb[obj.idx]  <- nrow(obj.info)
    object.df$outfilename[obj.idx] <- paste (date, "_", trimws(object.df$PID[obj.idx]),"_ref",object.df$ref.label[obj.idx],"_do", sep="")
    
  }
  
  #on n'a plus besoin de $cd
  object.df$cd <- NULL
   #on enlève les objet qui utilisent sur les meme fichiers temporaires
  object.tp <- sapply (object.df$temp,function(t) v <- paste(sort(unlist(strsplit(t,";"))), collapse=";"))
  object.df <- object.df[!(duplicated(object.tp)) | object.tp=="", ]
  
  #on trie par rapport en commençant par les anciens objet
  object.df <- object.df[order (object.df$SOP.label.nb, object.df$OK,object.df$acquisition.date,decreasing = TRUE), ]
  
  
  
  flag.duplicated.SOP <- duplicated(object.df$SOP.UID) 
  flag.duplicated.to.keep <- flag.duplicated.SOP & object.df$SOP.label.nb!="" & object.df$temp!="" 
  flag.duplicated.to.remove <- flag.duplicated.SOP & object.df$SOP.label.nb=="" & object.df$temp!="" 
  if (any(flag.duplicated.SOP)){
    
    if (any(flag.duplicated.to.keep))  warning("Some old Rdcm objects are duplicated and have been updated", call. = FALSE)
    if (any(flag.duplicated.to.remove)){
      if (ignore.duplicates){
        tempf <- file.path (tempqs.dir,  unlist(strsplit(object.df[flag.duplicated.to.remove, ]$temp,";")))
        object.df <- object.df[!flag.duplicated.to.remove,]
        file.remove(tempf, recursive = TRUE)
      } else {
        warning("Some objects are duplicated. Consider ignore.duplicates = TRUE", call. = FALSE)
      }
    }
  }
  #on trie par rapport au referentiel, puis par date, et de la modalité
  object.df <- object.df[order (object.df$ref.label, object.df$outfilename,
                                object.df$acquisition.date, tolower(object.df$modality)), ]
  
  
  #pour chaque ref, on met à jour le label de l'objet
  if ( any(object.df$SOP.label.nb=="")){
    for (ref.idx in 1:nrow(ref.tab)){
      do <-which(object.df$reference==ref.tab$reference [ref.idx])
      SOP.label.nb <- object.df$SOP.label.nb[do]
      label.possible <- "1": as.character(length (SOP.label.nb))
      vide <- which(SOP.label.nb=="")
      if (length(vide)> 0){
        label.possible <- label.possible[is.na(match(label.possible,SOP.label.nb[-vide]))]
        object.df$SOP.label.nb[do[vide]] <- label.possible[1:length(vide)]
      }
    }
  }
  
  #met à jour les noms  de sauvegarde de l'objet
  object.df$outfilename <- sapply (1:nrow(object.df), function (obj.idx){
    paste(object.df$outfilename[obj.idx], object.df$SOP.label.nb[obj.idx], "_",
          gsub("[ ]", "",tolower(object.df$modality[obj.idx])), sep="")
  })
  
  row.names(object.df) <- NULL
  L.list <- NULL
  
  object.df <- object.df[object.df$temp!="", ]
  if (nrow(object.df)==0) {
    warning("No new objects to convert.", call. = FALSE)
    return(NULL)
    }
  if (verbose) {
    pb <- progress_bar$new(format = " objects creation [:bar] :percent",
                           total = nrow(object.df), width= 60)
  }
  for (obj.idx in (1:nrow(object.df))) {
    modality <- castlow.str(object.df[obj.idx,]$modality)
    
    tempf <- file.path (tempqs.dir,  unlist(strsplit(object.df[obj.idx, ]$temp,";")))
    
    data <- lapply(tempf,function(f) qread (f))
    L <-.obj.save.by.modality (modality ,object.df[obj.idx, ], data, only.header, save.flag)
    
    if (is.null(L)) {
      fn <- paste(sapply(data,function(d) basename(d$filename)),collapse = ", ")
      if (nchar(fn)>200)  fn <- paste0(substr(fn,1,200),"...")
      dn <- paste(unique(sapply(data,function(d) basename(dirname(d$filename)))),collapse = ", ")
      warning (paste(fn,"in",dn, "directory can not be read by espadon. Please, contact maintainer."),
               call. = FALSE)
    } else {
      if (save.flag){
        for (L.idx in  1:length(L))  {
          .save.dicom.raw.data.to.Rdcm (L[[L.idx]], 
                                        file.path (save.dir, L[[L.idx]]$header$file.basename))
          L.list <- c (L.list, L[[L.idx]]$header$file.basename)
        }
      } else {
        L.list <- do.call(c, list(L.list, L))
      }
    }

    file.remove(tempf, recursive = TRUE)
    if (verbose ) pb$tick()
    
    
  }
  return (L.list)
}


################################################################################
.obj.save.by.modality <- function(modality ,object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  L <- tryCatch(switch (modality,
                        "rtstruct" = {.rtstruct.save (object.info, data, only.header, Rdcm.mode)},
                        "reg" = {.reg.save (object.info, data, only.header, Rdcm.mode)},
                        "rtdose" = {.rtdose.save (object.info, data, only.header, Rdcm.mode)},
                        "ct" = {.img3Dplan.save (object.info, data, only.header, Rdcm.mode)},
                        "ct1d" = {.img3Dplan.save (object.info, data, only.header, Rdcm.mode)},
                        "mr" = {.img3Dplan.save (object.info, data, only.header, Rdcm.mode)},
                        "pt" = {.img3Dplan.save (object.info, data, only.header, Rdcm.mode)},
                        "rtplan" = {.rtplan.save (object.info, data, only.header, Rdcm.mode)},
                        "rtimage" = {.rtimage.save (object.info, data, only.header, Rdcm.mode)},
                        
                        .other.save (object.info, data, only.header,Rdcm.mode)
  ), error = function(e) NULL)
  return(L)
}
################################################################################
.save.dicom.raw.data.to.Rdcm <- function (obj, Rdcm.filename) {
  h <- qserialize(c(obj$header,list(espadon.version=.espadon.version())))
  a <- numeric(0)
  if (!is.null(obj$address)) a <- qserialize(obj$address)
  d <- qserialize(obj$data)
  zz <-  file(Rdcm.filename, "wb")
  writeBin(length(h),zz,size=4, endian="little")
  writeBin(length(a),zz,size=4, endian="little")
  writeBin(length(d),zz,size=4, endian="little")
  writeBin(h,zz,endian="little")
  writeBin(a,zz,endian="little")
  writeBin(d,zz,endian="little")
  close(zz)
}
#############################################################################################
.rtimage.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  
  address.all <- lapply(data, function(l) l$address)
  data.all <- lapply(data, function(l) l$data)
  filename.all <- sapply (data, function(l) l$filename)
  nb.of.rtimg <- length(address.all)
  L <- list ()
  for (rtimg.idx in 1:length(address.all)) {
    address <- address.all[[rtimg.idx]]
    data <- data.all[[rtimg.idx]]
    filename <- filename.all [rtimg.idx]
    
    name <- names(data)
    header <- list()
    
    
    header$patient <- trimws (tryCatch (data[[grep("^[(]0010,0020[)]$",name)]],error = function (e) ""))
    header$patient.name <- tryCatch (data[[grep("^[(]0010,0010[)]$",name)]],error = function (e) "")
    if (is.na(header$patient.name)) header$patient.name <- ""
    header$patient.name <- trimws(header$patient.name)                               
    header$patient.bd <- tryCatch (data[[grep("^[(]0010,0030[)]$",name)]],error = function (e) "")
    header$patient.sex <- toupper(trimws (tryCatch (data[[grep("^[(]0010,0040[)]$",name)]],error = function (e) "")))
    header$patient.position <- toupper(trimws (tryCatch (data[[grep("^[(]0018,5100[)]$",name)]],error = function (e) "")))
    
    header$file.basename <- ""
    header$file.dirname <- ""
    header$object.name <- object.info$outfilename
    header$object.alias <- ""
    header$ref.object.alias <- ""
    
    header$object.info <- list ()
    header$ref.object.info <- list ()
    
    
    post.tag <- name[grep ("^[(]300C,0002[)] item[[:digit:]]+$", name)]
    db.name <- c( "SOP.ID", "SOP.label")
    dg.tag <- c("(0008,1150)", "(0008,1155)")
    num.tag <- numeric(0)
    reference.info <- .get.tag.info (db.name = db.name, 
                                     db.tag = dg.tag, 
                                     num.tag = num.tag, post.tag = post.tag,
                                     data = data, data.name = name)
    
    if (is.null(reference.info)){
      header$ref.object.info <- NULL
      header$ref.object.alias <- NULL
    } else {
      header$ref.object.info$SOP.ID <- as.character(unique (reference.info$SOP.ID))
      header$ref.object.info$SOP.label <- as.character(unique (reference.info$SOP.label))
    }
    
    
    
    header$object.info$SOP.ID <- object.info$SOP.ID
    header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
    header$object.info$implementation.ID <- object.info$implementation.ID
    header$object.info$SOP.type <- object.info$SOP.type
    header$object.info$study.ID <- object.info$study.ID
    header$object.info$study.UID <- object.info$study.UID
    header$object.info$serie.UID <- object.info$serie.UID
    header$object.info$scanning.sequence <- object.info$scanning.sequence
    header$object.info$SOP.label <- data[[grep("^[(]0008,0018[)]$",name)]]
    header$object.info$encoding <- tryCatch (data[[grep("^[(]0008,0005[)]$",name)]],error = function (e) "")
    if (Rdcm.mode) {
      header$object.info$dicom.file  <- sort(basename(filename))
    }
    
    header$object.info$nb.of.subobject <- nb.of.rtimg
    
    
    header$frame.of.reference <-  object.info$reference
    header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
    
    header$modality <- castlow.str (object.info$modality)
    header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
    
    
    # header$acq.date <- tryCatch (data[[grep("^[(]0008,0023[)]$",name)]],error = function (e) "")
    header$acq.date <- object.info$acquisition.date
    header$study.date <- tryCatch (data[[grep("^[(]0008,0020[)]$",name)]],error = function (e) "")
    header$creation.date <- tryCatch (data[[grep("^[(]0008,0012[)]$",name)]],error = function (e) "")
    header$study.time <- tryCatch (data[[grep("^[(]0008,0030[)]$",name)]],error = function (e) "")
    header$approval.status <-  tryCatch (data[[grep("^[(]300E,0002[)]$",name)]],error = function (e) "")
    
    header$number <- as.numeric(tryCatch (data[[grep("^[(]0020,0013[)]$",name)]],error = function (e) 1))
    header$object.alias <- paste(object.info$outfilename, header$number, sep="")
    if (Rdcm.mode) {
      header$file.basename <- paste(object.info$outfilename, header$number, ".Rdcm", sep="")
    } else {
      header$file.basename <- sort(basename(filename))
      header$file.dirname <- unique(dirname(filename))
    }
    
    info.name <- c ("conversion.type","exposure.time","exposure.duration",
                    "photometric.interpretation",
                    "label", "name","description",                                
                    "machine.name")
    info.tag <- c("(0008,0064)", "(0018,1150)",  "(0018,9328)", 
                  "(0028,0004)",
                  "(3002,0002)", "(3002,0003)", "(3002,0004)", 
                  "(3002,0020)")
    info.num.tag <- c(2, 4)
    pixel.name <- c ("sample.per.pixel","planar.configuration", "n.i", "n.j", 
                     "pixel.ratio", "bit.allocated","bit.stored","MSB",
                     "signed", "pixel.relationship","sign.relation.ship","pixel.spacing",
                     "min.value", "max.value", "intercept","slope",
                     "unit")
    pixel.tag <- c("(0028,0002)", "(0028,0006)", "(0028,0010)", "(0028,0011)",
                   "(0028,0034)", "(0028,0100)", "(0028,0101)", "(0028,0102)", 
                   "(0028,0103)", "(0028,1040)", "(0028,1041)","(3002,0011)",
                   "(0028,0106)", "(0028,0107)","(0028,1052)","(0028,1053)",
                   "(0028,1054)")
    pixel.num.tag <- c(1:9,11,13:16)
    beam.ref.name <- c ("diaphragm.position","collimator.angle")
    beam.ref.tag <- c("(3002,0034)", "(300A,0120")
    beam.ref.numtag <- c(2)
    patient.ref.name <-  c("machine.isocenter.position")
    patient.ref.tag <- c("(300A,012C)")
    patient.ref.numtag <-integer(0)
    fixed.ref.name <- c("gantry.angle", "patient.support.angle","gantry.pitch.angle")
    fixed.ref.tag <- c("(300A,011E)", "(300A,0122)", "(300A,014A)")
    fixed.ref.numtag <- 1:3
    table.ref.name <- c ("table.to.eccentric.axis.distance","eccentric.angle", 
                         "vertical.position", "longitudinal.position",
                         "lateral.position", "pitch.angle","roll.angle")
    table.ref.tag <- c("(300A,0124)", "(300A,0125)", 
                       "(300A,0128)", "(300A,0129)",
                       "(300A,012A)", "(300A,0140)", "(300A,0144)")
    table.numtag <- c(1:7)
    exposure.info.name <- c("referenced.frame.number", "KVP", "exposureTime", "XRay.tube.current",
                            "exposure.time.ms","XRay.tube.current.mA", "meterset.exposure")
    exposure.info.ref.tag <- c("(0008,1160)", "(0018,0060)", "(0018,1150)", "(0018,1151)",
                               "(0018,9328)","(0018,9330)","(3002,0032)")
    exposure.info.numtag  <- c(1:7)
    image.ref.name <- c ("orientation.to.beam","orientation","upper.left.position" )
    image.ref.tag <- c("(3002,000C)","(3002,0010)", "(3002,0012)")
    image.ref.num.tag <- numeric(0) #c()
    coll.name <- c("type", "source.distance", "nb","position.boundaries","position")
    coll.tag <- c("(300A,00B8)", "(300A,00BA)", "(300A,00BC)", "(300A,00BE)", "(300A,011C)")
    coll.numtag <- c(2:4)
    distance.name <- c ("source.gantry", "source.patient", "source.image","source.referenceObject")
    distance.tag <-  c("(3002,0022)", "(3002,0024)","(3002,0026)", "(3002,0028)")
    distance.num.tag <- 1:4
    gantry.ref.name <- c ("img.origin.position", "image.angle", "collimator.angle")
    gantry.ref.tag <- c( "(3002,000D)","(3002,000E)","(300A,0120)")
    gantry.ref.num.tag <- c(2:3)
    
    pixel.info <- .get.tag.info (pixel.name, pixel.tag, pixel.num.tag, data = data, data.name = name)
    
    image.ref <- .get.tag.info (db.name = image.ref.name, 
                                db.tag = image.ref.tag, 
                                num.tag = image.ref.num.tag, data = data, data.name = name)
    distance <- .get.tag.info (distance.name, distance.tag, 
                               distance.num.tag, data = data, 
                               data.name = name)
    
    header$unit <- pixel.info$unit
    
    
    
    if (is.null(pixel.info$n.i)) header$error <- c(header$error,"(0028,0011) number of column error")
    if (is.null(pixel.info$n.j)) header$error <- c(header$error,"(0028,0010) number of line error")
    header$n.ijk <- c (pixel.info$n.i[1], pixel.info$n.j[1], 1)
    
    ## slice.thickness 
    header$slice.thickness <-  0
    
    header$min.pixel <- NA
    header$max.pixel <- NA
    
    ##espacement des pixels en x, y
    if (is.null(pixel.info$pixel.spacing)) header$error <- c (header$error,"(3002,0011) xy-spacing error")
    # header$dxyz  <- c(as.numeric(unlist(strsplit(pixel.info$pixel.spacing[1],"\\\\"))),1e-6)
    # header$slice.thickness <- header$dxyz[3] <- min(header$dxyz[1:2]/1000)
    header$dxyz  <- c(as.numeric(unlist(strsplit(pixel.info$pixel.spacing[1],"\\\\"))),0)
    header$slice.thickness <- header$dxyz[3] <- 0
    
    if (image.ref$orientation.to.beam == "NORMAL") {
      idx.orientation <- which(colnames(image.ref) == "orientation")
      
      if (length(idx.orientation)==0) {header$orientation <- c(1,0,0,0,-1,0)
      } else  {header$orientation <- as.numeric(unlist(strsplit(image.ref[,idx.orientation],"\\\\")))}
    } else {
      header$orientation <- as.numeric(unlist(strsplit(image.ref[,idx.orientation],"\\\\")))
    }
    
    header$xyz0 <- matrix(0,ncol=3, nrow=1)
    header$xyz.from.ijk <- 0
    header$k.idx <- 0 
    header$missing.k.idx <- FALSE
    header$cube.idx  <- 0
    header$beam.source <-  matrix(c(0,0,0,1),ncol=4,byrow = TRUE)
    header$beam.orientation <- matrix(c(1,0,0,0,0,-1,0,0), ncol=2)
    header$beam.direction <-  matrix(c(0,0,-1,0),ncol=4,byrow = TRUE)
    header$beam.isocenter <- matrix(c(0,0,0,1),ncol=4,byrow = TRUE)
    if (!is.null(image.ref$upper.left.position)) 
      header$xyz0[1,1:2] <- as.numeric(unlist(strsplit(image.ref$upper.left.position,"\\\\")))
    if (!is.null(distance$source.gantry)){
      if (!is.null(distance$source.image)) header$xyz0[1,3] <- distance$source.gantry[1] - distance$source.image[1]
      header$beam.source[1,3] <- as.numeric(distance$source.gantry[1])
    }
    header$xyz.from.ijk <-.xyz.from.ijk.create (header$orientation, header$dxyz, header$xyz0[1, ])

    header$cube.idx <- matrix (c(0,0,0,1,
                                 header$n.ijk[1]-1, 0, 0, 1,
                                 header$n.ijk[1]-1, header$n.ijk[2]-1, 0, 1,
                                 0, header$n.ijk[2]-1, 0, 1,
                                 0, 0, header$k.idx[length (header$k.idx)], 1,
                                 header$n.ijk[1]-1, 0, header$k.idx[length (header$k.idx)], 1,
                                 header$n.ijk[1]-1, header$n.ijk[2]-1, header$k.idx[length (header$k.idx)], 1,
                                 0, header$n.ijk[2]-1, header$k.idx[length (header$k.idx)], 1), nrow=4, byrow= FALSE)
    
    row.names(header$cube.idx) <- c ("i","j","k","t")
    
    
    ##intercept

    if (is.null(pixel.info$intercept)) pixel.info$intercept <- 0
    ##pente

    if (is.null(pixel.info$slope)) pixel.info$slope <-1

    error.alloc <- FALSE
    if (any (is.na (match(c ("signed", "bit.stored", "bit.allocated", "MSB"), 
                          colnames(pixel.info))))) {
      header$error <- c(header$error,"voxel allocation error")
      error.alloc <- TRUE
    } else if (pixel.info$bit.stored >pixel.info$bit.allocated) {
      header$error <- c(header$error,"voxel allocation error")
      error.alloc <- TRUE
    }
    
    
    index.map <-grep("[(]7FE0,0010[)]$", name)
    if (length(index.map)!=0) {
      if (is.numeric(data[[index.map]])){
        range <-  range (data[[index.map]]*pixel.info$slope[1] +  pixel.info$intercept[1], na.rm=TRUE)
        range[range==Inf | range==-Inf] <- NA
        header$min.pixel <- range[1]
        header$max.pixel <- range[2]
      } else if (!error.alloc) {
        
        byte.nb <- length(data[[index.map]])/prod(header$n.ijk)
        if (byte.nb %in% c(1,2,4)){
          m <- data[[index.map]]
          endian <- address[index.map,3]
          data[[index.map]] <- readBin (m, what="integer", n= length(m)/byte.nb, 
                                       size = byte.nb, endian = endian)
          data[[index.map]] [is.na(data[[index.map]])] <- -2^(pixel.info$bit.stored[1])
          flag <- data[[index.map]] <0
          if (pixel.info$signed[1]==0) 
            data[[index.map]] [flag] <-  data[[index.map]] [flag] + 2^(pixel.info$bit.stored[1])
          range <-  range (data[[index.map]]*pixel.info$slope[1] +  pixel.info$intercept[1], na.rm=TRUE)
          range[range==Inf | range==-Inf] <- NA
          header$min.pixel <- range[1]
          header$max.pixel <- range[2]
        }
      }
      
      if (prod(header$n.ijk) != length(data[[index.map]])) {
        header$error <- c(header$error,"number of voxels error")
        error.alloc <- TRUE
      } 
      
    }
    
    header$pixeldecode <- list (MSB=pixel.info$MSB, bit.stored=pixel.info$bit.stored, 
                                signed=pixel.info$signed, slope=pixel.info$slope, 
                                intercept=pixel.info$intercept )
    header$machine.info <- .get.tag.info(db.name = info.name, 
                                         db.tag = info.tag, num.tag = info.num.tag, data = data, data.name = name)
    
    st.idx <- grep ("^[(]3002,0030[)] item[[:digit:]]+$", name)
    exposure.nb <- as.numeric(gsub("^[(]3002,0030[)] item","",name[st.idx]))
    if (length(exposure.nb)>0){
      header$exposure <- list()
      for(exp.nb in exposure.nb){
        header$exposure[[exp.nb]] <- list ()
        header$exposure[[exp.nb]]$info <- .get.tag.info (db.name = exposure.info.name, 
                                                         db.tag = exposure.info.ref.tag, 
                                                         num.tag = exposure.info.numtag, post.tag = name[st.idx[exp.nb]],
                                                         data = data, data.name = name)
        
        n.idx <- grep (paste0("^[(]3002,0030[)] item",exp.nb," [(]300A,00B6[)] item[[:digit:]]+$"),name)
        header$exposure[[exp.nb]]$leaf.jaw <- .get.tag.info (db.name = coll.name, 
                                                             db.tag = coll.tag, 
                                                             num.tag = coll.numtag, post.tag = name[n.idx],
                                                             data = data, data.name = name)
        
        header$exposure[[exp.nb]]$beam.ref <-  .get.tag.info (db.name = beam.ref.name, 
                                                              db.tag = beam.ref.tag, 
                                                              num.tag = beam.ref.numtag, post.tag = name[st.idx[exp.nb]],
                                                              data = data, data.name = name)
        
        
        header$exposure[[exp.nb]]$fixed.ref  <-.get.tag.info (db.name = fixed.ref.name, 
                                                              db.tag = fixed.ref.tag, 
                                                              num.tag = fixed.ref.numtag, post.tag = name[st.idx[exp.nb]],
                                                              data = data, data.name = name)
        
        header$exposure[[exp.nb]]$table.ref  <- .get.tag.info (db.name = table.ref.name, 
                                                               db.tag = table.ref.tag, 
                                                               num.tag = table.numtag, post.tag = name[st.idx[exp.nb]],
                                                               data = data, data.name = name)
        
        header$exposure[[exp.nb]]$gantry.ref <- .get.tag.info (gantry.ref.name, gantry.ref.tag, 
                                                               gantry.ref.num.tag, post.tag = name[st.idx[exp.nb]],
                                                               data = data, data.name = name)
      }
    }
    
    
    # header$distance <- distance 
    # header$gantry.ref <- .get.tag.info (gantry.ref.name, gantry.ref.tag, 
    #                              gantry.ref.num.tag, data = data, data.name = name)
    patient.ref <- .get.tag.info (db.name = patient.ref.name,
                                         db.tag = patient.ref.tag,
                                         num.tag = patient.ref.numtag, data = data,
                                         data.name = name)
    fixed.ref <- .get.tag.info (db.name = fixed.ref.name,
                                       db.tag = fixed.ref.tag,
                                        num.tag = fixed.ref.numtag, data = data, data.name = name)

    
    
    # changement de repère dans le référentiel patient
    theta <- 0
    if (!is.null(fixed.ref$gantry.angle)) theta <- fixed.ref$gantry.angle
    origin <- c(0,0,0)
    if (!is.null(patient.ref$machine.isocenter.position)) 
      origin = as.numeric(unlist(strsplit(patient.ref$machine.isocenter.position,"\\\\")))
    orientation <-switch (header$patient.position,
                          "HFS" = {c(1,0,0, 0,0,1)},
                          "HFP" = {c(-1,0,0, 0,0,1)},
                          "HFDL" = {c(0,-1,0, 0,0,1)},
                          "HFDR" = {c(0,1,0, 0,0,1)},
                          "FFS" = {c(-1,0,0, 0,0,-1)},
                          "FFP" = {c(1,0,0, 0,0,-1)},
                          "FFDL" = {c(0,1,0, 0,0,-1)},
                          "FFDR" = {c(0,-1,0, 0,0,-1)},
                          {c(1,0,0, 0,0,1)})
    

    M <- .ref.create (orientation = c(cos(theta * pi / 180),0, 
                                      -sin(theta * pi / 180), 0,1,0)) %*%
                 .ref.create (orientation = orientation, origin = origin)
    
    
    # header$beam.orientation <- -c(M[1:3,1],M[1:3,2])
    M <- solve(M)

    header$xyz.from.ijk <- M %*% header$xyz.from.ijk
    header$xyz.from.ijk[abs(header$xyz.from.ijk) < 1.0e-6] <- 0
    header$xyz0 <- matrix (matrix (c(0, 0, 0, 1),ncol=4, byrow = TRUE) %*% 
                             t(header$xyz.from.ijk)[ ,1:3],ncol=3)
    header$orientation <- c(header$xyz.from.ijk[1:3, 1]/header$dxyz[1], 
                            header$xyz.from.ijk[1:3, 2]/header$dxyz[2])
    
    header$beam.source <- matrix((header$beam.source %*% t(M))[1,1:3],ncol = 3)
    header$beam.orientation <- as.numeric((M  %*% header$beam.orientation)[1:3,])
    header$beam.direction <-  matrix((header$beam.direction %*% t(M))[1,1:3],ncol = 3)
    header$beam.isocenter <- matrix((header$beam.isocenter %*% t(M))[1,1:3],ncol = 3)
    
    d <- list()
    d[[1]]<- data
    a <- list()
    a [[1]] <- address
    
    if (only.header){ L[[header$object.alias]] <- list (header=header, from.dcm = TRUE)
    } else {L[[header$object.alias]] <- list (header=header, address=a, data=d, from.dcm = TRUE)}
    
  }
  return(L)
}
#############################################################################################
.rtstruct.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  # if (object.info$temp!= ""){
  #   tempf <- file.path (OUTDIR, paste ("temp", unlist(strsplit(object.info$temp,";")), sep = ""))
  #   data <- lapply(tempf,function(f) qread (f))
  address.all <- lapply(data, function(l) l$address)
  data.all <- lapply(data, function(l) l$data)
  filename.all <- sapply (data, function(l) l$filename)
  nb.of.struct <- length(address.all)
  L <- list ()
  for (struct.idx in 1:length(address.all)) {
    address <- address.all[[struct.idx]]
    data <- data.all[[struct.idx]]
    filename <- filename.all [struct.idx]
    
    name <- names(data)
    header <- list()
    
    header$patient <- trimws (tryCatch (data[[grep("^[(]0010,0020[)]$",name)]],error = function (e) ""))
    header$patient.name <- tryCatch (data[[grep("^[(]0010,0010[)]$",name)]],error = function (e) "")
    if (is.na(header$patient.name)) header$patient.name <- ""
    header$patient.name <- trimws(header$patient.name)                               
    header$patient.bd <- tryCatch (data[[grep("^[(]0010,0030[)]$",name)]],error = function (e) "")
    header$patient.sex <- toupper(trimws (tryCatch (data[[grep("^[(]0010,0040[)]$",name)]],error = function (e) "")))
     
    header$file.basename <- ""
    header$file.dirname <- ""
    header$object.name <- object.info$outfilename
    header$object.alias <- ""
    header$ref.object.alias <- ""
    
    header$object.info <- list ()
    header$ref.object.info <- list ()
    # header$ref.object.info$SOP.label <- sort(unique (unlist (data[grepl("^[(]3006,0010[)]",name) & grepl("[(]0008,1155[)]$",name) &
    #                                                                grepl("[(]3006,0016[)] item",name) & grepl("[(]3006,0014[)] item",name)])))
    # header$ref.object.info$SOP.ID  <- unique(as.character (sort (unlist (data[grepl("^[(]3006,0010[)]",name) & grepl("[(]0008,1150[)]$",name) &
    #                                                                             grepl("[(]3006,0016[)] item",name) & grepl("[(]3006,0014[)] item",name)]))))
    
    
    header$ref.object.info$SOP.ID <- unique(as.character (sort (unlist (data[grepl("^[(]3006,0010[)]",name) & grepl("[(]0008,1150[)]$",name) &
                                                                               grepl("[(]3006,0016[)] item",name) & grepl("[(]3006,0014[)] item",name)]))))
    if (length(header$ref.object.info$SOP.ID)>0) header$ref.object.info$SOP.ID <- header$ref.object.info$SOP.ID[header$ref.object.info$SOP.ID!=""]
    header$ref.object.info$SOP.label <- sort(unique (unlist (data[grepl("^[(]3006,0010[)]",name) & grepl("[(]0008,1155[)]$",name) &
                                                                    grepl("[(]3006,0016[)] item",name) & grepl("[(]3006,0014[)] item",name)])))
    if (length(header$ref.object.info$SOP.label)>0) header$ref.object.info$SOP.label <- header$ref.object.info$SOP.label[header$ref.object.info$SOP.label!=""]
    
    if (length(header$ref.object.info$SOP.ID)==0) {
      header$ref.object.info <- NULL
      header$ref.object.alias <- NULL
    }
    
    header$object.info$SOP.ID <- object.info$SOP.ID
    header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
    header$object.info$implementation.ID <- object.info$implementation.ID
    header$object.info$SOP.type <- object.info$SOP.type
    header$object.info$study.ID <- object.info$study.ID
    header$object.info$study.UID <- object.info$study.UID
    header$object.info$serie.UID <- object.info$serie.UID
    header$object.info$scanning.sequence <- object.info$scanning.sequence
    header$object.info$SOP.label <- data[[grep("^[(]0008,0018[)]$",name)]]
    header$object.info$encoding <- tryCatch (data[[grep("^[(]0008,0005[)]$",name)]],error = function (e) "")
    conv_idx <- integer(0)
    if(header$object.info$encoding!="") {
      conv_idx <- grep(paste0("^",tolower (gsub("[[:space:],_,-]", "", header$object.info$encoding)),"$"),
                       tolower (gsub("[[:space:],_,-]", "", iconvlist())))}
      
    if (Rdcm.mode) {
      header$object.info$dicom.file  <- sort(basename(filename))
      }
    header$object.info$nb.of.subobject <- nb.of.struct
    
    
    header$frame.of.reference <-  object.info$reference
    header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
    
    header$modality <- castlow.str (object.info$modality)
    header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")

    
    # header$acq.date <- tryCatch (data[[grep("^[(]3006,0008[)]$",name)]],error = function (e) "")
    header$acq.date <- object.info$acquisition.date
    header$study.date <- tryCatch (data[[grep("^[(]0008,0020[)]$",name)]],error = function (e) "")
    header$creation.date <- tryCatch (data[[grep("^[(]0008,0012[)]$",name)]],error = function (e) "")
    header$study.time <- tryCatch (data[[grep("^[(]0008,0030[)]$",name)]],error = function (e) "")
    header$approval.status <-  tryCatch (data[[grep("^[(]300E,0002[)]$",name)]],error = function (e) "")
    
    header$number <- as.numeric(tryCatch (data[[grep("^[(]0020,0013[)]$",name)]],error = function (e) 1))
    if (Rdcm.mode) {
      header$file.basename <- paste(object.info$outfilename, header$number, ".Rdcm", sep="")
    } else {
      header$file.basename <- sort(basename(filename))
      header$file.dirname <- unique(dirname(filename))
    }
    header$object.alias <- paste(object.info$outfilename, header$number, sep="")
    
    
    header$nb.of.roi <- sum(grepl("^[(]3006,0020[)] item",name) & grepl("[(]3006,0026[)]$",name), na.rm=TRUE)
    
    header$thickness <- 0
    header$roi.info  <- NULL
    header$ref.from.contour <- matrix (c (rep ( c (1.0, 0.0, 0.0, 0.0,0.0), 3), 1.0),nrow=4)
    
    if (header$nb.of.roi > 0) {
      
      
      
      #
      header$roi.info <- data.frame(matrix("",nrow=header$nb.of.roi,ncol=6), stringsAsFactors = FALSE)
      colnames (header$roi.info) <- c("number", "name", "description", "generation.algorithm", "color", "roi.pseudo")
      
      idx.scale <- grep ("^[(]3006,0020[)] item", name)
      L_ <- lapply (idx.scale, function (i) data[[i]])
      names(L_) <-  gsub ("^[(]3006,0020[)]", "", name[idx.scale])
      idx.roi <- as.numeric(sapply(names(L_), function(l) unlist (strsplit (l,"[ ]|item"))[3]))
      
      flag.idx <- grep("[(]3006,0022[)]$", names (L_))
      if (length (flag.idx)>0) header$roi.info$number[idx.roi[flag.idx]] <- as.character(unlist(L_[flag.idx]))
      
      flag.idx <- grep("[(]3006,0026[)]$", names (L_))
      if (length (flag.idx)) header$roi.info$name[idx.roi[flag.idx]] <- as.character(unlist(L_[flag.idx]))

      if (length(conv_idx)>0) header$roi.info$name <- iconv(header$roi.info$name,iconvlist()[conv_idx[1]])
      
      
      flag.idx <- grep("[(]3006,0028[)]$", names (L_))
      if (length (flag.idx)) header$roi.info$description[idx.roi[flag.idx]] <- as.character(unlist(L_[flag.idx]))
      
      flag.idx <- grep("[(]3006,0036[)]$", names (L_))
      header$roi.info$generation.algorithm[idx.roi[flag.idx]] <- as.character(unlist(L_[flag.idx]))
      
      
      idx.scale <- grep("^[(]3006,0039[)] item", name)
      L_ <- lapply( idx.scale, function (i) data[[i]])
      names(L_) <-  gsub("^[(]3006,0039[)]","" ,name[idx.scale])
      
      idx.roi <- as.numeric(sapply(names(L_), function(l) unlist(strsplit(l,"[ ]|item"))[3]))
      
      flag.idx <- grep("[(]3006,002A[)]$", names (L_))
      header$roi.info$color[idx.roi[flag.idx]] <- as.character(sapply (as.character(unlist(L_[flag.idx])), function(st) {
        l <-gregexpr("[+-]?([0-9]*[.])?[0-9]+", st)[[1]]
        if (l[1] == -1) return ("#000000")
        s <- paste(as.hexmode (as.numeric (sapply(1:length(l),
                                                 function(i) substring(st, l[i], l[i] + attr(l, "match.length")[i] - 1)))),
                  collapse="")
        return ( paste0("#",paste(rep("0", 6-nchar(s)), collapse =""),s))
      }))
      
      
      header$roi.info$roi.pseudo <- tolower (gsub("[[:space:],_]", "", iconv (header$roi.info$name,  to="ASCII//TRANSLIT")))
      
      
      
      db.name <- c ("(3006,0080)(3006,0082)","(3006,0080)(3006,0084)","(3006,0080)(3006,0085)",
                    "(3006,0080)(3006,0086) item1 (0008,0100)","(3006,0080)(3006,0086) item1 (0008,0102)",
                    "(3006,0080)(3006,0086) item1 (0008,0103)","(3006,0080)(3006,0086) item1 (0008,0104)",
                    "(3006,0080)(3006,00A4)","(3006,0080)(3006,00A")
      
      
      db.tag <- do.call(rbind.data.frame, lapply (db.name, function(l) unlist(strsplit(l,"[)][(]"))))
      colnames(db.tag) <- NULL   
      
      header$roi.obs <- array("", dim=c(header$nb.of.roi,length(db.name)))
      colnames(header$roi.obs) <- c("nb", "roi.nb", "label", "code.value", 
                                    "code.scheme", "code.scheme.v","code.meaning","type", "interpreter")
      rc.db <- expand.grid(1:header$nb.of.roi, 1:length(db.name))
      tag.list <- apply(rc.db,1, function(rc){paste0(db.tag[as.numeric(rc[2]),1],") item",
                                                     as.numeric(rc[1])," (",
                                                     db.tag[as.numeric(rc[2]),2])})
      ma <- match(tag.list, name)
      not.na <- which(!is.na(ma))
      header$roi.obs[as.matrix(rc.db[not.na,])] <- sapply (data[ma[not.na]],paste,collapse="\\")
      
      header$roi.obs <- data.frame(header$roi.obs,stringsAsFactors = F)
      header$roi.obs[is.na(header$roi.obs)] <- ""
      if (length(conv_idx)>0) {
        header$roi.obs$label <- iconv(header$roi.obs$label,iconvlist()[conv_idx[1]])
        header$roi.obs$code.meaning <- iconv(header$roi.obs$code.meaning,iconvlist()[conv_idx[1]])
        header$roi.obs$interpreter <- iconv(header$roi.obs$interpreter,iconvlist()[conv_idx[1]])
      }
    }
    d <- list()
    d[[1]]<- data
    a <- list()
    a [[1]] <- address
    
    if (only.header){ L[[header$object.alias]] <- list (header=header, from.dcm = TRUE)
    } else {L[[header$object.alias]] <- list (header=header, address=a, data=d, from.dcm = TRUE)}
  }
  return(L)
}
############################################################################################
.rtplan.save<- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  address.all <- lapply(data, function(l) l$address)
  data.all <- lapply(data, function(l) l$data)
  filename.all <- sapply (data, function(l) l$filename)
  nb.of.rtplan <- length(address.all)
  L <- list ()
  for (rtplan.idx in 1:length(address.all)) {
    address <- address.all[[rtplan.idx]]
    data <- data.all[[rtplan.idx]]
    filename <- filename.all [rtplan.idx]
    
    name <- names(data)
    header <- list()
    
    header$patient <- trimws (tryCatch (data[[grep("^[(]0010,0020[)]$",name)]],error = function (e) ""))
    header$patient.name <- tryCatch (data[[grep("^[(]0010,0010[)]$",name)]],error = function (e) "")
    if (is.na(header$patient.name)) header$patient.name <- ""
    header$patient.name <- trimws(header$patient.name)                               
    header$patient.bd <- tryCatch (data[[grep("^[(]0010,0030[)]$",name)]],error = function (e) "")
    header$patient.sex <- toupper(trimws (tryCatch (data[[grep("^[(]0010,0040[)]$",name)]],error = function (e) "")))
    
    header$file.basename <- ""
    header$file.dirname <- ""
    header$object.name <- object.info$outfilename
    header$object.alias <- ""
    header$ref.object.alias <- ""
    
    header$object.info <- list ()
    header$ref.object.info <- list ()
    
    header$object.info$SOP.ID <- object.info$SOP.ID
    header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
    header$object.info$implementation.ID <- object.info$implementation.ID
    header$object.info$SOP.type <- object.info$SOP.type
    header$object.info$study.ID <- object.info$study.ID
    header$object.info$study.UID <- object.info$study.UID
    header$object.info$serie.UID <- object.info$serie.UID
    header$object.info$scanning.sequence <- object.info$scanning.sequence
    header$object.info$SOP.label <- data[[grep("^[(]0008,0018[)]$",name)]]
    header$object.info$encoding <- tryCatch (data[[grep("^[(]0008,0005[)]$",name)]],error = function (e) "")
    if (Rdcm.mode) header$object.info$dicom.file  <- sort(basename(filename))
    header$object.info$nb.of.subobject <- nb.of.rtplan
    
    
    header$ref.object.info <- list ()
    post.tag <- c (name[grep ("^[(]300C,0002[)] item[[:digit:]]+$", name)], 
                   name[grep ("^[(]300C,0060[)] item[[:digit:]]+$", name)],
                   name[grep ("^[(]300C,0080[)] item[[:digit:]]+$", name)])
    db.name <- c( "SOP.ID", "SOP.label", "rtplan.relation.ship")
    dg.tag <- c("(0008,1150)", "(0008,1155)", "(300A,0055)")
    num.tag <- numeric(0)
    reference.info <- .get.tag.info (db.name = db.name, 
                                     db.tag = dg.tag, 
                                     num.tag = num.tag, post.tag = post.tag,
                                     data = data, data.name = name)
    if (is.null(reference.info)){
      header$ref.object.info <- NULL
      header$ref.object.alias <- NULL
    } else {
      header$ref.object.info$SOP.ID <- as.character(unique (reference.info$SOP.ID))
      header$ref.object.info$SOP.label <- as.character(unique (reference.info$SOP.label))
    }
    
    header$frame.of.reference <-  object.info$reference
    header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
    
    header$modality <- castlow.str (object.info$modality)
    header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
    
    
    # header$acq.date <- tryCatch (data[[grep("^[(]300A,0006[)]$",name)]],error = function (e) "")
    header$acq.date <- object.info$acquisition.date
    header$study.date <- tryCatch (data[[grep("^[(]0008,0020[)]$",name)]],error = function (e) "")
    header$creation.date <- tryCatch (data[[grep("^[(]0008,0012[)]$",name)]],error = function (e) "")
    header$study.time <- tryCatch (data[[grep("^[(]0008,0030[)]$",name)]],error = function (e) "")
    
    header$approval.status <-  tryCatch (data[[grep("^[(]300E,0002[)]$",name)]],error = function (e) "")
    header$number <- rtplan.idx
    if (Rdcm.mode) {
      header$file.basename <- paste(object.info$outfilename, header$number, ".Rdcm", sep="")
    } else {
      header$file.basename <- sort(basename(filename))
      header$file.dirname <- unique(dirname(filename))
    }
    header$object.alias <- paste(object.info$outfilename, header$number, sep="")

    header <- c(header, .rtplan.beam.field(data))
   
    d <- list()
    d[[1]]<- data
    a <- list()
    a [[1]] <- address
    
    if (only.header){ L[[header$object.alias]] <- list (header=header, from.dcm = TRUE)
    } else {L[[header$object.alias]] <- list (header=header, address=a, data=d, from.dcm = TRUE)}
  }
  return(L)  
}

############################################################################################
.rtplan.beam.field <- function(data) {
  name <- names(data)
  header <- list()
  header$plan.info <-  .reduc.tab (data.frame (
    label = tryCatch (data[[grep ("^[(]300A,0002[)]$", name)]],error = function (e) ""),
    plan.name = tryCatch (data[[grep ("^[(]300A,0003[)]$", name)]],error = function (e) ""),
    plan.description = tryCatch (data[[grep ("^[(]300A,0004[)]$", name)]],error = function (e) ""),
    tt.protocol = tryCatch (data[[grep ("^[(]300A,0004[)]$", name)]],error = function (e) ""),
    plan.intent = tryCatch (data[[grep ("^[(]300A,000A[)]$", name)]],error = function (e) ""),
    tt.site = tryCatch (data[[grep ("^[(]300A,000B[)]$", name)]],error = function (e) ""),
    geometry = tryCatch (data[[grep ("^[(]300A,000C[)]$", name)]],error = function (e) "")
  ))
  
  st.idx <- grep ("^[(]300A,0010[)] item[[:digit:]]+$", name)
  st.nb <- as.numeric(gsub("^[(]300A,0010[)] item","",name[st.idx]))
  header$presc.dose <- .reduc.tab (do.call(rbind , lapply (st.nb, function (st.idx){
    str <- paste0("^[(]300A,0010[)] item",st.idx)
    return(data.frame(
      ref.roi.nb = tryCatch (data[[grep (paste(str,"[(]3006,0084[)]$" ),name)]],error = function (e) ""),
      dose.ref.nb =tryCatch (data[[grep (paste(str,"[(]300A,0012[)]$" ),name)]], error = function (e) ""),
      dose.ref.id =tryCatch (data[[grep (paste(str,"[(]300A,0013[)]$" ),name)]], error = function (e) ""),
      struct.type =tryCatch (data[[grep (paste(str,"[(]300A,0014[)]$" ),name)]], error = function (e) ""),
      description = tryCatch (data[[grep (paste(str,"[(]300A,0016[)]$" ),name)]], error = function (e) ""),
      pt.coord = tryCatch (data[[grep (paste(str,"[(]300A,0018[)]$" ),name)]], error = function (e) ""),
      nominal.prior.dose= as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,001A[)]$" ),name)]], error = function (e) "")),
      dose.type = tryCatch (data[[grep (paste(str,"[(]300A,0020[)]$" ),name)]],  error = function (e) ""),
      constraint.weight = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0021[)]$" ),name)]], error = function (e) "")),
      deliv.warn.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0022[)]$" ),name)]], error = function (e) "")),
      deliv.max.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0023[)]$" ),name)]], error = function (e) "")),
      targ.min.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0025[)]$" ),name)]], error = function (e) "")),
      targ.presc.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0026[)]$" ),name)]], error = function (e) "")),
      targ.max.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0027[)]$" ),name)]], error = function (e) "")),
      targ.underdose.vol.frac = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0028[)]$" ),name)]], error = function (e) "")),
      org.risk.full.vol.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,002A[)]$" ),name)]], error = function (e) "")),
      org.risk.lim.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,002B[)]$" ),name)]], error = function (e) "")),
      org.risk.max.dose = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,002C[)]$" ),name)]], error = function (e) "")),
      org.risk.overdose.vol.frac = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,002D[)]$" ),name)]], error = function (e) ""))
    ))
  })))
  
  
  
  st.idx <- grep ("^[(]300A,0070[)] item[[:digit:]]+$", name)
  st.nb <- as.numeric(gsub("^[(]300A,0070[)] item","",name[st.idx]))
  fraction.info <- .reduc.tab(do.call (rbind ,lapply(st.nb, function (st.idx){
    str <- paste0("^[(]300A,0070[)] item",st.idx)
    data.frame (
      fraction.id = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0071[)]$"),name)]], error = function (e) "")),
      description = tryCatch (data[[grep (paste(str,"[(]300A,0072[)]$"),name)]], error = function (e) ""),
      nb.of.frac.planned = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0078[)]$"),name)]], error = function (e) "")),
      nb.of.frac.pattern.digit.per.day = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0079[)]$"),name)]], error = function (e) "")),
      repeat.frac.cycle.le = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,007A[)]$"),name)]], error = function (e) "")),
      frac.pattern = tryCatch (data[[grep (paste(str,"[(]300A,007B[)]$"),name)]], error = function (e) ""),
      nb.of.beam = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,0080[)]$"),name)]], error = function (e) "")),
      beam.dose.meaning = tryCatch (data[[grep (paste(str,"[(]300A,008B[)]$"),name)]], error = function (e) ""),
      nb.of.brachy.app = as.numeric(tryCatch (data[[grep (paste(str,"[(]300A,00A0[)]$"),name)]], error = function (e) "")))})))
  
  
  
  
  header$fraction.info<- fraction.info
  
  
  patient.setup  <- .get.tag.info(db.name = c("beam.nb","patient.position", "position.label"), 
                                  db.tag = c("(300A,0182)","(0018,5100)","(300A,0183)"), 
                                  num.tag = 1, 
                                  post.tag =  name[grep("^[(]300A,0180[)] item[[:digit:]]+$",name)],
                                  data = data, data.name = name)
  
  if (sum(header$fraction.info$nb.of.beam)>0){
    
    fraction.beam <- .reduc.tab(do.call(rbind,lapply(st.nb, function (st.idx){
      str <- paste0("^[(]300A,0070[)] item",st.idx)
      beam.idx <- grep (paste(str,"[(]300C,0004[)] item[[:digit:]]+$" ),name)
      beam.nb <- as.numeric(gsub(paste(str,"[(]300C,0004[)] item"),"",name[beam.idx]))
      beam.info <- lapply (beam.nb, function(beam.idx){
        str_ <- paste0(str," [(]300C,0004[)] item",beam.idx)
        data.frame (
          fraction.id = as.numeric(fraction.info$fraction.id[st.idx]),
          nb.of.frac.planned = as.numeric(fraction.info$nb.of.frac.planned[st.idx]),
          # beam.uid = tryCatch (data[[grep (paste(str_,"[(]300A,0083[)]$"),name)]], error = function (e) ""),
          beam.dose = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300A,0084[)]$"),name)]], error = function (e) "")),
          beam.specif.pt  = tryCatch (data[[grep (paste(str_,"[(]300A,0082[)]$"),name)]], error = function (e) ""),
          beam.meterset = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300A,0086[)]$"),name)]], error = function (e) "")),
          beam.type = tryCatch (data[[grep (paste(str_,"[(]300A,0090[)]$"),name)]], error = function (e) ""),
          alt.dose = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300A,0091[)]$"),name)]], error = function (e) "")),
          alt.type = tryCatch (data[[grep (paste(str_,"[(]300A,0092[)]$"),name)]], error = function (e) ""),
          duration.lim = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300A,00C5[)]$"),name)]], error = function (e) "")),
          beam.nb = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300C,0006[)]$"),name)]], error = function (e) ""))
        )})
      
      do.call(rbind ,beam.info)
    })))
    header$fraction.beam<- fraction.beam
    
    post.tag <- paste0("^[(]",unique(substring(name[grep("[(]300A,00C0[)]$", name)],2,10)),"[)]")
    post.tag <- name[grep (paste (post.tag,"item[[:digit:]]+$"), name)]
    
    
    db.name <- c ("beam.nb","beam.name","beam.description","beam.type","radiation.type",
                  "high.dose.technique.type","treatment.machine.name","device.serial.nb",                                 
                  "primary.dosimeter.unit", "referenced.tolerance.table.nb",                   
                  "src.axis.dist","referenced.patient.setup.nb","treatment.delivery.type","wedges.nb",                                      
                  "compensators.nb","total.compensator.tray.factor",                    
                  "boli.nb","blocks.nb", "total.block.tray.factor","final.cumul.meterset.weight",                      
                  "ctl.pts.nb","radiation.mass.nb", "radiation.atomic.nb","radiation.charge.state",                           
                  "scan.mode","modulated.scan.mode.type", "virtual.src.axis.dist","total.wedge.tray.water.equ.thickness",      
                  "total.compensator.tray.water.equ.thickness","total.block.tray.water.equ.thickness",      
                  "range.shifters.nb","lateral.spreading.devices.nb", "range.modulators.nb","fixation.light.azimuthal.angle",                   
                  "fixation.light.polar.angle")  
    
    db.tag <- c("(300A,00C0)", "(300A,00C2)", "(300A,00C3)", "(300A,00C4)", 
                "(300A,00C6)", "(300A,00C7)", "(300A,00B2)", "(0018,1000)",
                "(300A,00B3)", "(300C,00A0)", "(300A,00B4)", "(300C,006A)", 
                "(300A,00CE)", "(300A,00D0)", "(300A,00E0)", "(300A,00E2)",
                "(300A,00ED)", "(300A,00F0)", "(300A,00F2)", "(300A,010E)", 
                "(300A,0110)", "(300A,0302)", "(300A,0304)", "(300A,0306)",
                "(300A,0308)", "(300A,0309)", "(300A,030A)", "(300A,00D7)", 
                "(300A,02E3)", "(300A,00F3)", "(300A,0312)", "(300A,0330)",
                "(300A,0340)", "(300A,0356)", "(300A,0358)")
    
    num.tag <- c(1,10:12,14:24,27:35)
    rtbeam.info <- .get.tag.info(db.name = db.name, 
                                 db.tag = db.tag, num.tag = num.tag, post.tag =  post.tag,
                                 data = data, data.name = name)
    m <- unique(c(1,match (patient.setup$beam.nb,rtbeam.info$beam.nb)))
    m <- m[!is.na(m)]
    rtbeam.info$patient.position <- "" 
    rtbeam.info$patient.position[m] <- patient.setup$patient.position
    rtbeam.info[,c("beam.nb","patient.position")] <- .fill.db (rtbeam.info[,c("beam.nb","patient.position")])
    header$beam.info  <- rtbeam.info
    
    # m <- match (fraction.beam$beam.nb,rtbeam.info$beam.nb)
    # header$beam.info  <- cbind(fraction.beam,rtbeam.info[m,-1])
    
    
    
    ctl.pt.tag <- c("(0018,0060)",
                    "(300A,0112)","(300A,0114)","(300A,0115)","(300A,011E)",
                    "(300A,011F)","(300A,0120)","(300A,0121)","(300A,0122)",
                    "(300A,0123)","(300A,0124)","(300A,0125)","(300A,0126)",
                    "(300A,0128)","(300A,0129)","(300A,012A)","(300A,012C)",
                    "(300A,012E)","(300A,0130)","(300A,0132)","(300A,0133)",
                    "(300A,0134)","(300A,0140)","(300A,0142)","(300A,0144)",
                    "(300A,0146)","(300A,014A)","(300A,014C)",
                    "(300A,0151)","(300A,030D)","(300A,035A)",
                    "(300A,0390)","(300A,0392)","(300A,0394)","(300A,0395)",
                    "(300A,0396)","(300A,0398)","(300A,039A)")
    ctl.pt.name <- c("KVP",
                     "ctl.pt.idx","energy","dose.rate","gantry.angle","gantry.rot","beam.lim.angle","beam.lim.rot",
                     "pat.support.angle","pat.support.rot","table.top.ecc.dist","table.top.ecc.angle",
                     "table.top.ecc.rot","table.top.vert.pos","table.top.long.pos","table.top.lat.pos",
                     "isocenter","surf.xyz0","src.surf.dist","src.ext.ctr.dist","ext.ctr.entry.pt",
                     "cum.meterset.weight","table.top.pitch.angle","table.top.pitch.rot",
                     "table.top.roll.angle","table.top.roll.rot","gantry.pitch.angle","gantry.pitch.rot",
                     "chair.head.frame.position", "snout.position", "meterset.rate",
                     "scan.spot.id","scan.spot.pos.nb","scan.spot.pos.map","scan.spot.reordering.allowed",
                     "scan.spot.meterset.weight","scanning.spot.size","painting.nb")
    ctl.pt.num.tag <- c(1:5,7,9,11:12,14:16,19:23,25,27,30:31,33,36:38)
    
    coll.name <- c("type", "source.distance", "nb","position.boundaries","position")
    coll.tag <- c("(300A,00B8)", "(300A,00BA)", "(300A,00BC)", "(300A,00BE)", "(300A,011C)")
    coll.numtag <- c(2:4)
    
    header$beam.ctl.pt <-  lapply( 1:nrow(header$beam.info)   , function(i){
      p.t <- gsub(")","[)]",gsub("(","[(]",post.tag[i], fixed=TRUE), fixed=TRUE)
      p.t <- name[grepl (paste0 ("^",p.t," [(]300A,0111[)] item[[:digit:]]+$"), name)|  
                    grepl (paste0 ("^",p.t," [(]300A,03A8[)] item[[:digit:]]+$"), name)]
      if (length(p.t)==0) return(NULL)
      info <- .get.tag.info(db.name = ctl.pt.name, 
                            db.tag = ctl.pt.tag, num.tag = ctl.pt.num.tag, post.tag =  p.t,
                            data = data, data.name = name)
      info <- .fill.db(info)
      
      orientation <-switch (header$beam.info$patient.position[i],
                            "HFS" = {c(1,0,0, 0,0,1)},
                            "HFP" = {c(-1,0,0, 0,0,1)},
                            "HFDL" = {c(0,-1,0, 0,0,1)},
                            "HFDR" = {c(0,1,0, 0,0,1)},
                            "FFS" = {c(-1,0,0, 0,0,-1)},
                            "FFP" = {c(1,0,0, 0,0,-1)},
                            "FFDL" = {c(0,1,0, 0,0,-1)},
                            "FFDR" = {c(0,-1,0, 0,0,-1)},
                            {c(1,0,0, 0,0,1)})
      
      beam.source0 <-c(0,0,0,1)  
      if (!is.null(header$beam.info$src.axis.dist)) beam.source0[3] <- header$beam.info$src.axis.dist[i]
      
      L <- lapply(1:nrow(info),function(j){
        beam.source <- beam.source0
        beam.orientation <- matrix(c(1,0,0,0, 0,-1,0,0, 0,0,-1,0), ncol=3)
        beam.direction <-c(0,0,-1,0) 
        beam.isocenter <- c(0,0,0,1)  
        if (!is.null(info$snout.position)) if(info$snout.position[j]!="NA") 
          beam.source[3] <-info$snout.position[j]
        
        M <- .ref.create (orientation = c(cos(info$gantry.angle[j] * pi / 180),
                                          0, 
                                          -sin(info$gantry.angle[j] * pi / 180), 0,1,0))
        
        M <- M %*% .ref.create (orientation = c(cos(info$pat.support.angle[j] * pi / 180),
                                                -sin(info$pat.support.angle[j] * pi / 180),0,
                                                sin(info$pat.support.angle[j] * pi / 180),
                                                cos(info$pat.support.angle[j] * pi / 180),0))
        
        origin <- c(0,0,0)
        if (!is.null(info$isocenter)) if(info$isocenter[j]!="NA") 
          origin <-  as.numeric(unlist(strsplit(info$isocenter[j],"\\\\")))
        
        M <- M %*% .ref.create (orientation = orientation, origin = origin)
        M <- solve(M)
        beam.source <-(beam.source %*%  t(M)) [,1:3]
        beam.source[abs(beam.source) < 1e-10] <- 0
        beam.orientation <- as.numeric((M  %*% beam.orientation)[1:3,1:2])
        beam.orientation[abs(beam.orientation) < 1e-10] <- 0
        beam.direction <-(beam.direction %*%  t(M)) [,1:3]
        beam.direction[abs(beam.direction) < 1e-10] <- 0
        beam.isocenter <-(beam.isocenter %*%  t(M)) [,1:3]
        beam.isocenter[abs(beam.isocenter) < 1e-10] <- 0
        
        spot.map <- NULL
        if (!is.null(info$scan.spot.pos.map)) {
          
          spot.map<- cbind(matrix(as.numeric(unlist(strsplit(info$scan.spot.pos.map[j],",|c[(]|[)]"))[-1]), 
                                  ncol=2, byrow=TRUE),0,1)
          spot.map <-  (spot.map %*% t(M))
          spot.map <- cbind(rep(info$ctl.pt.idx[j],nrow(spot.map)),spot.map)
          colnames(spot.map) <- c("ctl.pt.idx","x","y","z","t")
        } 
        
        list(beam.orientation = beam.orientation, beam.source= beam.source, 
             beam.direction=beam.direction, beam.isocenter=beam.isocenter, 
             spot.map = spot.map)
      })
      
      beam.orientation <- as.matrix(do.call(rbind,lapply(L,function(l) l$beam.orientation)))
      beam.source <- as.matrix(do.call(rbind,lapply(L,function(l) l$beam.source)))
      colnames(beam.source) <- c("x","y","z")
      beam.direction <- as.matrix(do.call(rbind,lapply(L,function(l) l$beam.direction)))
      colnames(beam.direction) <- c("x","y","z")
      beam.isocenter <- as.matrix(do.call(rbind,lapply(L,function(l) l$beam.isocenter)))
      colnames(beam.isocenter) <- c("x","y","z")
      spot.map <- do.call(rbind,lapply(L,function(l) l$spot.map))
      if (!is.null(spot.map)) spot.map <- spot.map[,1:4]
      ########################
      p.t_ <- gsub(")","[)]",gsub("(","[(]",p.t, fixed=TRUE), fixed=TRUE)
      p.t2 <- paste0 ("^",p.t_," [(]300A,011A[)] item[[:digit:]]+$")
      leaf.jaw <- do.call(rbind,lapply(1:length(p.t2), function(i){
        db <- .get.tag.info(db.name = coll.name, 
                            db.tag = coll.tag, num.tag = coll.numtag, post.tag =  name[grepl (p.t2[i],name)],
                            data = data, data.name = name, adjust=FALSE)
        
        if (!is.null(.reduc.tab(db))) db$ctl.pt.idx <- info$ctl.pt.idx[i]
        db
      }))
      leaf.jaw <- .reduc.tab(leaf.jaw)
      if(!is.null(leaf.jaw)) {
        leaf.jaw <- leaf.jaw[ , c(ncol(leaf.jaw),1:(ncol(leaf.jaw)-1))]
        n <- unique(leaf.jaw$type)
        leaf.jaw <- lapply(n, function(type) leaf.jaw[leaf.jaw$type == type, c(1,3)])
        names(leaf.jaw) <- n
        return(list(info=info, leaf.jaw = leaf.jaw, beam.source = beam.source, 
                    beam.direction = beam.direction, 
                    beam.orientation = beam.orientation,
                    beam.isocenter = beam.isocenter))
      }
      
      
      ###########################
      if (!is.null(spot.map)) 
        return (list(info=info, beam.source = beam.source, beam.direction = beam.direction,
                     beam.orientation = beam.orientation, 
                     beam.isocenter = beam.isocenter,
                     spot.map = spot.map))
      return (list(info=info, beam.source = beam.source, beam.direction = beam.direction,
                   beam.orientation = beam.orientation, beam.isocenter = beam.isocenter))
    })
    if (!is.null(rtbeam.info$beam.name)) names(header$beam.ctl.pt) <- rtbeam.info$beam.name
    
  }
  
  
  if (sum(header$fraction.info$nb.of.brachy.app)>0){
    
    st.idx <- grep ("^[(]300A,0070[)] item[[:digit:]]+$", name)
    st.nb <- as.numeric(gsub("^[(]300A,0070[)] item","",name[st.idx]))
    header$fraction.brachy <- .reduc.tab(do.call (rbind, lapply(st.nb, function (st.idx){
      str <- paste0("^[(]300A,0070[)] item",st.idx)
      brach.idx <- grep (paste(str,"[(]300C,000A[)] item[[:digit:]]+$" ),name)
      brach.nb <- as.numeric(gsub(paste(str,"[(]300C,000A[)] item"),"",name[brach.idx]))
      brachy.app <- lapply (brach.nb, function(brach.idx){
        str_ <- paste0(str," [(]300C,000A[)] item",brach.idx)
        data.frame (
          fraction.id = as.numeric(fraction.info$fraction.id[st.idx]),
          nb.of.frac.planned = as.numeric(fraction.info$nb.of.frac.planned[st.idx]),
          brachy.dose = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300A,00A4[)]$"),name)]], error = function (e) "")),
          brachy.nb = as.numeric(tryCatch (data[[grep (paste(str_,"[(]300C,000C[)]$"),name)]], error = function (e) "")),
          brachy.specif.pt = tryCatch (data[[grep (paste(str_,"[(]300A,00A2[)]$"),name)]], error = function (e) "")
          
        )})
      do.call(rbind ,brachy.app)
    })))
  }
  return(header)
}


#############################################################################################

.other.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  # if (object.info$temp!= ""){
  #   tempf <- file.path (OUTDIR, paste ("temp", unlist(strsplit(object.info$temp,";")), sep = ""))
  #   data <- lapply(tempf,function(f) qread (f))
  address <- lapply(data, function(l) l$address)
  filename <- sapply (data, function(l) l$filename)
  data <- lapply(data, function(l) l$data)

  name <- names(data[[1]])
  
  header <- list()
  header$patient <- trimws (tryCatch (data[[1]][[grep("^[(]0010,0020[)]$",name)]],error = function (e) ""))
  header$patient.name <- tryCatch (data[[1]][[grep("^[(]0010,0010[)]$",name)]],error = function (e) "")
  if (is.na(header$patient.name)) header$patient.name <- ""
  header$patient.name <- trimws(header$patient.name)                               
  header$patient.bd <- tryCatch (data[[1]][[grep("^[(]0010,0030[)]$",name)]],error = function (e) "")
  header$patient.sex <- toupper(trimws (tryCatch (data[[1]][[grep("^[(]0010,0040[)]$",name)]],error = function (e) "")))
  
  header$file.basename <- ""
  if (Rdcm.mode) {
    header$file.basename <- paste(object.info$outfilename, ".Rdcm", sep="")
  } else {
    header$file.basename <- sort(basename(filename))
    header$file.dirname <- unique(dirname(filename))
  }
  header$file.dirname <- ""
  header$object.name <- object.info$outfilename
  header$object.alias <- object.info$outfilename
  header$ref.object.alias <- ""

  
  
  header$object.info <- list ()
  header$object.info$SOP.ID <- object.info$SOP.ID
  header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
  header$object.info$implementation.ID <- object.info$implementation.ID
  header$object.info$SOP.type <- object.info$SOP.type
  header$object.info$study.ID <- object.info$study.ID
  header$object.info$study.UID <- object.info$study.UID
  header$object.info$serie.UID <- object.info$serie.UID
  header$object.info$scanning.sequence <- object.info$scanning.sequence
  header$object.info$SOP.label <- unique(tryCatch (sapply(data, function (d) d[[grep("^[(]0008,0018[)]$",names(d))]]),error = function (e)  ""))
  header$object.info$encoding <- tryCatch (data[[1]][[grep("^[(]0008,0005[)]$",name)]],error = function (e) "")
  
  header$object.info$nb.of.subobject <- NA
  
  header$ref.object.info <- list ()
  header$ref.object.info$SOP.ID <- unique(unlist(lapply(data, function (d) tryCatch (d[grep("[(]0008,1150[)]$",names(d))],error = function (e)  ""))))
  if (length(header$ref.object.info$SOP.ID)>0) header$ref.object.info$SOP.ID <- header$ref.object.info$SOP.ID[header$ref.object.info$SOP.ID!=""]
  header$ref.object.info$SOP.label <- sort (unique(unlist(lapply(data, function (d) tryCatch (d[grep("[(]0008,1155[)]$",names(d))],error = function (e)  "")))))
  if (length(header$ref.object.info$SOP.label)>0) header$ref.object.info$SOP.label <- header$ref.object.info$SOP.label[header$ref.object.info$SOP.label!=""]
  
  if (length(header$ref.object.info$SOP.ID)==0) {
    header$ref.object.info <- NULL
    header$ref.object.alias <- NULL
  }
  
  header$frame.of.reference <-  object.info$reference
  header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
  
  header$modality <- castlow.str (object.info$modality)
  header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
  
  
  # header$acq.date <- unique(sapply(data, function (d) tryCatch (d[[grep("^[(]0008,0023[)]$",names(d))]],error = function (e)  "")))
  header$acq.date <- object.info$acquisition.date
  header$study.date <- unique(sapply(data, function (d) tryCatch (d[[grep("^[(]0008,0020[)]$",names(d))]],error = function (e)  "")))
  header$creation.date <- unique(sapply(data, function (d) tryCatch (d[[grep("^[(]0008,0012[)]$",names(d))]],error = function (e)  "")))
  header$study.time <- unique(sapply(data, function (d) tryCatch (d[[grep("^[(]0008,0030[)]$",names(d))]],error = function (e)  "")))
  header$approval.status <- unique(sapply(data, function (d) tryCatch (d[[grep("^[(]300E,0002[)]$",names(d))]],error = function (e)  "")))

  
  
  L <- list ()
  if (only.header){ L[[header$object.alias]] <- list (header=header, from.dcm = TRUE)
  } else {L[[header$object.alias]] <- list (header=header, address=address, data=data, from.dcm = TRUE)}
  return(L)
  
}

#############################################################################################
.reg.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE){
  # if (object.info$temp!= ""){
  #   tempf <- file.path (OUTDIR, paste ("temp", unlist(strsplit(object.info$temp,";")), sep = ""))
  #   data <- lapply(tempf,function(f) qread (f))
  address <- lapply(data, function(l) l$address)
  filename <- sapply (data, function(l) l$filename)
  data <- lapply(data, function(l) l$data)
  
  name <- names(data[[1]])
  
  header <- list()
  header$patient <- trimws (tryCatch (data[[1]][[grep("^[(]0010,0020[)]$",name)]],error = function (e) ""))
  header$patient.name <- tryCatch (data[[1]][[grep("^[(]0010,0010[)]$",name)]],error = function (e) "")
  if (is.na(header$patient.name)) header$patient.name <- ""
  header$patient.name <- trimws(header$patient.name)                               
  header$patient.bd <- tryCatch (data[[1]][[grep("^[(]0010,0030[)]$",name)]],error = function (e) "")
  header$patient.sex <- toupper(trimws (tryCatch (data[[1]][[grep("^[(]0010,0040[)]$",name)]],error = function (e) "")))
  
  header$file.basename <- ""
  if (Rdcm.mode) {
    header$file.basename <- paste(object.info$outfilename, ".Rdcm", sep="")
  } else {
    header$file.basename <- sort(basename(filename))
    header$file.dirname <- unique(dirname(filename))
  }
  header$file.dirname <- ""
  header$object.name <- object.info$outfilename
  header$object.alias <- object.info$outfilename
  
  header$object.info <- list ()
  header$object.info$SOP.ID <- object.info$SOP.ID
  header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
  header$object.info$implementation.ID <- object.info$implementation.ID
  header$object.info$SOP.type <- object.info$SOP.type
  header$object.info$study.ID <- object.info$study.ID
  header$object.info$study.UID <- object.info$study.UID
  header$object.info$serie.UID <- object.info$serie.UID
  header$object.info$scanning.sequence <- object.info$scanning.sequence
  header$object.info$SOP.label <- data[[1]][[grep("^[(]0008,0018[)]$",name)]]
  header$object.info$encoding <- tryCatch (data[[1]][[grep("^[(]0008,0005[)]$",name)]],error = function (e) "")
  if (Rdcm.mode) header$object.info$dicom.file <-sort(basename(filename))
  
  header$object.info$nb.of.subobject <- length(data)
  
  header$frame.of.reference <-  object.info$reference
  header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
  
  header$modality <- castlow.str (object.info$modality)
  header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
  
  
  # header$acq.date <- tryCatch (data[[1]][[grep("^[(]0008,0023[)]$",name)]],error = function (e) "")
  header$acq.date <- object.info$acquisition.date
  header$study.date <- tryCatch (data[[1]][[grep("^[(]0008,0020[)]$",name)]],error = function (e) "")
  header$creation.date <- tryCatch (data[[1]][[grep("^[(]0008,0012[)]$",name)]],error = function (e) "")
  header$study.time <- tryCatch (data[[1]][[grep("^[(]0008,0030[)]$",name)]],error = function (e) "")
  
  header$nb.of.ref <- sum(grepl("^[(]0070,0308[)] item",name) & grepl("[(]0020,0052[)]$",name), na.rm=TRUE)
  
  if (header$nb.of.ref>0){
    
    idx.scale <- grep("^[(]0070,0308[)] item", name)
    L_ <- lapply( idx.scale, function (i) data[[1]][[i]])
    names(L_) <-  gsub("^[(]0070,0308[)]","" ,name[idx.scale])
    
    idx.ref <- as.numeric(sapply(names(L_), function(l) unlist(strsplit(l,"[ ]|item"))[3]))
    
    src.idx <- grep("[(]0020,0052[)]$", names (L_))
    src <- as.character(unlist(L_[src.idx]))
    
    type.idx <- grep("[(]0070,030C[)]$", names (L_))
    type <- as.character(unlist(L_[type.idx]))
    
    mat.idx <-  grep("[(]3006,00C6[)]$", names (L_))
    mat<- as.character(unlist(L_[mat.idx]))
    sep.mat <-unique(unlist(strsplit(tolower(mat[1]),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep.mat <- sep.mat[sep.mat!=""]
    mat <- lapply(mat, function(st) {
      matrix (as.numeric(unlist(strsplit(st,sep.mat, fixed = TRUE))), ncol=4, byrow=TRUE)
    })
    
    header$ref.data <- list ()
    header$ref.data<- lapply (unique(idx.ref), function (idx) {
      idx.ref_ <- which(idx.ref[mat.idx] == idx)
      list( src = src[idx.ref_], type =  type[idx.ref_], matrix = mat[[idx.ref_]])
    })
    names (header$ref.data) <- src
    
    
  }
  L <- list ()
  if (only.header){ L[[header$object.alias]] <- list (header=header, from.dcm = TRUE)
  } else {  L[[header$object.alias]] <- list (header=header, address=address, data=data, from.dcm = TRUE)}
  return(L)
  
}

#############################################################################################
.rtdose.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE) {
  # if (object.info$temp!= ""){
  #   tempf <- file.path (OUTDIR, paste ("temp", unlist(strsplit(object.info$temp,";")), sep = ""))
  #   data <- lapply(tempf,function(f) qread (f))
  address <- lapply(data, function(l) l$address)
  filename <- sapply (data, function(l) l$filename)
  data <- lapply(data, function(l) l$data)
 
  
  header <- list()
  header$patient <- trimws (tryCatch (data[[1]][[grep("^[(]0010,0020[)]$",names (data[[1]]))]],error = function (e) ""))
  header$patient.name <- tryCatch (data[[1]][[grep("^[(]0010,0010[)]$",names (data[[1]]))]],error = function (e) "")
  if (is.na(header$patient.name)) header$patient.name <- ""
  header$patient.name <- trimws(header$patient.name)                               
  header$patient.bd <- tryCatch (data[[1]][[grep("^[(]0010,0030[)]$",names(data[[1]]))]],error = function (e) "")
  header$patient.sex <- toupper(trimws (tryCatch (data[[1]][[grep("^[(]0010,0040[)]$",names(data[[1]]))]],error = function (e) "")))
  
  header$file.basename <- ""
  header$file.dirname <- ""
  header$object.name <- object.info$outfilename
  header$object.alias <- ""
  header$ref.object.alias <- ""
  
  
  header$object.info <- list ()
  header$object.info$SOP.ID <- object.info$SOP.ID
  header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
  header$object.info$implementation.ID <- object.info$implementation.ID
  header$object.info$SOP.type <- object.info$SOP.type
  header$object.info$study.ID <- object.info$study.ID
  header$object.info$study.UID <- object.info$study.UID
  header$object.info$serie.UID <- object.info$serie.UID
  header$object.info$scanning.sequence <- object.info$scanning.sequence
  header$object.info$SOP.label <- ""
  header$object.info$encoding <- tryCatch (data[[1]][[grep("^[(]0008,0005[)]$",names(data[[1]]))]],error = function (e) "")
  if (Rdcm.mode) header$object.info$dicom.file <- ""
  
  image.nb <- as.numeric( sapply(data, function (d) tryCatch (d[[grep("[(]300C,0006[)]$",names(d))]],error = function (e) 1)))
  header$object.info$nb.of.subobject <- length(unique(image.nb))
  
  header$ref.object.info <- list ()
  
  
  header$frame.of.reference <-  object.info$reference
  header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
  
  header$modality <- castlow.str (object.info$modality)
  header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
  
  
  header$acq.date <- ""
  header$study.date <- ""
  header$creation.date <-  ""
  header$study.time <- ""
  
  
  L <- list ()
  for (img.idx in  unique(image.nb)) {
    header_ <- list ()
    selection <- image.nb==img.idx
    dat <- data [selection][[1]]
    add <- address [selection][[1]]
    fn <- filename [selection]
    
    name <- names (dat)
    
    header_$number <- img.idx
    if (Rdcm.mode) {
      header$file.basename <- paste(object.info$outfilename,img.idx, ".Rdcm", sep="")
      header$object.info$dicom.file  <- sort(basename(fn))
    } else {
      header$file.basename <- sort(basename(fn))
      header$file.dirname <- unique(dirname(fn))
    }
    header$object.alias <- paste(object.info$outfilename, img.idx, sep="")
    
    
    # header$acq.date <- tryCatch (dat[[grep("^[(]0008,0023[)]$",name)]],error = function (e) "")
    header$acq.date <- object.info$acquisition.date
    header$study.date <- tryCatch (dat[[grep("^[(]0008,0020[)]$",name)]],error = function (e) "")
    header$creation.date <- tryCatch (dat[[grep("^[(]0008,0012[)]$",name)]],error = function (e) "")
    header$study.time <- tryCatch (dat[[grep("^[(]0008,0030[)]$",name)]],error = function (e) "")
    
    header$object.info$SOP.label <- sort(unique(unlist(dat[grepl("^[(]0008,0018[)]$",name)])))
    
    
    # header$ref.object.info$SOP.ID <- unique (unlist(dat[grep("[(]0008,1150[)]$", name)]))
    # header$ref.object.info$SOP.label <- sort (unique (unlist(dat[grep("[(]0008,1155[)]$", name)])))
    
    header$ref.object.info$SOP.ID <- unique (unlist(dat[grep("[(]0008,1150[)]$", name)]))
    if (length(header$ref.object.info$SOP.ID)>0) header$ref.object.info$SOP.ID <- header$ref.object.info$SOP.ID[header$ref.object.info$SOP.ID!=""]
    header$ref.object.info$SOP.label <-  sort (unique (unlist(dat[grep("[(]0008,1155[)]$", name)])))
    if (length(header$ref.object.info$SOP.label)>0) header$ref.object.info$SOP.label <- header$ref.object.info$SOP.label[header$ref.object.info$SOP.label!=""]
    
    if (length(header$ref.object.info$SOP.ID)==0) {
      header$ref.object.info <- NULL
      header$ref.object.alias <- NULL
    }
    #
    beam.idx <- grep("[(]300C,0006[)]$", name)
    if (length(beam.idx)>0) {header_$beam.number <- as.numeric(unlist(dat[beam.idx]))
    } else {header_$beam.nb <- NA}
    
    # Dose_Summation_Type
    Dose.Summation.Type <-  tryCatch (dat[[grep("[(]3004,000A[)]$", name)]],error = function (e) "")
    if (Dose.Summation.Type=="") header$warning <- c(header$warning,"no (3004,000A) for Dose Summation Type")
    
    #unit
    header_$unit <- tryCatch (dat[[grep("^[(]3004,0002[)]$", name)]],error = function (e) "")
    if (header_$unit=="") {header$warning <- c(header$warning,"no specified unit")
    } else {
      if (tolower (trimws (header_$unit))!="gy") {
        header$warning <- c(header$warning,"(3004,0002) should be Gy for map")
      } else { 
        header_$unit <- "Gy"
      }
    }
    physical <- tryCatch (dat[[grep("^[(]3004,0004[)]$", name)]],error = function (e) "")
    if (physical=="") {header$warning <- c(header$warning,"not Physical dose")
    }else {if (tolower (trimws (physical))!="physical") header$warning <- c(header$warning,"(3004,0004) should be Physical dose for map")
    }
    
    vol.exist <-    tryCatch (any(grepl("^[(]7FE0,0010[)]$", name)),error = function (e) FALSE)
    
    nz <- as.numeric(tryCatch (dat[[grep("^[(]0028,0008[)]$", name)]],error = function (e) 0))
    if (nz==0) header$error <- c(header$error,"(0028,0008) number of slice error")
    ## nombre de lignes et colonnes
    ny <-  as.numeric(tryCatch (dat[[grep("^[(]0028,0010[)]$", name)]],error = function (e) 0))
    if (ny==0) header$error <- c(header$error,"(0028,0010) number of line error")
    nx <- as.numeric(tryCatch (dat[[grep("^[(]0028,0011[)]$", name)]],error = function (e) 0))
    if (nx==0) header$error <- c(header$error,"(0028,0011) number of column error")
    header_$n.ijk <- c (nx, ny, nz)
    
    ## slice.thickness 
    header_$slice.thickness <-  0
    
    header_$min.pixel <- NA
    header_$max.pixel <- NA
    
    ##espacement des pixels en x, y
    header_$dxyz <- tryCatch (dat[[grep("[(]0028,0030[)]$", name)]],error = function (e) "")
    if (header_$dxyz =="")  header$error <- c (header$error,"(0028,0030) xy-spacing error")
    sep <-unique(unlist(strsplit(tolower(header_$dxyz),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    if (length(sep)==1)  {header_$dxyz <- c(as.numeric (unlist(strsplit(header_$dxyz,sep, fixed=TRUE))),header_$slice.thickness)
    } else {header_$dxyz <- c (0, 0, header_$slice.thickness)}
    ## referientiel patient
    header_$orientation <- tryCatch (dat[[grep("[(]0020,0037[)]$", name)]],error = function (e) "")
    if (header_$orientation =="")  header$error <- c(header$error,"(0020,0037) patient orientation error")
    sep <-unique(unlist(strsplit(tolower(header_$orientation),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    ref.machine <- NULL
    if (length(sep)==1){
      header_$orientation <- as.numeric (unlist(strsplit(header_$orientation,sep, fixed=TRUE)))
      ref.machine <-  (.ref.create (header_$orientation))
    } else {header_$orientation <- c(1,0,0,0,1,0)}
    
    ## z coordinates
    z <- tryCatch (dat[[grep("[(]3004,000C[)]$", name)]],error = function (e) "0")
    if (z == "0" & header_$n.ijk[3]>1)  {
      header$error <- c(header$error,"no (3004,000C) definition")
    }
    sep <-unique(unlist(strsplit(tolower(z),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    if (length(sep)==1)  {z <- as.numeric (unlist(strsplit(z,sep, fixed=TRUE)))
    } else {z=0}
    
    if (length(z)>1) {
      dz <- unique(round(diff(z),3))
      dz <- dz[dz!=0]
      if (length (dz)>0) {
        header_$dxyz[3] <- dz [ which(abs(dz)==min(abs(dz)))]
        header_$slice.thickness <- abs(header_$dxyz[3])
      }
    }
    
    ## position de l'image dans le référentiel patient
    header_$xyz0 <- tryCatch (dat[[grep("[(]0020,0032[)]$", name)]],error = function (e) "")
    if (any(header_$xyz0 == "")) header$error <- c(header$error,"(0020,0032) pos. in referentiel error")
    sep <- unique(unlist(strsplit(tolower(header_$xyz0[1]),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    if (length(sep)==1) {
      header_$xyz0 <- as.numeric(unlist(sapply (header_$xyz0, function(st) unlist(strsplit(st,sep,fixed=TRUE)))))
    } else {
      header_$xyz0 <- c(0,0,0)
    }                                                     
    
    header_$xyz0  <-  matrix (rep (header_$xyz0, length(z)), 
                                      ncol= 3, 
                                      byrow = TRUE,
                                      dimnames = list(NULL,c("x0","y0","z0")))
    # if (z[1]==0) {
    header_$xyz0[, 3] <- header_$xyz0[, 3] + z
    # } else {
    #   header_$xyz0[, 3] <-  z
    # }  
    
    
    header_$xyz.from.ijk <-.xyz.from.ijk.create (header_$orientation, header_$dxyz, header_$xyz0[1, ])
    if (header_$dxyz[3]!=0) {header_$k.idx <- round ((z-z[1])/header_$dxyz[3], 0 )
    } else {header_$k.idx <- 0 }
    header_$missing.k.idx <- FALSE
    if (header_$n.ijk[3]>0) header_$missing.k.idx <- any ((1:header_$n.ijk[3]) != (header_$k.idx+1))
    header_$cube.idx <- matrix ( c(0,0,0,1,
                                   header_$n.ijk[1]-1, 0, 0, 1,
                                   header_$n.ijk[1]-1, header_$n.ijk[2]-1, 0, 1,
                                   0, header_$n.ijk[2]-1, 0, 1,
                                   0, 0, header_$k.idx[length (header_$k.idx)], 1,
                                   header_$n.ijk[1]-1, 0, header_$k.idx[length (header_$k.idx)], 1,
                                   header_$n.ijk[1]-1, header_$n.ijk[2]-1, header_$k.idx[length (header_$k.idx)], 1,
                                   0, header_$n.ijk[2]-1, header_$k.idx[length (header_$k.idx)], 1), nrow=4, byrow= FALSE)
    
    row.names(header_$cube.idx) <- c ("i","j","k","t")
    
    
    ##intercept
    intercept <- 0
    
    ##pente
    slope <-  as.numeric (tryCatch (dat[[grep("[(]3004,000E[)]$", name)]], error = function (e) NA))
    if (is.na(slope)) {
      header$warning <- c(header$warning,"no slope specification")
      slope <- 1
    }
    
    error.alloc <- FALSE
    bit.allocated <- unique(as.numeric (tryCatch (dat[[grep("[(]0028,0100[)]$", name)]],error = function (e) 0)))
    bit.stored <- unique(as.numeric (tryCatch (dat[[grep("[(]0028,0101[)]$", name)]],error = function (e) 0)))
    signed <- unique(as.numeric (tryCatch (dat[[grep("[(]0028,0103[)]$", name)]],error = function (e) 0)))
    MSB <-  unique(as.numeric (tryCatch (dat[[grep("[(]0028,0102[)]$", name)]],error = function (e) 0)))
    if (length(signed)!=1){
      header$error <- c(header$error,"(0028,0103) voxel sign error")
      error.alloc <- TRUE
    }
    if (bit.stored >bit.allocated | length(bit.allocated)!=1 | length(bit.stored)!=1) {
      header$error <- c(header$error,"(0028,0100) or (0028,0101) voxel allocation error")
      error.alloc <- TRUE
    }
    if (bit.stored-1 != MSB | length(bit.stored)!=1 | length(MSB)!=1) {
      header$error <- c(header$error,"(0028,0101) or (0028,0102) voxel storage error")
      error.alloc <- TRUE
    }
    
    index.map <-grep("^[(]7FE0,0010[)]$", name)
    if (length(index.map)==0) index.map <- NA
    if (length(dat[[index.map]])==0) index.map <- NA
    if (is.na(index.map)) {
      header$error <- c(header$error,"pixel TAG (7FE0,0010) error")
      error.alloc <- TRUE
    } else {
      
      if (is.numeric(dat[[index.map]])){
        range <-  range (dat[[index.map]]*slope +  intercept, na.rm=TRUE)
        range[range==Inf | range==-Inf] <- NA
        header_$min.pixel <- range[1]
        header_$max.pixel <- range[2]
      } else if (!error.alloc) {
        
        byte.nb <- length(dat[[index.map]])/prod(header_$n.ijk)
        if (byte.nb %in% c(1,2,4)){
          m <- dat[[index.map]]
          endian <- add[index.map,3]
          dat[[index.map]] <- readBin (m, what="integer", n= length(m)/byte.nb, 
                                       size = byte.nb, endian = endian)
          dat[[index.map]] [is.na(dat[[index.map]])] <- -2^bit.stored
          flag <- dat[[index.map]] <0
          if (signed==0) dat[[index.map]] [flag] <-  dat[[index.map]] [flag] + 2^bit.stored 
          range <-  range (dat[[index.map]]*slope +  intercept, na.rm=TRUE)
          range[range==Inf | range==-Inf] <- NA
          header_$min.pixel <- range[1]
          header_$max.pixel <- range[2]
        }
      }
        
      # if (is.raw(dat[[index.map]])){
      #   byte.nb <- length(dat[[index.map]])/prod(header_$n.ijk)
      #   if (byte.nb %in% c(1,2,4)){
      #     m <- dat[[index.map]]
      #     VR <- c ("OB", "OW", "SL") [which(byte.nb==c(1,2,4))]
      #     endian <- add[index.map,3]
      #     dat[[index.map]] <- dicom.tag.parser (1, length(m), VR, endian, m)
      #   }
      # }
      
      
      # if (prod(header_$n.ijk) != length(dat[[index.map]])) {
      #   byte.nb <- length(dat[[index.map]])/prod(header_$n.ijk)
      #   if (byte.nb %in% c(1,2,4)){
      #     m <- as.raw (dat[[index.map]])
      #     VR <- c ("OB", "OW", "SL") [which(byte.nb==c(1,2,4))]
      #     endian <- add[index.map,3]
      #     dat[[index.map]] <- dicom.tag.parser (1, length(m), VR, endian, m)
      #   }
      # }
      
      if (prod(header_$n.ijk) != length(dat[[index.map]])) {
        header$error <- c(header$error,"number of voxels error")
        error.alloc <- TRUE
      } 
      
      # #min max
      # if (!error.alloc) {
      #   
      #   pixel <- tryCatch (as.numeric (dat[[index.map]]), error = function (e) NA)
      #   flag <- pixel > 2^MSB
      #   flag[is.na(flag)] <- FALSE
      #   if (signed==1) pixel [flag] <-  pixel [flag] - 2^bit.stored 
      #   dat[[index.map]] <- pixel
      #   
      #   
      #   range <-  range (dat[[index.map]]*slope +  intercept, na.rm=TRUE)
      #   range[range==Inf | range==-Inf] <- NA
      #   header_$min.pixel <- range[1]
      #   header_$max.pixel <- range[2]
      # }
    }
    
    header_$pixeldecode <- list (MSB=MSB, bit.stored=bit.stored, signed=signed, slope=slope, intercept=intercept )
    
    
    
    
    
    #DVH
    
    
    # header_$roi.info <- data.frame(matrix("",nrow=header$nb.of.roi,ncol=6), stringsAsFactors = FALSE)
    # colnames (header$roi.info) <- c("number", "name", "description", "generation.algorithm", "color", "roi.pseudo")
    
    idx.dvh.part <- grep ("^[(]3004,0050[)] item",  name)
    L_ <- lapply (idx.dvh.part, function (i) dat[[i]])
    names(L_) <-  gsub ("^[(]3004,0050[)]", "", name[idx.dvh.part])
    idx.dvh <- as.numeric(sapply(names(L_), function(l) unlist (strsplit (l,"[ ]|item"))[3]))
    
    if (length(idx.dvh)>0) {
      
      header_$dvh <- list(nb.of.dvh = length (unique (idx.dvh)))
      header_$dvh$ref.object.alias  <- ""
      header_$dvh$ref.object.info <- list (SOP.ID =tryCatch ( dat[[which(grepl("^[(]300C,0060[)]", name) & grepl("[(]0008,1150[)]$", name))]],error = function (e) ""),
                                           SOP.label = tryCatch ( dat[[which(grepl("^[(]300C,0060[)]", name) & grepl("[(]0008,1155[)]$", name))]],error = function (e) ""))
      
      
      header_$dvh$normalization.point <-  tryCatch ( dat[[grep("^[(]3004,0040[)]$", name) ]],error = function (e) "")
      l <-gregexpr("[+-]?([0-9]*[.])?[0-9]+", header_$dvh$normalization.point)[[1]]
      header_$dvh$normalization.point <- as.numeric (sapply(1:length(l),function(i) substring(header_$dvh$normalization.point, l[i], l[i] + attr(l, "match.length")[i] - 1)))
      
      header_$dvh$normalization.dose.value<-  tryCatch ( as.numeric(dat[[grep("^[(]3004,0042[)]$", name) ]]),error = function (e) 1)
      
      header_$dvh$info <- data.frame(matrix("",nrow = header_$dvh$nb.of.dvh, ncol = 12), stringsAsFactors = FALSE)
      colnames (header_$dvh$info ) <- c("number","roi.name", "contribution", "dvh.type", "dose.unit", "dose.type", "dose.scaling",
                                        "volume.unit","min.dose","max.dose","mean.dose","pixel.nb")
      
      flag.idx <- grep("[(]3006,0084[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$number[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0062[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$contribution[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0001[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$dvh.type[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0002[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$dose.unit[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0004[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$dose.type[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0052[)]$", names (L_))
      header_$dvh$info$dose.scaling <- 1
      if (length (flag.idx)>0) header_$dvh$info$dose.scaling[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0054[)]$", names (L_))
      if (length (flag.idx)>0) header_$dvh$info$volume.unit[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0070[)]$", names (L_))
      header_$dvh$info$min.dose <- 0
      if (length (flag.idx)>0) header_$dvh$info$min.dose[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0072[)]$", names (L_))
      header_$dvh$info$max.dose <- 0
      if (length (flag.idx)>0) header_$dvh$info$max.dose[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0074[)]$", names (L_))
      header_$dvh$info$mean.dose <- 0
      if (length (flag.idx)>0) header_$dvh$info$mean.dose[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      flag.idx <- grep("[(]3004,0056[)]$", names (L_))
      header_$dvh$info$pixel.nb <- 0
      if (length (flag.idx)>0) header_$dvh$info$pixel.nb[idx.dvh[flag.idx]] <- unlist(L_[flag.idx])
      
      header_$dvh$info[ ,c(7,9:12)]<- data.frame(lapply( header_$dvh$info[ ,c(7,9:12)], as.numeric), stringsAsFactors=FALSE)
      
      
      
      pt.idx <-  grep("[(]3004,0058[)]$", names (L_))
      pt<- as.character(unlist(L_[pt.idx]))
      
      idx.max <- which.max(nchar(pt))
      sep.pt <- unique(unlist(strsplit(tolower(pt[idx.max]),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
      sep.pt <- sep.pt[sep.pt!=""]
      if (length(pt)>0) {
        header_$dvh$data <- lapply(header_$dvh$info$number, function(i) NULL)
        for (idx in 1:length(pt.idx)){ 
          p <- as.numeric(unlist(strsplit(pt[idx],sep.pt, fixed = TRUE)))
          if (header_$dvh$info$pixel.nb[idx] != length(p)/2) {
            warning ("dvh loading error : number of points is not a multiple of 2." , call.=FALSE)
            
          } else {
            header_$dvh$data [[idx.dvh[pt.idx[idx]]]] =  data.frame (matrix (p, byrow=TRUE, ncol=2, dimnames = list(NULL,c ("x","y"))),
                                                                     stringsAsFactors = FALSE)
          }
        }
      }
      
      # 
      
    }
    
    ########### Save ########
    d <- list()
    d[[1]]<- dat
    a <- list()
    a [[1]] <- add
    
    if (only.header){ L[[header$object.alias]] <- list (header= do.call (c, list(header,header_)), from.dcm = TRUE)
    } else {  L[[header$object.alias]] <- list (header= do.call (c, list(header,header_)), address=a, data=d, from.dcm = TRUE)}
    
  }
  return(L)
  # }
}

#############################################################################################
.img3Dplan.save <- function (object.info, data, only.header=FALSE, Rdcm.mode=FALSE) {
  
  # if (object.info$temp!= ""){
  #   tempf <- file.path (OUTDIR, paste ("temp", unlist(strsplit(object.info$temp,";")), sep = ""))
  #   data <- lapply(tempf,function(f) qread (f))
  address <- lapply(data, function(l) l$address)
  filename <- sapply (data, function(l) l$filename)
  data <- lapply(data, function(l) l$data)
  # name <- names(data[[1]])
  
  header <- list()
  header$patient <- trimws (tryCatch (data[[1]][[grep("^[(]0010,0020[)]$",names (data[[1]]))]],error = function (e) ""))
  header$patient.name <-  tryCatch (data[[1]][[grep("^[(]0010,0010[)]$",names (data[[1]]))]],error = function (e) "")
  if (is.na(header$patient.name)) header$patient.name <- ""
  header$patient.name <- trimws(header$patient.name)
  
  header$patient.bd <- tryCatch (data[[1]][[grep("^[(]0010,0030[)]$",names(data[[1]]))]],error = function (e) "")
  header$patient.sex <- toupper(trimws (tryCatch (data[[1]][[grep("^[(]0010,0040[)]$",names(data[[1]]))]],error = function (e) "")))
  
  header$file.basename <- ""
  header$file.dirname <- ""
  header$object.name<- object.info$outfilename
  header$object.alias <- ""
  # header$ref.object.alias <- ""
  
  header$object.info <- list ()
  header$object.info$SOP.ID <- object.info$SOP.ID
  header$object.info$transfer.syntax.UID <- object.info$transfer.syntax.UID
  header$object.info$implementation.ID <- object.info$implementation.ID
  header$object.info$SOP.type <- object.info$SOP.type
  header$object.info$study.ID <- object.info$study.ID
  header$object.info$study.UID <- object.info$study.UID
  header$object.info$serie.UID <- object.info$serie.UID
  header$object.info$scanning.sequence <- object.info$scanning.sequence
  header$object.info$SOP.label <- ""
  header$object.info$encoding <- tryCatch (data[[1]][[grep("^[(]0008,0005[)]$",names(data[[1]]))]],error = function (e) "")
  if (Rdcm.mode) header$object.info$dicom.file  <- ""
  image.nb <- sapply(data, function (d) tryCatch (d[[grep("^[(]0020,0100[)]$",names(d))]],error = function (e) tryCatch (d[[grep("[(]0020,0037[)]$",names(d))]],error = function (e) 1)))
  image.nb <- match(image.nb, unique(image.nb))
  
  header$object.info$nb.of.subobject <- length(unique(image.nb))
  # header$ref.object.info <- list ()
  
  header$frame.of.reference <-  object.info$reference
  header$ref.pseudo <- paste("ref", object.info$ref.label, sep="")
  
  header$modality <- castlow.str (object.info$modality)
  header$description <- paste(object.info$study.description, object.info$serie.description, sep="|")
  
  
  header$acq.date <- ""
  header$study.date <- ""
  header$creation.date <- ""
  header$study.time <- ""
  
  L <- list ()
  for (img.idx in unique(image.nb)) {
    header_ <- list ()
    selection <- image.nb==img.idx
    dat <- data[selection]
    add <- address[selection]
    fn <- filename[selection]
    
    header_$number <- img.idx
    if (Rdcm.mode) {
      header$file.basename <- paste(object.info$outfilename,img.idx, ".Rdcm", sep="")
      header$object.info$dicom.file  <- basename(fn)
    } else {
      header$file.basename <- basename(fn)
      header$file.dirname <- unique(dirname(fn))
    }
    
    header$object.alias <- paste(object.info$outfilename, img.idx, sep="")
    
    header$object.info$SOP.label <- tryCatch (sapply(dat, function (d) d[[grep("^[(]0008,0018[)]$",names(d))]]),error = function (e)  "")
    
    # 
    # header$ref.object.info$SOP.ID <- unique(unlist(lapply(dat, function (d) 
    #   tryCatch (d[[which(grepl("^[(]0008,1120[)]", names(d)) & grepl("[(]0008,1150[)]$", names(d)))]],error = function (e)  ""))))
    # if (length(header$ref.object.info$SOP.ID)>1) header$ref.object.info$SOP.ID <- header$ref.object.info$SOP.ID[header$ref.object.info$SOP.ID!=""]
    # header$ref.object.info$SOP.label <- sort (unique(unlist(lapply(dat, function (d) 
    #   tryCatch (d[[which(grepl("^[(]0008,1120[)]", names(d)) & grepl("[(]0008,1155[)]$", names(d)))]],error = function (e)  "")))))
    # if (length(header$ref.object.info$SOP.label)>1) header$ref.object.info$SOP.label <- header$ref.object.info$SOP.label[header$ref.object.info$SOP.label!=""]

    # header$acq.date <- unique(sapply(dat, function (d) tryCatch (d[[grep("^[(]0008,0023[)]$",names(d))]],error = function (e)  "")))
    header$acq.date <- object.info$acquisition.date
    header$study.date <- unique(sapply(dat, function (d) tryCatch (d[[grep("^[(]0008,0020[)]$",names(d))]],error = function (e)  "")))
    header$creation.date <- unique(sapply(dat, function (d) tryCatch (d[[grep("^[(]0008,0012[)]$",names(d))]],error = function (e)  "")))
    header$study.time <- unique(sapply(dat, function (d) tryCatch (d[[grep("^[(]0008,0030[)]$",names(d))]],error = function (e)  "")))
 
    header_$unit <- ""
    
    nz <-  sum (sapply(dat, function (d) any(grepl("^[(]7FE0,0010[)]$",names(d)))), na.rm=TRUE)
    ## nombre de lignes et colonnes
    ny <-  as.numeric (tryCatch (dat[[1]][[grep("^[(]0028,0010[)]$",names(dat[[1]]))]],error = function (e) 0))
    if (ny==0) header$error <- c(header$error,"(0028,0010) number of line error")
    nx <-  as.numeric (tryCatch (dat[[1]][[grep("^[(]0028,0011[)]$",names(dat[[1]]))]],error = function (e) 0))
    if (nx==0) header$error <- c(header$error,"(0028,0011) number of column error")
    header_$n.ijk <- c (nx, ny, nz)
    
    ## slice.thickness 
    header_$slice.thickness <-  as.numeric(tryCatch (dat[[1]][[grep("[(]0018,0050[)]$",names(dat[[1]]))]],error = function (e) 0))
    
    header_$min.pixel <- NA
    header_$max.pixel <- NA
    
    ##espacement des pixels en x, y
    header_$dxyz <- tryCatch (dat[[1]][[grep("[(]0028,0030[)]$",names(dat[[1]]))]],error = function (e) "")
    if (header_$dxyz =="")  header$error <- c (header$error,"(0028,0030) xy-spacing error")
    sep <-unique(unlist(strsplit(header_$dxyz,"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    if (length(sep)==1)  {header_$dxyz <- c(as.numeric (unlist(strsplit(header_$dxyz,sep, fixed=TRUE))),header_$slice.thickness)
    } else { header_$dxyz <- c (0, 0, header_$slice.thickness)}
    
    ## referientiel patient
    header_$orientation <- tryCatch (dat[[1]][[grep("[(]0020,0037[)]$",names(dat[[1]]))]],error = function (e) "")
    if (header_$orientation[1] =="")  header$error <- c (header$error,"(0020,0037) patient orientation error")
    sep <-unique(unlist(strsplit(tolower(header_$orientation),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    ref.machine <- NULL
    if (length(sep)==1){
      header_$orientation <- as.numeric (unlist(strsplit(header_$orientation,sep, fixed=TRUE)))
      ref.machine <-  (.ref.create (header_$orientation))
    } else {header_$orientation <- c(1,0,0,0,1,0,0)}
    
    ## position de l'image dans le référentiel patient
    header_$xyz0 <- sapply(dat, function (d) tryCatch (d[[grep("[(]0020,0032[)]$",names(d))]],error = function (e) ""))
    if (any(header_$xyz0 == "")) header$error <- c(header$error,"(0020,0032) pos. in referentiel error")
    sep <- unique(unlist(strsplit(tolower(header_$xyz0[1]),"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep <- sep[sep!=""]
    if (length(sep)==1) {
      header_$xyz0 <- matrix (as.numeric(unlist(sapply (header_$xyz0, 
                                                                function(st) unlist(strsplit(st,sep,fixed=TRUE))))), 
                                      ncol= 3, 
                                      byrow = TRUE,
                                      dimnames = list(NULL,c("x0","y0","z0")))
      if (!is.null(ref.machine)) {
        machine.xyz0 <- header_$xyz0 %*% t (ref.machine[1:3,1:3])
        if (nrow(machine.xyz0)>1) {
          Z.order <- order(machine.xyz0[,3])
          machine.xyz0 <- machine.xyz0[Z.order,]
          
          header_$xyz0 <- header_$xyz0[Z.order, ]
          header$object.info$SOP.label  <- header$object.info$SOP.label[Z.order] 
          if (Rdcm.mode) {
            header$object.info$dicom.file  <- header$object.info$dicom.file[Z.order] 
          } else {
            header$file.basename <-  header$file.basename[Z.order] 
          }
          
          #on corrige l'ordre des DICOM
          dat <- dat[Z.order]
          add <- add[Z.order]
          dz <- machine.xyz0[2:nrow(machine.xyz0),3]- machine.xyz0[1:(nrow(machine.xyz0)-1),3]
          dz0 <- min(abs(dz))
          header_$dxyz[3] <- mean(dz[round(abs(dz)/dz0,0)==1])
          
        }
        header_$xyz.from.ijk <-.xyz.from.ijk.create (header_$orientation, header_$dxyz, header_$xyz0[1, ])
        if (nrow(machine.xyz0)>1){ header_$k.idx <- round((machine.xyz0[ ,3] - machine.xyz0[1,3])/header_$dxyz[3],0)
        } else {header_$k.idx <- 0}
        header_$missing.k.idx <- FALSE
        if (header_$n.ijk[3]>0) header_$missing.k.idx <- any ((1:header_$n.ijk[3]) != (header_$k.idx+1))
        
        header_$cube.idx <- matrix ( c(0,0,0,1,
                                       header_$n.ijk[1]-1, 0, 0, 1,
                                       header_$n.ijk[1]-1, header_$n.ijk[2]-1, 0, 1,
                                       0, header_$n.ijk[2]-1, 0, 1,
                                       0, 0, header_$k.idx[length (header_$k.idx)], 1,
                                       header_$n.ijk[1]-1, 0, header_$k.idx[length (header_$k.idx)], 1,
                                       header_$n.ijk[1]-1, header_$n.ijk[2]-1, header_$k.idx[length (header_$k.idx)], 1,
                                       0, header_$n.ijk[2]-1, header_$k.idx[length (header_$k.idx)], 1), nrow=4, byrow= FALSE)
        
        row.names(header_$cube.idx) <- c ("i","j","k","t")
      }
    }
    
    ##intercept
    intercept <-  as.numeric(sapply(dat, function (d) tryCatch (d[[grep("[(]0028,1052[)]$",names(d))]],error = function (e) 0)))
    
    ##pente
    slope <-  as.numeric(sapply(dat, function (d) tryCatch (d[[grep("[(]0028,1053[)]$",names(d))]],error = function (e) 1)))
    
    error.alloc <- FALSE
    bit.allocated <- unique(as.numeric(sapply(dat, function (d) tryCatch (d[[grep("^[(]0028,0100[)]$",names(d))]],error = function (e) 0))))
    bit.stored <- unique(as.numeric(sapply(dat, function (d) tryCatch (d[[grep("^[(]0028,0101[)]$",names(d))]],error = function (e) 0))))
    signed <- unique(as.numeric(sapply(dat, function (d) tryCatch (d[[grep("^[(]0028,0103[)]$",names(d))]],error = function (e) 0))))
    MSB <-  unique(as.numeric(sapply(dat, function (d) tryCatch (d[[grep("^[(]0028,0102[)]$",names(d))]],error = function (e) 0))))
    if (length(signed)!=1){
      header$error <- c(header$error,"(0028,0103) voxel sign error")
      error.alloc <- TRUE
    }
    if (bit.stored >bit.allocated | length(bit.allocated)!=1 | length(bit.stored)!=1) {
      header$error <- c(header$error,"(0028,0100) or (0028,0101) voxel allocation error")
      error.alloc <- TRUE
    }
    if (bit.stored-1 != MSB | length(bit.stored)!=1 | length(MSB)!=1) {
      header$error <- c(header$error,"(0028,0101) or (0028,0102) voxel storage error")
      error.alloc <- TRUE
    }
    
    index.map <- sapply(1:length(dat), function (i){
      id <- grep("^[(]7FE0,0010[)]$",names(dat[[i]]))
      if (length(id)==0) return(NA)
      if (length(dat[[i]][[id]])==0) return(NA)
      return(id)
    })
    
    if (any (is.na(index.map))) {
      header$error <- c(header$error,"pixel TAG (7FE0,0010) error")
      error.alloc <- TRUE
    } else {
      
      if (is.numeric(dat[[1]][[index.map[1]]])){
        range <- range(unlist(lapply(1:length(dat), function (i){
          r <- range (dat[[i]][[index.map[i]]], na.rm=TRUE)
          r[r==Inf | r==-Inf] <- NA
          r* slope[i] +   intercept[i]
        })), na.rm = TRUE)
        range[range==Inf | range==-Inf] <- NA
        header_$min.pixel <- range[1]
        header_$max.pixel <- range[2]
      } else if (!error.alloc) {
        byte.nb <- sapply(1:length(dat), function (i) {length(tryCatch (dat[[i]][[index.map[i]]], error = function (e) NA))})/prod(header_$n.ijk[1:2])
        
        for(map.idx in 1:length(dat)) {
          if (byte.nb[map.idx] %in% c(1,2,4)){
            m <- dat[[map.idx]][[index.map[map.idx]]]
            endian <- add[[map.idx]][index.map[map.idx],3]
            dat[[map.idx]][[index.map[map.idx]]] <- readBin (m, what="integer", 
                                                             n= length(m)/byte.nb[map.idx], 
                                                             size = byte.nb, endian = endian)
            dat[[map.idx]][[index.map[map.idx]]] [is.na(dat[[map.idx]][[index.map[map.idx]]] )] <- -2^bit.stored
          flag <- dat[[map.idx]][[index.map[map.idx]]]  <0
          if (signed==0) dat[[map.idx]][[index.map[map.idx]]][flag] <-  dat[[map.idx]][[index.map[map.idx]]][flag] + 2^bit.stored 
          
          }
        }
        
        if (prod(header_$n.ijk) != 
            sum(sapply(1:length(dat), function (i) {length(tryCatch (as.numeric (dat[[i]][[index.map[i]]]), error = function (e) NA))}))) {    
          header$error <- c(header$error,"number of voxels error")
          error.alloc <- TRUE
        } 
      
        range <- range(unlist(lapply(1:length(dat), function (i){
          r <- range (dat[[i]][[index.map[i]]], na.rm=TRUE)
          r[r==Inf | r==-Inf] <- NA
          r* slope[i] +   intercept[i]
        })), na.rm = TRUE)
        range[range==Inf | range==-Inf] <- NA
        header_$min.pixel <- range[1]
        header_$max.pixel <- range[2]
      }

    }
    header_$pixeldecode <- list (MSB=MSB, bit.stored=bit.stored, signed=signed, slope=slope, intercept=intercept )
    
    if (only.header){ L[[header$object.alias]] <- list (header= do.call (c, list(header,header_)), from.dcm = TRUE)
    } else {  L[[header$object.alias]] <- list (header= do.call (c, list(header,header_)), address=add, data=dat, from.dcm = TRUE)}
  }
  return (L)
  # }
}
#############################################################################################

.ref.create <- function(orientation, origin= c (0, 0, 0)){
  TM <- NULL
  
  if (length(orientation)==6){
    TM <- as.matrix (cbind (c (orientation[1:3],0),c (orientation[4:6], 0),
                  c (vector.product (orientation[1:3], orientation[4:6]),0),c(0,0,0,1)))
  }
  if (length(orientation)==9){
    TM <- as.matrix (cbind (c (orientation[1:3],0), c (orientation[4:6], 0),
                              c (orientation[7:9],0), c (0,0,0,1)))
  }
  if (is.null(TM)) return (NULL)
  TM <- solve(TM)
  TM[1:3,4] <- -(TM  %*% c(origin,1))[1:3,1]
  return(TM)
}
#############################################################################################
.xyz.from.ijk.create <- function(orientation_vect, dxyz, corner.pt){
  orientation <- matrix(orientation_vect,ncol=2, byrow=FALSE)
  
  as.matrix(cbind(c(orientation[,1],0) * dxyz[1],
                  c(orientation[,2],0) * dxyz[2],
                  c(vector.product(orientation[,1],orientation[,2]),0) * dxyz[3],
                  c(corner.pt,1)),dimnames = list(NULL,NULL))
  
}


#############################################################################################
# T.MAT fonctions
#############################################################################################

.load.T.MAT.by.reglist <- function (reg.list){
  
  .transfert.ref.tab.create  <- function (reg.list) {
    db <-as.data.frame (matrix(nrow = 0, ncol=6,
                               dimnames=list (NULL, c("file", "dest.ref", "src.ref", "type", "matrix.index", "creation.date"))) ,stringsAsFactors=FALSE)
    if (is.null(reg.list)) return (db)
    if (length(reg.list)==0) return (db)
    reg.ref<-  lapply (1:length(reg.list), function(idx) {
      n.target <- reg.list[[idx]]$nb.of.ref
      if (length(n.target)==0) return (NULL)
      if (n.target==0) return (NULL)
      db <- do.call(rbind.data.frame, lapply (reg.list[[idx]]$ref.data, function (lr) c (lr$src, trimws (lr$type))))
      colnames(db) <- c("src.ref",  "type")
      db <- cbind( data.frame (file=rep (names(reg.list)[idx],nrow(db)), dest.ref = rep(reg.list[[idx]]$frame.of.reference,nrow(db))),
                   db, data.frame (matrix.index= 1:nrow(db)), creation.date = rep (reg.list[[idx]]$creation.date,nrow(db)))

      db$path <- sub("^[/]","",file.path(reg.list[[idx]]$file.dirname,reg.list[[idx]]$file.basename))
      db[] <- data.frame (lapply(db [], as.character), stringsAsFactors=FALSE)
      return (db)
    })
    reg.ref <- do.call (rbind, reg.ref)
    # reg.ref <- reg.ref[reg.ref$dest.ref!=reg.ref$src.ref,]
    rownames(reg.ref) <- NULL
    
    return(reg.ref)
  }
  ref.info <- unique (do.call (rbind.data.frame, 
                               lapply(reg.list, 
                                      function(l) list(ref.pseudo=l$ref.pseudo,ref=l$ref,
                                                       patient =l$patient,
                                                       patient.name = l$patient.name,
                                                       patient.bd =l$patient.bd,
                                                       patient.sex=l$patient.sex)
                                      )))

  if (nrow(ref.info)==0) return(NULL)
  ref.info[] <- data.frame(lapply(ref.info[], as.character), stringsAsFactors=FALSE)
  ref.info <- ref.info [order(ref.info$ref.pseudo), ]
  row.names(ref.info) <- NULL
  
  ref.info_ <- unique(ref.info[,1:2])
  if (nrow(ref.info_) != nrow(ref.info)) warning("Check the uniqueness of the patient : different PID, name, birthday or sex.", call.=FALSE)
  
  
  reg.list <- lapply (reg.list, function (l) l$reg.list)
  reg.list <- reg.list[which(!sapply(reg.list, is.null))]
  names(reg.list) <- sapply (reg.list,function (l) l$object.alias)
  ref.tab <- NULL
  if (length (reg.list)>0){
    ref.tab <- .transfert.ref.tab.create (reg.list)
    ref.tab$dest.ref <-sapply(ref.tab$dest.ref, function(r) ref.info_$ref.pseudo[ref.info_$ref==r])
    ref.tab$src.ref <-sapply(ref.tab$src.ref, function(r) ref.info_$ref.pseudo[ref.info_$ref==r])
  }
  
  reg.info <- list()
  reg.info$patient <- unique(ref.info[,3:6])
  reg.info$file <- data.frame(t=character(0),path = character(0))
  
  if (!is.null(ref.tab)){
    ref.tab <- ref.tab[order (ref.tab$creation.date, decreasing = TRUE),  ]
    dupli.flag <- !duplicated(ref.tab[, 2:3])
    ref.tab <- ref.tab[dupli.flag, ]
    rownames(ref.tab) <- NULL
    ref.tab$t <- paste(ref.tab$dest.ref,ref.tab$src.ref,sep="<-")
    reg.info$file <- ref.tab[ref.tab$dest.ref!=ref.tab$src.ref,c(8,7)]
    ref.tab <- ref.tab[ , 1:5 ]
  }

 
  comb <- unlist(lapply(ref.info_$ref.pseudo,function(r) paste(paste(ref.info_$ref.pseudo,r,sep="<-"),r,ref.info_$ref.pseudo, sep=";")))
  
  if (length(comb) == 0) return (NULL)
  
  comb <- do.call(rbind.data.frame,strsplit(comb,";"))
  names(comb)  <- c("t","src","dest")
  comb$type=""
  comb [] <- data.frame (lapply(comb [], as.character), stringsAsFactors=FALSE)
  list.matrix <-lapply(comb$t,function (r) return(NULL))
  names(list.matrix) <- comb$t
  
  # identite <- matrix (c (rep ( c (1.0, 0.0, 0.0, 0.0, 0.0), 3), 1.0),nrow=4)
  
  # on remplit le tableau
  for (comb.idx in 1:nrow(comb)) {
    if (comb$src[comb.idx] ==  comb$dest[comb.idx]) {
      list.matrix[[comb.idx]] <- diag (4)
      comb$type[comb.idx] <- 'RIGID'
    } else if (!is.null (ref.tab)){
      tab.idx <- which(ref.tab$src.ref==comb$src[comb.idx] & ref.tab$dest.ref==comb$dest[comb.idx] )
      if (length(tab.idx)>0){
        label  <- ref.info_$ref [ref.info_$ref.pseudo==comb$src[comb.idx]]
        comb$type[comb.idx] <- ref.tab$type[tab.idx]
        list.matrix[[comb.idx]] <- (reg.list[[ref.tab$file[tab.idx]]]$ref.data[[label]]$matrix)
        
      }
      tab.idx <- which(ref.tab$dest.ref==comb$src[comb.idx] & ref.tab$src.ref==comb$dest[comb.idx] )
      if (length(tab.idx)>0){
        label  <- ref.info_$ref [ref.info_$ref.pseudo == comb$dest[comb.idx]]
        comb$type[comb.idx] <- ref.tab$type[tab.idx]
        list.matrix[[comb.idx]] <- solve (reg.list[[ref.tab$file[tab.idx]]]$ref.data[[label]]$matrix)
        
      }
    }
  }
  
  change <- TRUE
  while(change){
    idx.left <- which(comb$type=="RIGID" & comb$dest != comb$src)
    change <- FALSE
    for (idx in idx.left){
      new.link <- which(comb$src == comb[idx,]$dest & comb$type=="RIGID" &
                          comb$dest != comb[idx,]$src & comb$dest != comb$src )
      if (length(new.link)>0) for (new.link.idx in 1:length(new.link)){
        to.complete.idx <- which(comb$src == comb[idx,]$src & comb$dest==comb[new.link[new.link.idx],]$dest)
        if (comb[to.complete.idx, ]$type == ""){
          
          comb[to.complete.idx, ]$type <- "RIGID"
          list.matrix[[comb[to.complete.idx, ]$t]] <- list.matrix[[comb[new.link[new.link.idx],]$t]] %*% list.matrix[[comb[idx, ]$t]]
          change <- TRUE
        }
      }
    }
  }
  
  l <- list (ref.info = ref.info_, reg.info=reg.info, matrix.description = comb, matrix.list = list.matrix)
  class(l) <- "t.mat"
  return (l)
}

#############################################################################################
# load fonctions
#############################################################################################

.load.object <- function(Lobj, data=FALSE, nb = NULL, raw.data.list=NULL) {
  modality  <- Lobj$header$modality
  # fname <- file.path (Lobj$header$file.dirname, Lobj$header$file.basename)
  switch( modality,
          "rtstruct" = {
            return (.load.rtstruct (Lobj, nb, raw.data.list))
          },
          "rtimage" = {
            return (.load.vol.object(Lobj, raw.data.list))
          },
          "rtdose" = {
            return (.load.vol.object(Lobj, raw.data.list))
          },
          "ct" = {
            return (.load.vol.object (Lobj))
          },
          "ct1" = {
            return (.load.vol.object (Lobj))
          },
          "mr" = {
            return (.load.vol.object (Lobj))
          },

          "pt" = {
            return (.load.vol.object (Lobj))
          },
          
          "binary" = {
            return (.load.vol.object (Lobj))
          },
          
          "reg" = {
            return (.load.reg (Lobj))
          },
          
          "mesh" = {
            return (.load.mesh (Lobj))
          },
          
          "histo" = {
            return (.load.histo (Lobj))
          },   
          
          "dvh" = {
            return (.load.histo (Lobj))
          }, 
          
          "histo2D" = {
            return (.load.histo2D (Lobj))
          },
          "rtplan" = {
            L <- .load.other (Lobj, raw.data.list)
            class (L) <- "rtplan"
            return (L)
          },
          return (.load.other (Lobj, raw.data.list))
  )
  
}

#############################################################################

.struct.moreinfo <- function(roi.data, ref.from.contour,thickness){
  if (is.null(roi.data)) return (NULL)
  db <- lapply (roi.data , function (l2) {
    if (length(l2)==0) return (list(min.x=NA, max.x=NA,
                                    min.y=NA, max.y = NA,
                                    min.z=NA, max.z = NA,
                                    vol = NA, Gx = NA, Gy = NA, Gz = NA, continue=NA))
    
    pt <- lapply(l2,function(l) l$pt)
    pt <- do.call(rbind,pt)
    pt$t <- 1
    pt <-as.data.frame(t(ref.from.contour %*% t(as.matrix(pt))))


    if (nrow(pt)==0) return (list(min.x=NA, max.x=NA,
                                   min.y=NA, max.y = NA,
                                   min.z=NA, max.z = NA,
                                   vol = NA, Gx = NA, Gy = NA, Gz = NA, continue=NA))
    
    r <- t(matrix (c(range (pt[,1]), range (pt[,2]), range (pt[,3])),
                 ncol=2, byrow = T))
   
    if (sign(thickness)==-1) {
      z.breaks <- sort(unique(unlist (sapply(l2,function(l) l$pt$z))), decreasing = TRUE)
    } else { 
      z.breaks <- sort(unique(unlist (sapply(l2,function(l) l$pt$z))))
    }
    
    if (length(z.breaks)==1)  {
      continue <- TRUE
    } else {
      dz <- unique( round (diff(z.breaks),6))
      continue = FALSE
      if (length(dz)==1)  if (dz==round(thickness,6)) continue <- TRUE
      
    }
    
    

    
    if (setequal (r[1,], r[2,])) return (list (min.x = r[1,1], max.x = r[2,1],
                                               min.y = r[1,2], max.y = r[2,2],
                                               min.z = r[1,3], max.z = r[2,3],
                                               vol = 0,
                                               Gx = r[1,1], Gy = r[1,2], Gz = r[1,3], continue=continue))
    
    dum <- lapply(l2, function (l) {
      # if (length(unique(round(l$pt[ ,3],3)))!=1) return( list(A= NA, Gx=NA, Gy=NA, Gz= NA))
      if (!(castlow.str(l$type) %in% c("closedplanar","point"))) return( list(A= NA, Gx=NA, Gy=NA, Gz= NA))
      if (nrow(l$pt)==1){
        G_ <- as.numeric(c(l$pt[1,],1)) %*% t(ref.from.contour)
        return( list(A= 0, Gx=G_[1], Gy=G_[2], Gz= G_[3]))
      }
      idx <- 2:nrow(l$pt)
      A <- sum (l$pt[idx-1,1] * l$pt[idx,2] - l$pt[idx,1] * l$pt[idx-1,2]) / 2
      G <- c (sum ((l$pt[idx-1,1] + l$pt[idx,1]) * (l$pt[idx-1,1] * l$pt[idx,2] - l$pt[idx,1] * l$pt[idx-1,2])) / (6 * A),
              sum ((l$pt[idx-1,2] + l$pt[idx,2]) * (l$pt[idx-1,1] * l$pt[idx,2] - l$pt[idx,1] * l$pt[idx-1,2])) / (6 * A),
              l$pt[1,3], 1 )
      G_ <- G %*% t(ref.from.contour)
      A <- ifelse (l$level%% 2 == 0, abs(A), -abs(A))
      
      return( list(A= A, Gx=G_[1], Gy=G_[2], Gz= G_[3]))
      
    })
    levels <- sapply (l2, function (l) l$level)
    A <- sapply(dum, function(l) l$A)
    sumA <- sum(A)
    vol <- round (sumA * thickness, 3)
    Gx <- sapply(dum, function(l) l$Gx)
    Gy <- sapply(dum, function(l) l$Gy)
    Gz <- sapply(dum, function(l) l$Gz)
    if (is.na(sumA)){
      G <- c(NA, NA, NA)
    } else if (sumA != 0) {
      G <- c(sum(A * Gx), sum(A * Gy), sum(A * Gz)) / sumA
    } else {
      G <- c(r[1,1],r[1,2], r[1,3])
    }
    return (list (min.x = r[1,1], max.x = r[2,1],
                  min.y = r[1,2], max.y = r[2,2],
                  min.z = r[1,3], max.z = r[2,3],
                  vol = vol/1000, Gx = G[1], Gy = G[2], Gz = G[3], continue=continue))
  })
  n <- names(db[[1]])
  db <- do.call(rbind.data.frame,db)
  db [,1:10] <- data.frame(lapply(db [,1:10], as.numeric), stringsAsFactors=FALSE)
  # colnames(db) <- n
  row.names(db) <- NULL
  return(db)
}

#######################################################################################
.load.rtimage<- function (Lobj, roi.nb = NULL, raw.data.list=NULL) {
}
#######################################################################################
#' @importFrom stats median
#' @importFrom sp point.in.polygon
.load.rtstruct <- function (Lobj, roi.nb = NULL, raw.data.list=NULL) {
  
  # Lobj <- load.Rdcm.raw.data (filename, data=TRUE)
  L <- Lobj$header
  # L$info <- NULL
  class(L) <- "struct"
  
  if (!Lobj$from.dcm)  {
    if (!is.null(Lobj$data)) {
      L$roi.info <- Lobj$data$roi.info
      L$roi.data <- Lobj$data$roi.data
    }
    return (L)
  }
  
  
  Lobj$header <- NULL
  
  if (is.null(raw.data.list)){
    lf <- list.files(L$file.dirname, pattern="[.]Rdcm", full.names = T)
    dicomlist <-lapply (lf,function(f) {
      d <- tryCatch(load.Rdcm.raw.data (f, data=FALSE, address = FALSE), error = function (e) NULL)
      if (length(L$ref.object.info$SOP.ID)!=0){ 
        if(is.null(d$header$object.info$SOP.ID)) return(NULL)
        if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
      if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
      return(d)
    })
  } else {
    dicomlist <- lapply (raw.data.list,function(d) {
      if (length(L$ref.object.info$SOP.ID)!=0){
        if(is.null(d$header$object.info$SOP.ID)) return(NULL)
        if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
      if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
      return(list(header=d$header))
    })
  }
  dicomlist <- dicomlist [which(!sapply (dicomlist, is.null))]
  if (length(dicomlist)==0) dicomlist <- NULL
  ref.object.f <- FALSE
  
  for(Rdcmf in dicomlist){
    if (any(!is.na(match(Rdcmf$header$object.info$SOP.label,L$ref.object.info$SOP.label)))) {
      L$ref.object.alias <- c(L$ref.object.alias, Rdcmf$header$object.alias)
      L$ref.object.alias <- L$ref.object.alias[L$ref.object.alias!=""]
      ref.object.f <- FALSE
      if(!is.null(Rdcmf$header$orientation)) {
        ori <- Rdcmf$header$orientation #unlist(Rdcmf[[1]][[grep("^[(]0020,0037[)]$",names(Rdcmf[[1]]))]])
        if (!any(is.na(ori)))
          L$ref.from.contour <- matrix(c (ori[1:3],0,ori[4:6],0,vector.product(ori[1:3],ori[4:6]),
                                          0,0,0,0,1), ncol=4, byrow = FALSE)
        ref.object.f <- TRUE
      }
      if (!is.null(Rdcmf$header$slice.thickness)) L$thickness <- Rdcmf$header$slice.thickness #  abs(as.numeric(Rdcmf[[1]][[grep("^[(]0018,0050[)]$",names(Rdcmf[[1]]))]]))
      # rm(Rdcmf)
      # break;
    }
    
  }
  
  
  
  rm(dicomlist)
  if (is.null(Lobj$data))  return (L)
  Lobj<- Lobj$data[[1]]
  
  if (L$nb.of.roi>0 & !is.null(Lobj)) {
    
    
    idx.scale <- grep("^[(]3006,0039[)] item", names (Lobj))
    L_ <- lapply( idx.scale, function (i) Lobj[[i]])
    names(L_) <-  gsub("^[(]3006,0039[)]","" ,names (Lobj)[idx.scale])
    rm(Lobj)
    idx.roi <- as.numeric(sapply(names(L_), function(l) unlist(strsplit(l,"[ ]|item"))[3]))
    
    roi.used <- unique(idx.roi)
    if (!is.null (roi.nb)) {
      roi.used <- sort(roi.nb)
      roi.used <- match(roi.used, unique(idx.roi))
      roi.used <- roi.used[!is.na(roi.used)]
    }
    
    if (length(roi.used)==0) return (L)
    
    nb.pt.idx <- grep("[(]3006,0046[)]$", names (L_))
    nb.pt <- as.numeric(unlist(L_[nb.pt.idx]))
    
    type.idx <- grep("[(]3006,0042[)]$", names (L_))
    type <- as.character(unlist(L_[type.idx]))
    
    pt.idx <-  grep("[(]3006,0050[)]$", names (L_))
    pt<- as.character(unlist(L_[pt.idx]))
    idx.max <- which.max(nchar(pt))
    sep.pt <-unique(unlist(strsplit(pt[idx.max],"[[:digit:]]|[ ]|[.]|[e]|[+]|[-]")))
    sep.pt <- sep.pt[sep.pt!=""]
    ptmax <- as.numeric(unlist(strsplit(pt[idx.max],sep.pt, fixed = TRUE))) 
    
    # si on a pas trouvé de ref.object, on crée le référentiel patient, en prenant la RoI qui a le plus données.
    if (!ref.object.f & length(ptmax)>6){
      ptmax <-data.frame (matrix (ptmax, ncol = 3, byrow = TRUE,
                                  dimnames=list(NULL,c("x","y","z"))),
                          stringsAsFactors = FALSE)
      
      if (length(unique(round(ptmax[,3], 3))) != 1){
        vA <- as.numeric(ptmax[round(nrow(ptmax)*0.333), ]-ptmax[1,])
        vB <- as.numeric(ptmax[round(nrow(ptmax)*0.666), ]-ptmax[1,])
        v1 <- round(vA/ as.numeric(sqrt(vA%*%vA)),9)
        v3 <- vector.product(v1, vB)
        v3 <- round(v3/ as.numeric(sqrt(v3%*%v3)),9)
        if (!setequal(v3[1:2],c(0,0))){
          v2 <- round(vector.product(v3, v1),9)
          L$ref.from.contour<- matrix(c(v1,0,v2,0,v3,0,0,0,0,1), ncol=4, byrow = FALSE)
        }
      }
      rm(ptmax) 
    }
    
    M <- t(solve(L$ref.from.contour))
    L$roi.data <-lapply(L$roi.info$name, function(st) return(NULL))
    # names(L$roi.data) <- L$roi.info$name
    for(ru in roi.used){
      idx.roi.used <- which(idx.roi[pt.idx] ==ru)
      if (length(idx.roi.used)>0) L$roi.data[[ru]] <- lapply(idx.roi.used, function(idx){
        
        p <- as.numeric(unlist(strsplit(pt[idx],sep.pt, fixed = TRUE)))
        if (nb.pt[idx] != length(p)/3) {
          warning ("struct loading error : number of points is not a multiple of 3.", call.=FALSE)
          return(NULL)
        }
        if (any(is.na(p)))  return (NULL)
        lpt <- list(type = type[idx])
        lpt$pt <- data.frame (matrix (round ((cbind (matrix (p, byrow=TRUE, ncol=3), 1) %*% M)[, 1:3],6),
                                      byrow=FALSE, ncol=3, dimnames = list(NULL,c ("x","y","z"))),
                              stringsAsFactors = FALSE)
        
        if (castlow.str (lpt$type) == "closedplanar")lpt$pt <- rbind(lpt$pt,lpt$pt[1,])
        
        return(lpt)
      })
    }
    
    #rectification thickness
    if (L$thickness==0){
      nb.plane<- as.numeric(sapply( L$roi.data,function(l) length(l)))
      
      z.diff <- lapply(L$roi.data[which(nb.plane>=2)],function(l) {
        z <- sort(unique(sapply(l,function(l_) median(l_$pt[,3]))))
        if(length(z)<2) return(NULL)
        unique(sort(diff (z)))
      })
      z.diff <- z.diff[!sapply(z.diff,is.null)]
      if (length(z.diff)>0) {
        L$thickness <- unlist(z.diff)
        L$thickness <- tryCatch(round(median(L$thickness [L$thickness >1e-2]),3), error = function (e) 0)
      }
    }
    
    # on vérifie si les contours sont des contours inscrits ou non.
    for (i in roi.used) {
      roi.all.z<- sapply(L$roi.data[[i]], function(li)  li$pt[1,3])
      if (length(roi.all.z)>0) {     
        kz <- rep(0,length(roi.all.z))
        if (L$thickness>0) kz <- round((roi.all.z -roi.all.z[1])/L$thickness)
        
        for (k.value in unique(kz)){
          same.k.roi <- which (kz ==k.value)
          if (length(same.k.roi)>1) {
            for (j in same.k.roi){
              ptj<- L$roi.data[[i]][[j]]$pt
              roi.index.k <-same.k.roi[same.k.roi!=j]
              # if (length(roi.index.z)!=0) {
              r <- unique (sapply (roi.index.k, function (k) {
                ptk <- L$roi.data[[i]][[k]]$pt
                keep <- point.in.polygon (ptj[ ,1], ptj[ ,2],
                                              ptk[ ,1], ptk[ ,2]) > 0.5
                return (ifelse (any(keep), k,NA))}))
              r <- r[!is.na (r)]
              L$roi.data[[i]][[j]]$level <- ifelse (length(r)!=0, length(r), 0)
            } #else L$roi.data[[i]][[j]]$level <- 0
          } else L$roi.data[[i]][[same.k.roi]]$level <- 0
        }
      }
    }
    
    
    db <- .struct.moreinfo (L$roi.data, L$ref.from.contour,L$thickness)
    L$roi.info  <- cbind (L$roi.info, db)
  }
  
  return (L)
}



#######################################################################################
.load.other <- function (Lobj, raw.data.list=NULL) {
  
  # Lobj <- load.Rdcm.raw.data (filename)
  L <- Lobj$header
  if (!is.null (L$ref.object.alias) & (length(L$ref.object.info)!=0)) {
    if (is.null(raw.data.list)){
      lf <- list.files(L$file.dirname, pattern="[.]Rdcm", full.names = T)
      dicomlist <-lapply (lf,function(f) {
        d <- tryCatch(load.Rdcm.raw.data (f, data=FALSE, address = FALSE), error = function (e) NULL)
        if (length(L$ref.object.info$SOP.ID)!=0){
          if(is.null(d$header$object.info$SOP.ID)) return(NULL)
          if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
        if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
        return(d)
      })
    } else {
      dicomlist <- lapply (raw.data.list,function(d) {
        if (length(L$ref.object.info$SOP.ID)!=0){
          if(is.null(d$header$object.info$SOP.ID)) return(NULL)
          if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
        # if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
        return(list(header=d$header))
      })
    }
    
    dicomlist <- dicomlist [which(!sapply (dicomlist, is.null))]
    if (length(dicomlist)==0) dicomlist <- NULL
    
    for(Rdcmf in dicomlist){
      if (any(!is.na(match(Rdcmf$header$object.info$SOP.label, L$ref.object.info$SOP.label)))) {
         L$ref.object.alias <- c(L$ref.object.alias, Rdcmf$header$object.alias)
         L$ref.object.alias <- L$ref.object.alias[L$ref.object.alias!=""]
      }
    }
    
  }
  # L$info <- NULL
  class(L) <- "undef"
  return (L)
}

#######################################################################################
.load.reg <- function (Lobj) {
  
  # Lobj <- load.Rdcm.raw.data (filename)
  L <- Lobj$header
  # L$info <- NULL
  class(L) <- "reg"
  return (L)
}

#######################################################################################

.load.mesh <- function (Lobj) {
  
  L <- Lobj$header
  class(L) <- "mesh"
  if (!is.null(Lobj$data)) {
    L$mesh <- Lobj$data
  }
  return (L)

}
#######################################################################################

.load.histo2D <- function (Lobj) {
  
  L <- Lobj$header
  class(L) <- "histo2D"
  if (!is.null(Lobj$data)) {
    L$density.map <- Lobj$data
  }
  return (L)
  
}

#######################################################################################

.load.histo <- function (Lobj) {
  
  L <- Lobj$header
  class(L) <- "histo"
  if (!is.null(Lobj$data)) {
    L<- c (L,Lobj$data)
  }
  return (L)
  
}

#######################################################################################
.load.vol.object <- function (Lobj, raw.data.list=NULL) {
  
  # Lobj <- load.Rdcm.raw.data (filename, data=TRUE)
  L <- Lobj$header
  if (!is.null (L$ref.object.alias) & (length(L$ref.object.info)!=0)) {
    if (is.null(raw.data.list)){
      lf <- list.files(L$file.dirname, pattern="[.]Rdcm", full.names = T)
      dicomlist <-lapply (lf,function(f) {
        d <- tryCatch(load.Rdcm.raw.data (f, data=FALSE, address = FALSE), error = function (e) NULL)
        if (length(L$ref.object.info$SOP.ID)!=0){
          if(is.null(d$header$object.info$SOP.ID)) return(NULL)
          if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
        # if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
        return(d)
      })
    } else {
      dicomlist <- lapply (raw.data.list,function(d) {
        if (length(L$ref.object.info$SOP.ID)!=0){
          if(is.null(d$header$object.info$SOP.ID)) return(NULL)
          if (!(d$header$object.info$SOP.ID %in% L$ref.object.info$SOP.ID)) return(NULL)}
        # if (d$header$ref.pseudo!=L$ref.pseudo) return(NULL)
        return(list(header=d$header))
      })
    }
    
    dicomlist <- dicomlist [which(!sapply (dicomlist, is.null))]
    if (length(dicomlist)==0) dicomlist <- NULL
	
    
    for(Rdcmf in dicomlist){
      if (any(!is.na(match(Rdcmf$header$object.info$SOP.label, L$ref.object.info$SOP.label)))) {
        L$ref.object.alias <- c(L$ref.object.alias, Rdcmf$header$object.alias)
        L$ref.object.alias <- L$ref.object.alias[L$ref.object.alias!=""]
        # break
      }
    }
  }
  
  # L$info <- NULL
  class(L) <-  "volume"
  Lobj$header <- NULL
  
  if (Lobj$from.dcm){
    if (!is.null(Lobj$data) & is.null(L$error)) {
    
      L$vol3D.data <-unlist(lapply(1:length(Lobj$data), function (i){
        pixel <- tryCatch (as.numeric (Lobj$data[[i]][[grep("[(]7FE0,0010[)]$",names(Lobj$data[[i]]))]]), error = function (e) NA)
        # flag <- pixel > 2^ L$pixeldecode$MSB
        # flag[is.na(flag)] <- FALSE
        # if (L$pixeldecode$signed==1) pixel [flag] <-  pixel [flag] - 2^ L$pixeldecode$bit.stored 
        pixel * L$pixeldecode$slope[i] +    L$pixeldecode$intercept[i]
      }))
      L$vol3D.data <- round(array(L$vol3D.data , dim=c(L$n.ijk[1], L$n.ijk[2], L$n.ijk[3])),3)
    }
    L$pixeldecode <- NULL
    L$dvh <- NULL
  } else {
    if (!is.null(Lobj$data)) L$vol3D.data <- Lobj$data
  }
  if (L$missing.k.idx) 
  warning(paste("some imaging plans of",L$object.alias, 
                "are missing. Some results of espadon calculations are unpredictable. Consider using vol.repair() function to fix it."), call.=FALSE )
  return (L)
}

##########################################################################################
.load.dcm <- function(dcm.filenames, data=FALSE,  ignore.duplicates = FALSE,verbose = TRUE,
                      tag.dictionary = dicom.tag.dictionary ()){
  flag <-file.exists(dcm.filenames) & !dir.exists(dcm.filenames) &
    !grepl("[.]doc$|[.]docx$|[.]dot$|[.]dotx$|[.]odt$|[.]xml$|[.]txt$",dcm.filenames) &
    !grepl("[.]pdf$|[.]ppt$|[.]pptx$|[.]tif$|[.]png$|[.]jpg$|[.]mp3$",dcm.filenames) &
    !grepl("[.]mp4$|[.]wmv$|[.]bmp$|[.]gif$|[.]svg$|[.]exe$|[.]bat$",dcm.filenames) &
    !grepl("[.]hex$|[.]csv$|[.]xls$|[.]xlsx$|[.]qs$|[.]Rdcm$",dcm.filenames)
  dcm.filenames <- dcm.filenames[flag]
  if (length(dcm.filenames)==0) return(NULL)
  dicomlist <- .create.Rdcm.object (dcm.filenames = dcm.filenames, existing.obj = NULL, 
                                    verbose = verbose, update = TRUE, 
                                    save.flag = FALSE, save.dir = NULL, 
                                    only.header = !data,  ignore.duplicates= ignore.duplicates,
                                    tag.dictionary = tag.dictionary)
  
  
  dicomlist <- dicomlist [which(!sapply (dicomlist, is.null))]
  
  return (dicomlist)
}


################################################################################

#' @importFrom methods is
.get.ijkt.from.index <- function (idx, vol) {
  
  .get.ijk.from.index <- function (idx, vol) {
    N.map <- prod (vol$n.ijk[1:2])
    idx <- (idx-1)
    k <- floor (idx/N.map)
    rest <- idx - k*N.map
    j <- floor ((rest/vol$n.ijk[1]))
    i <- rest - j*vol$n.ijk[1]
    t <- rep (1,length(idx))
    k_<- vol$k.idx[k+1]
    matrix (c(i,j,k_,t), byrow=FALSE, ncol=4)
  }
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.", call.=FALSE)
    return (NULL)
  }
  
  if (length(idx)==0) {
    warning ("idx is not defined in vol.", call.=FALSE)
    return (NULL)
  }
  if (is.null(dim(idx))){
    idx_ <- .get.ijk.from.index (idx,vol)
  } else {
    if (dim(idx)[2]!=3) {
      warning ("idx has to be a vector or 3-columns matrix.", call.=FALSE)
      return(NULL)
    } else {
      idx[, 1:2] <- idx[, 1:2]-1
      idx[,3] <- vol$k.idx[idx[,3]]
      idx_ <- cbind(idx,t=1)
    }
  }
  #idx_ <- data.frame(idx_)
  colnames(idx_) <- c ("i","j","k","t")
  #idx_ <- as.matrix (idx_)
  return(idx_)
}


###############################################################################################

.display.select.struct.by.z <- function (struct, list.roi.idx, z, dz){
  new.struct <- lapply(struct$roi.info$roi.pseudo,function (r) return(NULL))
  names(new.struct) <- struct$roi.info$roi.pseudo
  for (roi.idx in list.roi.idx){
    if (length(struct$roi.data[[roi.idx]])!=0){
      cont.z <- round(sapply (struct$roi.data[[roi.idx]], function(co) co$pt$z[1]),3)
      concerned.z <- unique(cont.z)
      near.z <- concerned.z[which.min(abs(concerned.z-z))]
      
      if (near.z >= z-dz/2 && near.z < z+dz/2) {
        c.idx <- which(cont.z==near.z[1])
        new.struct[[roi.idx]] <-  lapply(1:length(c.idx), function (r) return(NULL))
        for (i in 1:length(c.idx)){
          new.struct[[roi.idx]][[i]] <- struct$roi.data[[roi.idx]][[c.idx[i]]]
          colnames (new.struct[[roi.idx]][[i]]$pt) <- c ("x", "y", "z")
        }
        
      }
      # c.idx  <- which (sapply(cont.z, function(zc){zc >= z-dz/2 && zc < z+dz/2}))
      # if (length(c.idx)>0){
      #   new.struct[[roi.idx]] <-  lapply(1:length(c.idx), function (r) return(NULL))
      #   for (i in 1:length(c.idx)){
      #     new.struct[[roi.idx]][[i]] <- struct$roi.data[[roi.idx]][[c.idx[i]]]
      #     colnames (new.struct[[roi.idx]][[i]]$pt) <- c ("abs", "ord", "map")
      #   }
      # }
    }
  }
  return(new.struct)
}
###############################################################################################

#' @importFrom stats ecdf
#' @importFrom grDevices contourLines
.display.roi.data.from.bin <- function (bin){#, pt000=matrix(c(0,0,0), ncol=3)) {
  bin <- .vol.border.tuning(bin, pre.nijk=c(1,1,1), post.nijk=c(1,1,1))
  bin$vol3D.data[is.na(bin$vol3D.data)] <- FALSE
  roi.data <- list()
  # pt000 <- matrix(pt000, ncol=3)
  for(k in 1:length(bin$k.idx)) {
    nb.one <- which(bin$vol3D.data[,,k]==1, arr.ind = TRUE)
    if (nrow(nb.one)==0){
      l <- NULL
    } else if (nrow(nb.one)==1){
      xyz <- as.matrix (cbind (nb.one-1, bin$k.idx[k], 1)) %*% t(bin$xyz.from.ijk)
      nb.one <- as.numeric(nb.one)-1
      pt <- as.data.frame (matrix(c (nb.one[1]-0.4999, nb.one[2]-0.4999, bin$k.idx[k], 1,
             nb.one[1]+0.4999, nb.one[2]-0.4999, bin$k.idx[k], 1,
             nb.one[1]+0.4999, nb.one[2]+0.4999, bin$k.idx[k], 1,
             nb.one[1]-0.4999, nb.one[2]+0.4999, bin$k.idx[k], 1,
             nb.one[1]-0.4999, nb.one[2]-0.4999, bin$k.idx[k], 1), 
             ncol=4,byrow=TRUE) %*% t(bin$xyz.from.ijk))
      colnames(pt) = c("x","y","z","t")
      l <- list ()
      l[[1]] <- list (type="CLOSED_PLANAR", pt=pt[,1:3])
    } else if (nrow(nb.one) > 1){
      xgrid <- 0:(bin$n.ijk[1]-1)
      ygrid <- 0:(bin$n.ijk[2]-1)
      
      l <- contourLines (x = xgrid, y = ygrid, z=bin$vol3D.data[,,bin$k.idx[k]+1], levels =0.5)
      for (i in 1:length(l)) {
        l[[i]]$type <- "CLOSED_PLANAR"
        l[[i]]$level <- NULL
        # l[[i]]$pt <- data.frame (x=bin$dxyz[1]*l[[i]]$x + pt000[1],
        #                          y=bin$dxyz[2]*l[[i]]$y + pt000[2],
        #                          z=bin$dxyz[3]*bin$k.idx[k] + pt000[3])
        
        l[[i]]$pt <- as.matrix(data.frame (x=l[[i]]$x ,y=l[[i]]$y ,z=bin$k.idx[k], t=rep(1,length(l[[i]]$x))))
        l[[i]]$pt  <- l[[i]]$pt %*% t(bin$xyz.from.ijk)
        l[[i]]$pt <- as.data.frame(l[[i]]$pt[ ,1:3])
        colnames (l[[i]]$pt ) <- c("x","y","z")
        l[[i]]$x <- NULL
        l[[i]]$y <- NULL
      }
    }
    roi.data <- c(roi.data,l)
  }
  if (length(roi.data)==0) return (NULL)
  return (roi.data)
}

###############################################################################################

.pixel.scale <- function (min.pixel,max.pixel,scale.nb){
  r <- c(min.pixel,max.pixel)
  if (all(!is.na(r))) {
    if (r[1]!=r[2]) { b <- seq (r[1],r[2], length.out = scale.nb + 1)
    } else { b <- seq (r[1],r[1]+1, length.out = scale.nb + 1)}
  } else if (all(is.na(r))) { b <- NA
  } else if (is.na(r[1])) { b <- seq (r[2]-1,r[2], length.out = scale.nb + 1)
  } else { b <- seq (r[1],r[1]+1, length.out = scale.nb + 1)}
  return (b)
}

###############################################################################################

####################################################################
#' @importFrom graphics abline axis legend lines mtext par points
#' polygon rasterImage rect text
#' @importFrom stats quantile
#' @importFrom grDevices hcl.colors
.display.histo <- function (histo, y, new.ylab=NULL, MC.y=NULL,
                            xgrid = TRUE, ygrid = TRUE,
                            MC.plot=FALSE, MC.col = grey.colors(4, rev = TRUE), add = FALSE, ...) {
  
  
  display.polygon <- function (lb, hb, mids, col, ...) {
    polygon(c(mids, rev(mids)), c(lb, rev (hb)), col=col, ...)
  }
  
  if (!is.null(histo)) {
    args <- list(...)
    
    if (is.null(args[['xlab']])){
      if (castlow.str(histo$mids.unit)==""){
        xlab <- "voxel unit"
      } else {
        xlab <- histo$mids.unit
      }
    } else {xlab <- args[['xlab']]}
    
    
    type <-  ifelse (is.null(args[['type']]), "l", args[['type']])
    
    main <- ifelse (is.null(args[['main']]),  histo$object.alias, args[['main']])
    log <- ifelse (is.null(args[['log']]) , "", tolower(args[['log']]))
    col <- ifelse (is.null(args[['col']]) , "#000000", args[['col']])
    lwd <- ifelse (is.null(args[['lwd']]) , 1, args[['lwd']])
    lty <- ifelse (is.null(args[['lty']]) , 1, args[['lty']])
    
    
    xlim <- c (0,0)
    xlim [1] <-  ifelse (is.null(args[['xlim']]) , min(histo$mids), args[['xlim']][1])
    xlim [2] <-  ifelse (is.null(args[['xlim']]) , max(histo$mids), args[['xlim']][2])
    
    if (grepl("x",log)) {
      x_ <-histo$mids[histo$mids>0 & histo$mids>=xlim[1]& histo$mids <=xlim[2]]
      x.axis <- 10^(floor(log10 (min(x_))):ceiling(log10 (max(x_))))
      xlim <- c(max(xlim[1],min(x.axis)), min(xlim[2],max(x.axis)))
      
    } else {
      dx <-10^round(log10((max(xlim)-min(xlim))/5))
      x.axis <- (floor(min(xlim)/dx): ceiling(max(xlim)/dx))*dx
    }
    
    ylim <- c (0,0)
    ylim [1] <-  ifelse (is.null(args[['ylim']]) , min(y), args[['ylim']][1])
    ylim [2] <-  ifelse (is.null(args[['ylim']]) , max(y), args[['ylim']][2])
    
    if (grepl("y",log)) {
      y_ <-y[y>0 & y>=ylim[1]& y <=ylim[2]]
      y.axis <- 10^(floor(log10 (min(y_))):ceiling(log10 (max(y_))))
      ylim <- c(max(ylim[1],min(y.axis)), min(ylim[2],max(y.axis)))
      
    } else {
      if (max(ylim)!=min(ylim)) {
        dy <-10^round(log10((max(ylim)-min(ylim))/5))
        y.axis <- (floor(min(ylim)/dy): ceiling(max(ylim)/dy))*dy
      } else {
        y.axis <- seq(min(ylim)-1,min(ylim)+1,1)
      }
    }
    
    if (!add){
      plot(xlim,ylim, type= "n",main=main, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, log=log, xlab="",ylab="")
      mtext (xlab, line=2 + 0.5*par()$cex, side =1, cex= par()$cex)
      mtext (new.ylab, line=2 + 0.5*par()$cex, side =2, cex= par()$cex)

      if (grepl("x",log) & xgrid) for (i in 1:9) abline(v=i*x.axis, col="grey80", lty=3)
      if (!grepl("x",log) & xgrid) abline(v=x.axis, col="grey80", lty=3)
      if (grepl("y",log) & ygrid) for (i in 1:9) abline(h=i*y.axis, col="grey80", lty=3)
      if (!grepl("y",log) & ygrid) abline(h=y.axis, col="grey80", lty=3)
      axis(side=1, at=x.axis, labels = x.axis)
      axis(side=2, at=y.axis, labels = y.axis)
    }
    
    if (!is.null(MC.y) & MC.plot) {
      Qmat <- apply (MC.y, 2, function (obs) {quantile (obs, probs = c(0, .025, .25, .5, .75, .975, 1))})

      display.polygon(Qmat[1, ], Qmat[7, ], histo$mids, col= MC.col[1], border = NA)
      display.polygon(Qmat[2, ], Qmat[6, ], histo$mids, col= MC.col[2], border = NA)
      display.polygon(Qmat[3, ], Qmat[5, ], histo$mids, col= MC.col[3], border = NA)
      lines (histo$mids, Qmat[4, ], col=MC.col[4])
    }
    if ((type=="h") | type=="l") lines (histo$mids, y, col=col, type=type, lwd=lwd, lty=lty)
    if ((type=="p") ) points (histo$mids, y)
  }
}


#####################################################################################################
.kernel <- function(Vb,radius){
  cut.v <- function (v) {c(v[-(1:(length(v)/2+0.5))], v[1:(length(v)/2+0.5)])}
  i.c <- cut.v (1:Vb$n.ijk[1] -1)
  j.c <- cut.v (1:Vb$n.ijk[2] -1)
  k.c <- cut.v (Vb$k.idx)
  
  gr <- as.matrix (expand.grid (i.c, j.c, k.c, 1)) %*% t(Vb$xyz.from.ijk)
  gc <- (c(mean(i.c) + 1, mean(j.c) + 1, mean(k.c) +1, 1) %*% t(Vb$xyz.from.ijk)) [1,1:3]
  idx <- which ((gr[,1] - gc[1])^2 + (gr[,2] - gc[2])^2 + (gr[,3] - gc[3])^2 <= radius^2)
  
  kernel <- Vb$vol3D.data
  kernel [,,] <- FALSE
  kernel [idx] <- TRUE
  return( kernel)

}


################################################################################
#' @importFrom Rvcg setRays vcgRaySearch nfaces
# @export 
.mesh.pt.distance <- function(mesh, xyz.pt, unit.vector) {
  .get.pt.flag<- function(xyz.c, xyz.pt= c(0,0,0), xyz.test = c(">=",">=",">=")) {
    .test.decript <- function(vect, value=0, test=">="){
      t <- as.character(match(test,c(">=", ">","<=","<")))
      if (is.na(t[1])) return(rep(FALSE,length(vect)))
      switch (t,
              "1" = return(vect>=value),
              "2" = {return(vect>value)},
              "3" = {return(vect<=value)},
              "4" = {return(vect<value)}
      )
    }
    
    rf <- rep(TRUE,nrow(xyz.c))
    for (ci in 1:ncol(xyz.c)) rf <- rf & .test.decript(xyz.c[,ci], test = xyz.test[ci], 
                                                       value = xyz.pt[ci])
    # xyz.c[rf,]
    rf
  }
  .get.mesh.selection <- function (mesh, xyz.pt= c(0,0,0), xyz.test = c(">=",">=",">=")){
    
    vb.idx <- which(.get.pt.flag(xyz.c=t(mesh$vb)[,1:3],xyz.pt=xyz.pt, xyz.test = xyz.test))
    #on cherche les triangles impactés
    f1 <- !is.na(match(mesh$it[1,],vb.idx)) | !is.na(match(mesh$it[2,],vb.idx)) | 
      !is.na(match(mesh$it[3,],vb.idx))
    pt.idx1 <- sort(unique(as.numeric(mesh$it[,f1])))
    
    m1 <-mesh
    m1$vb <- mesh$vb[,pt.idx1]
    m1$it <- mesh$it[,f1]
    m1$it[1,] <- match(mesh$it[1,f1],pt.idx1)
    m1$it[2,] <- match(mesh$it[2,f1],pt.idx1)
    m1$it[3,] <- match(mesh$it[3,f1],pt.idx1)
    m1$normals <- mesh$normals[,pt.idx1]
    m1$remface <- mesh$remface[f1]
    return(m1)
  }
  
  mesh$vb <- mesh$vb-t(matrix(rep(c(xyz.pt,0),ncol(mesh$vb)), 
                              ncol=4, byrow=TRUE))
  
  mesh$vb[abs(mesh$vb)<1e-6 ]  <- 0
  
  test.matrix.uv <- expand.grid( c(">=","<"), c(">=","<"), c(">=","<"))
  test.matrix.uv[] <- lapply(test.matrix.uv[] ,as.character)
  test.matrix <- expand.grid( c(">=","<="), c(">=","<="), c(">=","<="))
  sign.matrix <- expand.grid( c(1,-1), c(1,-1), c(1,-1))
  test.matrix[] <- lapply(test.matrix[] ,as.character)
  sign.matrix[] <- lapply(sign.matrix[] ,as.numeric)
  
  dist <- rep(NA,nrow(unit.vector))
  for (zone.idx in 1:nrow(test.matrix)){
    xyz.test.uv <- as.character(test.matrix.uv[zone.idx,])
    xyz.test <- as.character(test.matrix[zone.idx,])
    f.zone <- .get.pt.flag(unit.vector, xyz.test = xyz.test.uv)
    s <- as.numeric(sign.matrix[zone.idx,])
    m <- .get.mesh.selection(mesh, xyz.pt= rep(-1,3)*s,xyz.test = xyz.test)
    if (nfaces(m)>0){
      ray <- setRays (matrix(rep(c(0,0,0),nrow(unit.vector[f.zone,])), ncol=3, byrow= T), 
                      as.matrix(unit.vector[f.zone,]))
      result <- vcgRaySearch(ray, m, mindist = TRUE)
      fq <- result$quality==1
      dist[f.zone][fq] <- result$distance[fq]
    
      # open3d()
      # par3d(windowRect=wr)
      # wire3d(m)
      # na.pt <- 50* (unit.vector[f.zone,])[!fq,]
      # segments3d( matrix(t(cbind(matrix(c(0,0,0),nrow(na.pt),byrow=T,ncol=3), na.pt)),
      #                    byrow=T,ncol=3),col="red")
      # segments3d( matrix(t(cbind(matrix(c(0,0,0),nrow(unit.vector[f.zone,]),byrow=T,ncol=3),
      #                            dist[f.zone] *unit.vector[f.zone,])),
      #                    byrow=T,ncol=3),col="red")
      
    }
  }
  return(dist)
}

# @importFrom rgl mesh3d
# @export
# .sphere.mesh <- function(center=c(0,0,0),rayon=1, iteration =7){
#   dist<- function(vect) {sqrt(sum(vect^2))}
#   
#   pt <- rayon*rbind(diag(3),-diag(3))
#   pt <- pt[order(pt[,1],pt[,2]),]
#   row.names(pt) <- NULL
#   it <- expand.grid(1:nrow(pt),2:nrow(pt),3:nrow(pt))
#   it <- it[sapply(1:nrow(it), function (r.idx) length(unique(match(it[r.idx,],it[r.idx,])))==3),]
#   it <- unique(do.call(rbind,lapply(1:nrow(it), function (r.idx) sort(as.numeric(it[r.idx,])))))
#   AB <- apply(pt[it[,2],]-pt[it[,1],],1,dist)
#   AC <- apply(pt[it[,3],]-pt[it[,1],],1,dist)
#   BC <- apply(pt[it[,3],]-pt[it[,2],],1,dist)
#   l <- min(c(AB,BC,AC))
#   it <- it [which(AB==l & AC==l & BC==l), ]
#   
#   # pt <- matrix(c(sqrt(2/3),0,-sqrt(1/3), 
#   #                -sqrt(2/3),0,-sqrt(1/3),
#   #                0,sqrt(2/3), sqrt(1/3),
#   #                0,-sqrt(2/3), sqrt(1/3)), ncol=3 ,byrow = T)
#   # 
#   # it <- matrix(c(1,2,3, 1,2,4,1,3,4,2,3,4), ncol=3 ,byrow = T)
#   
#   row.names(it) <- NULL
#   row.names(pt) <- NULL
#   
#   iterat <- 1
#   while (iterat<iteration) {
#     
#     mid12 <- (pt[it[,1],]+ pt[it[,2],])/2
#     mid23 <- (pt[it[,2],]+ pt[it[,3],])/2
#     mid13 <- (pt[it[,1],]+ pt[it[,3],])/2
#     mid <- unique(rbind(mid12,mid23,mid13))
#     txtmid <- apply(mid,1,paste,collapse=";")
#     vb12 <- match(apply(mid12,1,paste,collapse=";"),txtmid)+ nrow(pt)
#     vb23 <- match(apply(mid23,1,paste,collapse=";"),txtmid)+ nrow(pt)
#     vb13 <- match(apply(mid13,1,paste,collapse=";"),txtmid)+ nrow(pt)
#     pt <- rbind(pt,rayon*mid/apply(mid,1,dist))
#     row.names(pt) <- NULL
#     it <- do.call(rbind,lapply(1:nrow(it), function(r.idx) {
#       matrix(c(it[r.idx,1],vb12[r.idx],vb13[r.idx],
#                vb12[r.idx],it[r.idx,2],vb23[r.idx],
#                vb12[r.idx],vb23[r.idx],vb13[r.idx],
#                vb13[r.idx],vb23[r.idx],it[r.idx,3]), ncol=3, byrow=T)
#     }))
#     iterat <- iterat+1
#   }
#   
#   pt <- sweep(pt,2,-center)
#  m <-  mesh3d(vertices=t(pt), triangles=t(it))
#  m
# }

.get.tag.info <- function(db.name, db.tag, num.tag = numeric(), post.tag="",
                          data,  data.name=names(data), adjust = TRUE){
  if (length(post.tag)==0) return (NULL)
  db <- array("", dim=c(length(post.tag),length(db.name)))
  colnames(db) <- db.name
  rc.db <- expand.grid(1:length(post.tag),1:length(db.name))
  tag.list <- trimws(apply(expand.grid(post.tag, db.tag), 1, paste, collapse=" "))
  ma <- match(tag.list, data.name)
  not.na <- which(!is.na(ma))
  if (length(not.na)==0) return(NULL)
  db[as.matrix(rc.db[not.na,])] <- unlist(trimws(data[ma[not.na]]))
  db <- data.frame(trimws(db),stringsAsFactors = F)
  if (length(num.tag)!=1) {
    db[ ,num.tag] <- suppressWarnings(lapply(db[,num.tag], as.numeric))
  } else {
    if ((length(num.tag)==1)) db[ ,num.tag] <- as.numeric(db[ ,num.tag])
  }
  if (adjust) db <- .reduc.tab(db)
  db
}

.reduc.tab <- function (db) { 
  if (is.null(db)) return (NULL)
  f <- sapply(colnames(db),function(c) any(db[,c]!="") & any(!is.na(db[,c])))
  if (sum(f)==0 | nrow(db)==0) return (NULL)
  n <- colnames(db)[f]
  rn <- NULL
  if (nrow(db)>0) rn <- row.names(db)
  if (sum(f)==1) {db <- data.frame(n=db[,f]) 
  } else {db <- db[,f]} 
 
  if (!is.null(rn)){
    f <- sapply(row.names(db),function(c) any(db[c,]!="") & any(!is.na(db[c,])))
    if(ncol(db)==1){ db <- data.frame(n=db[f,]) 
    } else { db <- db[f,]}
    row.names(db) <- rn[f]
  }
  colnames(db) <- n
  db
}

.fill.db <- function(db) {
  n <- colnames(db)
  num.tag <- which(as.logical(sapply(db,is.numeric)))
  nr <- nrow(db)
  db_ <- data.frame(do.call(cbind,lapply(1: ncol(db),function(col.nb){
    idx <-c(which(!(is.na(db[,col.nb]) | db[,col.nb]=="")), nr+1)
    if (length(idx)==1) return(db[,col.nb]) 
    d.idx <- diff(idx)
    postdb <- NULL
    if (idx[1]>1) postdb <- db[1:(idx[1]-1),col.nb]
    return (c(postdb, 
              db[unlist(lapply(1:length(d.idx),function (i) rep(idx[i],d.idx[i]))),col.nb]))
  })))
  colnames (db_) <- n
  if (length(num.tag)>1) {
    db_[ ,num.tag] <- suppressWarnings(lapply(db_[,num.tag], as.numeric))
  } else {
    if (length(num.tag)==1) db_[,num.tag] <- as.numeric(db_[,num.tag])
  }
  db_
}

