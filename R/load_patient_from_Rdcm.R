#' Loading patient data from *.Rdcm files
#' @description The \code{load.patient.from.Rdcm} function is used to load or 
#' pre-load in memory all patient objects converted in *.Rdcm files.
#' @param dirname Full paths of the directories of a single patient, or vector 
#' of full.path of Rdcm.files.
#' @param data Boolean. If \code{data = TRUE}, the voxels value of the "volume" 
#' class objects, or the coordinates of the RoI (region of interest)
#' of the \code{struct} class objects, are loaded into memory.
#' @param dvh Boolean. if \code{dvh = TRUE} and if they exist, patient DVH are 
#' loaded, for convenience. They are not used as is in \pkg{espadon} package.
#' @param upgrade.to.latest.version Boolean. If \code{TRUE}, the function attempts 
#' to upgrade to the latest version, parsing the DICOM data. It may take longer 
#' to load the data. Consider using the\link[espadon]{Rdcm.upgrade} function.
#' @param ignore.duplicates Boolean. If \code{TRUE}, the function ignores duplicated objects.
#' @return Returns an \pkg{espadon} object of class "patient", describing the 
#' information contained in \code{dirname}. See \link[espadon]{espadon.class} for a 
#' description of the "patient" class.

#' @seealso \link[espadon]{dicom.to.Rdcm.converter}, \link[espadon]{load.patient.from.dicom}, 
#' \link[espadon]{load.obj.data}, \link[espadon]{load.obj.from.dicom}, 
#' \link[espadon]{load.obj.from.Rdcm} and \link[espadon]{load.T.MAT}.
#' @examples
#' # First, save toy patient objects to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = c (4, 4, 4))
#' save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
#' save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
#' save.T.MAT (patient$T.MAT, dirname = pat.dir)
#' # Rdcm files in pat.dir
#' list.files(pat.dir)
#' 
#' # loading patient from Rdcm files with data: 
#' new.patient <- load.patient.from.Rdcm (pat.dir, data = TRUE)
#' str (new.patient, max.level = 2 )
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @export
load.patient.from.Rdcm <- function (dirname, data = FALSE, dvh = FALSE,
                                    upgrade.to.latest.version = FALSE,
                                    ignore.duplicates = FALSE){
  
  
  if (length(dirname)==0) {
    warning ("no files to load.")
    return (NULL)
  } 
  flag <- dir.exists(dirname)
  Rdcm.dir <- dirname[flag]
  dcm.filenames1 <-  list.files(Rdcm.dir,pattern = "[.]Rdcm",recursive = TRUE,full.names = TRUE)
  dcm.filenames2 <- dirname[!flag]
  dcm.filenames2 <- dcm.filenames2[grepl("[.]Rdcm$",dcm.filenames2)]
  dcm.filenames2 <- dcm.filenames2[file.exists(dcm.filenames2)]
  lf <- c(dcm.filenames1,dcm.filenames2)
  
  if (length(lf)==0) {
    warning ("no patient found \n")
    return(NULL)
  }
  
  dicomlist <-lapply (lf,function(f) {
    d <- tryCatch(load.Rdcm.raw.data (f, data=data, address=FALSE,
                                      upgrade.to.latest.version = upgrade.to.latest.version),
                  error = function (e) list(stopmessage= e$message))
    if  (!is.null(d$stopmessage))  {
      warning(d$stopmessage)
      return(NULL)
    }         
    return(d)
  })
  
  f <- !sapply(dicomlist, is.null)
  dicomlist <- dicomlist[f]
  if (length(dicomlist) == 0) return(NULL)
  lf <- lf[f]
  
  ok <- sapply(dicomlist, function (l) {
    if (is.null(l)) return(FALSE)
    return (is.null(l$header$error)) 
  })
  
  sort.idx <- order(ok, decreasing = TRUE)
  
  lf <- lf[sort.idx]
  ok <- ok[sort.idx] 
  dicomlist <- dicomlist[sort.idx]
  SOP <- sapply(dicomlist, function (l) {
    SOP_ <- l$header$object.info$SOP.label
    return(paste0(sort(SOP_), collapse = ";"))
  })

  exist.SOP <- SOP!=""
  SOP.f <- rep(FALSE,length(SOP))
  SOP.f[exist.SOP] <- duplicated(SOP[exist.SOP])
  if (any(SOP.f)){
    if (ignore.duplicates){
      lf <- lf[!SOP.f]
      dicomlist <-lapply (lf,function(f) {
        d <- tryCatch(load.Rdcm.raw.data (f, data=data, address=FALSE,
                                          upgrade.to.latest.version = upgrade.to.latest.version),
                      error = function (e) NULL)
        return(d)
      })
    } else {
      warning("Some objects are duplicated. Consider ignore.duplicates = TRUE", call. = FALSE)
    }
  }
  
  update.list <- sapply(dicomlist, function(l) l$update.needed)
  
  if (any(update.list)){
    if (upgrade.to.latest.version) {
      warning("Patient file versions have been upgraded, consider using Rdcm.upgrade() for faster loading.")
    } else {
      warning("Patient file versions are not up to date, consider using Rdcm.upgrade().")
    }
  }
    
  names(dicomlist) <- sapply(dicomlist, function(l) l$header$object.alias)
  
  base.n <- do.call(rbind.data.frame, lapply(dicomlist, function(l) {
    nb <- switch(l$header$modality, "rtstruct" =  l$header$nb.of.roi,  
                 "rtdose" = l$header$n.ijk[3], "ct" = l$header$n.ijk[3],
                 "ct1" = l$header$n.ijk[3], "mr" = l$header$n.ijk[3], 
                 "pt" = l$header$n.ijk[3], "binary" = l$header$n.ijk[3],
                 "reg" = l$header$nb.of.ref, "mesh" = l$header$nb.faces, 
                 "histo" = l$header$nb.MC,"dvh" = l$header$nb.MC,
                 "histo2D" = l$header$nb.pixels, 
				 "rtplan"=sum(l$header$fraction.info[1,c("nb.of.beam","nb.of.brachy.app")],na.rm=T), 
				 NA)
    subobj.nb <- NA
    idx <- grep ("nb[.]of[.]subobj", names(l$header$object.info))
    if (length(idx)>0) subobj.nb <-l$header$object.info[[idx]]
    max.pix <- NA
    idx <- grep ("max[.]pixel", names(l$header))
    if (length(idx)>0)  max.pix <- l$header[[idx]]
    c(l$header$patient, l$header$patient.name, as.character(l$header$patient.bd),
      l$header$patient.sex, l$header$modality, l$header$object.name, l$header$ref.pseudo,
      subobj.nb, l$header$description, nb,
      max.pix,
      l$header$object.alias, l$header$file.basename)
  }))
  
  colnames(base.n) <- c ("PIN", "name","birth.date","sex", "modality", "obj", 
                         "ref.pseudo", "nb.of.subobject" ,"description", "nb", 
                         "max","object.alias", "file.basename")
  base.n <- base.n[order(base.n$PIN,base.n$ref.pseudo,base.n$modality),]
  base.n$max<- suppressWarnings(as.character(round(as.numeric(base.n$max),3)))
  base.n$nb<- suppressWarnings(as.numeric(base.n$nb))
  row.names(base.n) <- NULL
  
  l <- list()
  l$patient <- unique (base.n[,c ("PIN", "name", "birth.date", "sex")])
  row.names(l$patient) <- NULL
  if (nrow(l$patient) != 1) 
    warning("Check the uniqueness of the patient : different PID, name, birthday or sex.", call.=FALSE)
  
  l$pat.pseudo <- l$patient[1,1]
  
  # db <- base.n[base.n[,1]==patient[patient.idx],2:5 ]
  
  l$description <- unique(base.n[,c (1, 5:9)])
  row.names (l$description) <- NULL
  l$description$nb <- sapply (l$description$obj, function (obj) paste(base.n$nb[which(base.n$obj==obj)], collapse = ";"))
  l$description$max <- sapply (l$description$obj, function (obj) paste(base.n$max[which(base.n$obj==obj)], collapse = ";"))
  l$description$object.alias <- sapply (l$description$obj, function (obj) paste(base.n$object.alias[which(base.n$obj==obj)], collapse = ";"))
  
  
  l$description.by.reg <- list ()
  
  l$T.MAT <- suppressWarnings(load.T.MAT (lf))
  
  modality <- sort (unique (base.n$modality[base.n$modality!="reg"]))
  obj <- lapply(modality, function (m) {
    obj.flag <- base.n$modality==m
    
    alias <-  base.n$object.alias[obj.flag]
    tab.L <- strsplit(alias,paste0("[_]ref|[_]do|[_]", m))
    tab.L <- do.call(rbind.data.frame,lapply(tab.L,function(v) c(v,rep("",4))[1:4]))
    tab.L[2:4] <-lapply(tab.L[2:4], as.numeric)
    alias <- alias[order(tab.L[,1],tab.L[,2],tab.L[,3],tab.L[,4])]
    
    # fname <- file.path(dirname,base.n$file.basename[obj.flag])
    match.index <- match(alias, sapply (dicomlist, function(l) l$header$object.alias))
    lobj.l <- list()
    for (idx in 1:length(alias)) lobj.l[[idx]] <- .load.object (Lobj = dicomlist[[match.index[idx]]], data=data, raw.data.list=dicomlist)
    # lobj.l <- lapply(1:length(alias), function(idx) 
    # .load.object (dicomlist[[match.index[idx]]], fname[idx], data=data))
    names(lobj.l) <- alias
    return(lobj.l)
  })
  names(obj) <- modality
  l <- do.call(c, list(l, obj))
  if (dvh){
    l$dicom.dvh <- lapply(dicomlist, function(obj) return(obj$header$dvh))
    l$dicom.dvh <- l$dicom.dvh [which(!sapply (l$dicom.dvh, is.null))]
    
    if (length( l$dicom.dvh)>0){
      for (dvh.idx in 1:length(l$dicom.dvh)){
        ref.object <- names(dicomlist)[(sapply(dicomlist,function(dl) dl$header$object.info$SOP.ID)==l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.ID) & 
                                         (sapply(lapply(dicomlist, function(dl) dl$header$object.info$SOP.label), function(v) l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.label %in% v))]
        if (length(ref.object)>0) l$dicom.dvh[[dvh.idx]]$ref.object.alias <- ref.object
        
        if (l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]!=""){
          ma <- match(l$dicom.dvh[[dvh.idx]]$info$number, dicomlist[[l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]]]$header$roi.info$number)
          l$dicom.dvh[[dvh.idx]]$info$roi.name[!is.na(ma)] <- 
            dicomlist[[l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]]]$header$roi.info$name[ma[!is.na(ma)]]
          # n <- names(l$dicom.dvh[[dvh.idx]]$data)
          # n[!is.na(ma)] <- l$dicom.dvh[[dvh.idx]]$info$roi.name[!is.na(ma)]
          # names(l$dicom.dvh[[dvh.idx]]$data) <- n
        }
      }
    } else {l$dicom.dvh <- NULL}
  }
  
  l.reg <- strsplit(names(l$T.MAT$matrix.list)[sapply(l$T.MAT$matrix.list, function(m) !is.null(m))],"<-")
  l.reg <- lapply(l.reg,function(li) sort(unique(li)))
  l.reg <- l.reg[!duplicated(l.reg)]
  l.reg_ <- list ()
  reg.idx <- 1
  for (ref in sort(unique(unlist(l.reg)))) {
    v <- sort(unique(unlist(l.reg[sapply(l.reg,function(l) ref %in% l)])))
    if (!is.null(v)){
      l.reg_ [[reg.idx]]<- v
      l.reg <- l.reg[!sapply(l.reg, function(li) any(!is.na(match(li,l.reg_[[reg.idx]]))))]
      if (length(l.reg)==0) break
      reg.idx <- reg.idx + 1
    }
  }
  l$description.by.reg <- lapply (l.reg_, function (li) l$description[l$description$ref.pseudo %in% li,])
  
  return (l)
}

