#' Loading patient data from *.Rdcm files
#' @description The \code{load.patient.from.Rdcm} function is used to load or 
#' pre-load in memory all patient objects converted in *.Rdcm files.
#' @param dirname Full paths of the directories of a single patient.
#' @param data Boolean. If \code{data = TRUE}, the voxels value of the "volume" 
#' class objects, or the coordinates of the RoI (region of interest)
#' of the \code{struct} class objects, are loaded into memory.
#' @param dvh Boolean. if \code{dvh = TRUE} and if they exist, patient DVH are 
#' loaded, for convenience. They are not used as is in \pkg{espadon} package.
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
load.patient.from.Rdcm <- function (dirname, data = FALSE, dvh = FALSE){
  
  lf <- list.files (dirname, pattern = "[.]Rdcm", full.names =  TRUE) 
  if (length(lf)==0) {
    warning ("no patient found \n")
    return(NULL)
  }
  
  dicomlist <-lapply (lf,function(f) {
    d <- tryCatch(load.Rdcm.raw.data (f, data=data, address=FALSE), error = function (e) NULL)
    return(d)
  })
  
  
  dicomlist <- dicomlist [which(!sapply (dicomlist, is.null))]
  names(dicomlist) <- sapply(dicomlist, function(l) l$header$object.alias)
  
  # base.n <- do.call(rbind.data.frame, lapply(dicomlist, function(l) {
  #   nb <- switch(l$header$modality, "rtstruct" =  l$header$nb.of.roi,  
  #                "rtdose" = l$header$n.ijk[3], "ct" = l$header$n.ijk[3],
  #                "ct1" = l$header$n.ijk[3], "mr" = l$header$n.ijk[3], 
  #                "pt" = l$header$n.ijk[3], "binary" = l$header$n.ijk[3],
  #                "reg" = l$header$nb.of.ref,  NA)
  #   c(l$header$patient,as.character(l$header$patient.bd),l$header$patient.sex, l$header$modality, l$header$object.name, l$header$ref.pseudo,
  #     tryCatch(l$header$object.info[[grep ("nb[.]of[.]subobj", names(l$header$object.info))]], error = function(e) NA),
  #     l$header$description, nb,
  #     tryCatch(l$header[[grep ("max[.]pixel", names(l$header))]], error = function(e) NA),
  #     l$header$object.alias, l$header$file.basename)
  # }))
  
  base.n <- do.call(rbind.data.frame, lapply(dicomlist, function(l) {
    nb <- switch(l$header$modality, "rtstruct" =  l$header$nb.of.roi,  
                 "rtdose" = l$header$n.ijk[3], "ct" = l$header$n.ijk[3],
                 "ct1" = l$header$n.ijk[3], "mr" = l$header$n.ijk[3], 
                 "pt" = l$header$n.ijk[3], "binary" = l$header$n.ijk[3],
                 "reg" = l$header$nb.of.ref, "mesh" = l$header$nb.faces, 
                 "histo" = l$header$nb.MC,"dvh" = l$header$nb.MC,
                 "histo2D" = l$header$nb.pixels, 
				 "rtplan"=sum(l$header$fraction.info[1,c("beam.nb","brachy.app.nb")],na.rm=T), 
				 NA)
    subobj.nb <- NA
    idx <- grep ("nb[.]of[.]subobj", names(l$header$object.info))
    if (length(idx)>0) subobj.nb <-l$header$object.info[[idx]]
    max.pix <- NA
    idx <- grep ("max[.]pixel", names(l$header))
    if (length(idx)>0)  max.pix <- l$header[[idx]]
    c(l$header$patient,as.character(l$header$patient.bd),l$header$patient.sex, l$header$modality, l$header$object.name, l$header$ref.pseudo,
      subobj.nb, l$header$description, nb,
      max.pix,
      l$header$object.alias, l$header$file.basename)
  }))
  
  colnames(base.n) <- c ("PIN", "birth.date","sex", "modality", "obj", "ref.pseudo", "nb.of.subobject" ,"description", "nb", "max","object.alias", "file.basename")
  base.n <- base.n[order(base.n$PIN,base.n$ref.pseudo,base.n$modality),]
  base.n$max<- suppressWarnings(as.character(round(as.numeric(base.n$max),3)))
  base.n$nb<- suppressWarnings(as.numeric(base.n$nb))
  row.names(base.n) <- NULL
  
  l <- list()
  l$patient <- unique (base.n[,c ("PIN", "birth.date", "sex")])
  l$patient$birth.date <- l$patient$birth.date
  row.names(l$patient) <- NULL
  l$pat.pseudo <- l$patient[1,1]
  
  # db <- base.n[base.n[,1]==patient[patient.idx],2:5 ]
  
  l$description <- unique(base.n[,c (1, 4:8)])
  row.names (l$description) <- NULL
  l$description$nb <- sapply (l$description$obj, function (obj) paste(base.n$nb[which(base.n$obj==obj)], collapse = ";"))
  l$description$max <- sapply (l$description$obj, function (obj) paste(base.n$max[which(base.n$obj==obj)], collapse = ";"))
  l$description$object.alias <- sapply (l$description$obj, function (obj) paste(base.n$object.alias[which(base.n$obj==obj)], collapse = ";"))
  
  
  l$description.by.reg <- list ()
  
  l$T.MAT <- load.T.MAT (dirname)
  
  modality <- sort (unique (base.n$modality[base.n$modality!="reg"]))
  obj <- lapply(modality, function (m) {
    obj.flag <- base.n$modality==m
    alias <-  base.n$object.alias[obj.flag]
    fname <- file.path(dirname,base.n$file.basename[obj.flag])
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
  
  l$dicom.dvh <- lapply(dicomlist, function(obj) return(obj$header$dvh))
  l$dicom.dvh <- l$dicom.dvh [which(!sapply (l$dicom.dvh, is.null))]
  
  if (length( l$dicom.dvh)>0){
    for (dvh.idx in 1:length(l$dicom.dvh)){
      ref.object <- names(dicomlist)[(sapply(dicomlist,function(dl) dl$header$object.info$SOP.ID)==l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.ID) & 
                                       (sapply(lapply(dicomlist, function(dl) dl$header$object.info$SOP.label), function(v) l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.label %in% v))]
      if (length(ref.object)>0) l$dicom.dvh[[dvh.idx]]$ref.object.name <- ref.object
      
      if (l$dicom.dvh[[dvh.idx]]$ref.object.name!=""){
        ma <- match(l$dicom.dvh[[dvh.idx]]$info$number, dicomlist[[l$dicom.dvh[[dvh.idx]]$ref.object.name]]$header$roi.info$number)
        l$dicom.dvh[[dvh.idx]]$info$roi.name[!is.na(ma)] <- 
          dicomlist[[l$dicom.dvh[[dvh.idx]]$ref.object.name]]$header$roi.info$name[ma[!is.na(ma)]]
        # n <- names(l$dicom.dvh[[dvh.idx]]$data)
        # n[!is.na(ma)] <- l$dicom.dvh[[dvh.idx]]$info$roi.name[!is.na(ma)]
        # names(l$dicom.dvh[[dvh.idx]]$data) <- n
      }
    }
  } else {l$dicom.dvh <- NULL}
  
  
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

