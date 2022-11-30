#' Loading patient data from DICOM files
#' @description The \code{load.patient.from.dicom} function is used to load or 
#' pre-load in memory all patient objects from DICOM files.
#' @param dcm.files String vector, representing the list of the full names of the  
#' DICOM files of the same patient, or its directories.
#' @param data Boolean. If \code{data = TRUE}, the voxels value of the "volume" 
#' class objects, or the coordinates of the RoI (region of interest)
#' of the \code{struct} class objects, are loaded into memory.
#' @param dvh Boolean. if \code{dvh = TRUE} and if they exist, patient DVH are 
#' loaded, for convenience. They are not used as is in \pkg{espadon} package.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @param verbose Boolean. If \code{TRUE}, a progress bar indicates the progress 
#' of the conversion.
#' @return Returns an \pkg{espadon} object of class "patient", describing the 
#' information from \code{dcm.files}. See \link[espadon]{espadon.class} for a 
#' description of the "patient" class.

#' @seealso  \link[espadon]{dicom.to.Rdcm.converter}, \link[espadon]{load.patient.from.Rdcm}, 
#' \link[espadon]{load.obj.data},  \link[espadon]{load.obj.from.dicom}, 
#' \link[espadon]{load.obj.from.Rdcm} and \link[espadon]{load.T.MAT}.
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "toy_dccm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = pat.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' 
#' # loading patient. Here the toy patient ha only a unique rt-plan object
#' patient <- load.patient.from.dicom (pat.dir, data = FALSE)
#' str (patient, max = 2)
#' # description of object
#' patient$description
#' # transfer matrices :
#' patient$T.MAT
#' # rt-plan object
#' str (patient$rtplan[[1]])

#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @export
load.patient.from.dicom <- function (dcm.files, data = FALSE, dvh = FALSE, 
                                     tag.dictionary = dicom.tag.dictionary (),
                                     verbose = TRUE){
  
  if (length(dcm.files)==0) {
    warning ("no files to load.")
    return (NULL)
  } 
  flag <- dir.exists(dcm.files)
  dcm.dir <- dcm.files[flag]
  dcm.filenames1 <-  list.files(dcm.dir,recursive = TRUE,full.names = TRUE)
  dcm.filenames2 <- dcm.files[!flag]
  dcm.filenames2 <- dcm.filenames2[file.exists(dcm.filenames2)]
  dcm.filenames <- c(dcm.filenames1,dcm.filenames2)
  if (length(dcm.filenames)==0) {
    warning ("no files to load.")
    return (NULL)
  } 
  obj.list <- .load.dcm (dcm.filenames, data=data, tag.dictionary = tag.dictionary,
                         verbose = verbose)
  
  if (is.null(obj.list)) return (NULL)
  
  base.n <- do.call (rbind.data.frame, lapply(obj.list, function(l) {
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
    c(l$header$patient,as.character(l$header$patient.bd),l$header$patient.sex, l$header$modality, l$header$object.name, l$header$ref.pseudo,
      subobj.nb, l$header$description, nb,
      max.pix,
      l$header$object.alias, paste(l$header$file.basename,collapse=";"))
  }))
  
  colnames(base.n) <- c ("PIN", "birth.date","sex", "modality", "obj", "ref.pseudo", "nb.of.subobject" ,"description","nb", "max","object.alias", "file.basename")
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
  
  
  
  reg.list <-lapply (obj.list, function (fn){
    h <- fn$header
    if (h$modality!="reg") return(list(ref.pseudo =  h$ref.pseudo,
                                       ref = h$frame.of.reference, 
                                       patient =h$patient,
                                       patient.bd =h$patient.bd,
                                       patient.sex=h$patient.sex))
    return(list(ref.pseudo =  h$ref.pseudo,
                ref = h$frame.of.reference, reg.list=h,
                patient =h$patient,
                patient.bd =h$patient.bd,
                patient.sex=h$patient.sex))
  })
  l$T.MAT <- .load.T.MAT.by.reglist (reg.list)
  modality <- sort (unique (base.n$modality[base.n$modality!="reg"]))
  obj <- lapply(modality, function (m) {
    alias <-  sort(base.n$object.alias[(base.n$modality==m)])
    match.index <- match(alias, sapply (obj.list, function(l) l$header$object.alias))
    lobj.l <- list()
    for (idx in 1:length(alias)) lobj.l[[idx]] <- .load.object (Lobj = obj.list[[match.index[idx]]], data=data, raw.data.list=obj.list)
    names (lobj.l) <- alias
    return (lobj.l)
  })
  names(obj) <- modality
  l <- do.call(c, list(l, obj))
  
  if (dvh){
    l$dicom.dvh <- lapply(obj.list, function(obj) return(obj$header$dvh))
    l$dicom.dvh <- l$dicom.dvh [which(!sapply (l$dicom.dvh, is.null))]
    
    if (length( l$dicom.dvh)>0){
      for (dvh.idx in 1:length(l$dicom.dvh)){
        ref.object <-  names(obj.list)[(sapply(obj.list,function(dl) dl$header$object.info$SOP.ID)==l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.ID) & 
                                         (sapply(lapply(obj.list, function(dl) dl$header$object.info$SOP.label), function(v) l$dicom.dvh[[dvh.idx]]$ref.object.info$SOP.label %in% v))]
        
        if (length(ref.object)>0) l$dicom.dvh[[dvh.idx]]$ref.object.alias <-ref.object
        
        if (l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]!=""){
          ma <- match(l$dicom.dvh[[dvh.idx]]$info$number, obj.list[[l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]]]$header$roi.info$number)
          l$dicom.dvh[[dvh.idx]]$info$roi.name[!is.na(ma)] <- 
            obj.list[[l$dicom.dvh[[dvh.idx]]$ref.object.alias[1]]]$header$roi.info$name[ma[!is.na(ma)]]
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