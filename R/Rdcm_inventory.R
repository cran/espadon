#' Inventory of \pkg{espadon} objects from Rdcm files
#' @description The \code{Rdcm.inventory} function creates, from Rdcm files in a 
#' patient's directory, a dataframe describing objects.
#' @param dirname Character string, representing the full name of patient 
#' directory, including Rdcm files.
#' @param upgrade.to.latest.version Boolean. If \code{TRUE}, the function attempts 
#' to upgrade to the latest version, parsing the DICOM data. It may take longer 
#' to load the data. Consider using the \link[espadon]{Rdcm.upgrade} function.
#' @return Returns a dataframe, providing information of DICOM objects.
#' @export
#'
#' @examples
#' # First, save toy patient objects to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = c (4, 4, 4))
#' save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
#' save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
#' save.T.MAT (patient$T.MAT, dirname = pat.dir)
#' 
#'
#' Rdcm.inventory (pat.dir)
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

Rdcm.inventory <- function (dirname, upgrade.to.latest.version = FALSE) {
  
  lf <- list.files (dirname, pattern = "[.]Rdcm", full.names =  TRUE) 
  
  info <-lapply (lf,function(f) {
    d <- tryCatch(load.Rdcm.raw.data (f, data=FALSE, address=FALSE,
                                      upgrade.to.latest.version = upgrade.to.latest.version)$header, 
                  error = function (e) NULL)
    if (is.null(d)) return(d)
    SOP_label <- unlist(strsplit(d$file.basename,"_"))
    SOP_label <- tryCatch(gsub("^do","",SOP_label [grepl("^do",SOP_label)]),error=function(e) "")
    if (length(SOP_label)==0) SOP_label <- gsub("[.]Rdcm","",d$file.basename)
    l <-list (ref.label = gsub("ref","",d$ref.pseudo), 
              reference = d$frame.of.reference,
              acquisition.date = d$acq.date, 
              study.date = d$study.date,
              study.time = d$study.time,
              creation.date = d$creation.date, 
              modality = d$modality, 
              SOP.label.nb=SOP_label)
    l[sapply(l, is.null)] <- ""
    # if (is.null(d$object.info)){ 
    #   l1 <- list(SOP.ID = "",transfer.syntax.UID="", implementation.ID = "", SOP.type = "", study.ID  = "", 
    #              study.UID = "", serie.UID = "", scanning.sequence = "")
    # } else {
      l1 <- list(SOP.ID = d$object.info$SOP.ID,
                 transfer.syntax.UID=d$object.info$transfer.syntax.UID,
                 implementation.ID = d$object.info$implementation.ID, 
                 SOP.type =  d$object.info$SOP.type, 
                 study.ID  = d$object.info$study.ID, 
                 study.UID = d$object.info$study.UID, 
                 serie.UID = d$object.info$serie.UID, 
                 scanning.sequence = d$object.info$scanning.sequence)
      l1[sapply(l1, is.null)] <- ""
    # }
    
    l <- do.call(c,list (l,l1 ))
    description <- unlist(strsplit (d$description, "[|]"))
    l[['study.description']] <- description[1]
    l[['serie.description']] <- description[2]
    l[['PID']] <- d$patient
    l[['name']] <- d$patient.name
    l[['birthday']] <- d$patient.bd
    l[['sex']] <- d$patient.sex
    l[['outfilename']] <- d$file.basename
    l[['input.file.nb']] <- length(d$object.info$SOP.label)
    return (l)
  })
  
  info <- info [which(!sapply (info, is.null))]
  if (length (info)==0) return(NULL)
  
  n <- names(info[[1]])
  info <- do.call (rbind.data.frame, info)
  names (info) <- n
  info[is.na (info)] <- ""
  
  info <- info [,c("reference",
                   "acquisition.date", "study.date","study.time","creation.date", 
                   "modality", 
                   "SOP.ID", "transfer.syntax.UID", "implementation.ID",
                   "SOP.type",
                   "scanning.sequence", "study.description", "serie.description", 
                   "study.ID", "study.UID","serie.UID",
                   "PID", "name", "birthday","sex",
                   "ref.label", 
                   "SOP.label.nb", "outfilename", "input.file.nb")]
  
  
  
  
  info <- info[order(info$ref.label,info$SOP.label.nb,info$outfilename), ]
  
  duplicated.index <- which(duplicated(info[, c("ref.label","SOP.label.nb")]))
  info$SOP.label.nb <- suppressWarnings(as.numeric(info$SOP.label.nb))
  info$SOP.label.nb[is.na(info$SOP.label.nb)] <- 0
  info$input.file.nb <- as.numeric(info$input.file.nb)
  final <- info
  if (length(duplicated.index)>0){
    dupli<-unique(info[duplicated.index, c("ref.label","SOP.label.nb")])
    final <- final[-duplicated.index,]
    for (idx in 1:nrow(dupli)) {
      final[final$ref.label==dupli[idx,"ref.label"] & final$SOP.label.nb==dupli[idx,"SOP.label.nb"], 
            "input.file.nb"] <- sum(info[info$ref.label==dupli[idx,"ref.label"] & info$SOP.label==dupli[idx,"SOP.label.nb"], 
                                         "input.file.nb"])
      final[final$ref.label==dupli[idx,"ref.label"] & final$SOP.label.nb==dupli[idx,"SOP.label.nb"], 
            "outfilename"] <- paste(info[info$ref.label==dupli[idx,"ref.label"] & info$SOP.label==dupli[idx,"SOP.label.nb"], 
                                         "outfilename"], collapse = ";")
    }
  } 
  
  return(final)
  
}
