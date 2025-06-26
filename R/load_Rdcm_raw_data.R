############################################################################
#' Loading a *.Rdcm file
#' @description the \code{load.Rdcm.raw.data} function loads the content of a *.Rdcm file.
#' @param Rdcm.filename Character string, representing the full name of a *.Rdcm 
#' file created by \link[espadon]{dicom.to.Rdcm.converter}.
#' @param address Boolean. If TRUE, a dataframe with the address of the tags in 
#' the raw DICOM data is returned.
#' @param data Boolean. If TRUE, the DICOM information are returned as an R list.
#' @param upgrade.to.latest.version Boolean. If \code{TRUE}, the function attempts 
#' to upgrade to the latest version, parsing the DICOM data. It may take longer 
#' to load the data. Consider using the \link[espadon]{Rdcm.upgrade} function.
#' @return Returns a list containing the information, converted by \pkg{espadon}, of a 
#' DICOM object..
#' @seealso \link[espadon]{dicom.to.Rdcm.converter}, \link[espadon]{load.obj.from.Rdcm}. 
#' @examples
#' # For testing, save first toy.dicom.raw () raw data to a temporary file, and
#' # convert it in Rdcm fie
#' pat.src.dir <- file.path (tempdir(), "PM_dcm") 
#' dir.create (pat.src.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "PM_rtplan", tmpdir = pat.src.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' pat.dir <- file.path (tempdir(), "PM_Rdcm")
#' dicom.to.Rdcm.converter (pat.src.dir, pat.dir, update = TRUE)
#' lf <- list.files (pat.dir, pattern = "[.]Rdcm$", full.names = TRUE)
#' lf
#' 
#' # Inspect Rdcm raw data
#' L <- load.Rdcm.raw.data (lf[1])
#' str (L, max.level =3)
#' 
#' # Cleaning  temporary directory
#' unlink (pat.src.dir, recursive = TRUE)  
#' @export
# @importFrom qs  qdeserialize 
#' @importFrom qs2 qs_deserialize
load.Rdcm.raw.data <- function (Rdcm.filename, address= TRUE, data=TRUE, 
                                upgrade.to.latest.version = FALSE) {
  if (!file.exists(Rdcm.filename)) return (NULL)
  zz <-  file(Rdcm.filename, "rb")
  on.exit(close(zz))
  l <- readBin(zz,what="int",size=4, n=3, endian="little")

  h_raw <- readBin(zz, what = "raw", n = l[1])
  h <- tryCatch(qs_deserialize(h_raw), error = function(e) NULL)
  old_qs <- is.null(h)
  if (old_qs){
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop(paste0("To read ",basename(Rdcm.filename),", install 'qs' : install.packages('qs')."))
    } else {
      warning(paste("the serialization of",basename(Rdcm.filename),
      "is no longer up-to-date, consider updating your Rdcm files with the Rdcm.upgrade function."),
      call. = FALSE)
      h <- qs::qdeserialize (h_raw)
      }
  }
  
  
  a <- NULL
  d <- NULL
  
  if (is.null(h$espadon.version)){
    espadon.version <- "0.0.0"
  } else {
    espadon.version <- h$espadon.version
    h$espadon.version <- NULL
  }
  
  from.dcm <- l[2]>0
  version <- as.numeric(strsplit(espadon.version,"[.]")[[1]])
  #correction à apporter
  correction <- list(
    version0 = espadon.version=="0.0.0",
    rtplan = espadon.version=="0.0.0" & h$modality == "rtplan",
    nopatient = is.null(h$patient.name),
    acq.date = (version[1]<1) ||  (version[1]==1 & version[2]<3) || (version[1]==1 & version[2]==3 & version[3]<0)
  )
  ################################
  update.needed <- any(unlist(correction))
  if (from.dcm & (address | data | (upgrade.to.latest.version & update.needed))){
    a_raw <- readBin(zz, what = "raw", n = l[2])
    if (old_qs){
      a <- qs::qdeserialize (a_raw)
    } else{
      a <- qs_deserialize (a_raw)
      }
    }
    
  if (data | (update.needed & from.dcm & upgrade.to.latest.version)) {
    d_raw <- readBin(zz, what = "raw", n = l[3])
    if (old_qs){
      d <- qs::qdeserialize (d_raw)
    } else{
      d <- qs_deserialize (d_raw)
    }
  }
  
  
  if (update.needed){
    if (correction$version0) {
      n <- names(h)
      idx <- which(n=="ref.object.name")
      if (length(idx)!=0) n[idx] <- "ref.object.alias"
      idx <- which(n=="patient.xyz0")  
      if (length(idx)!=0) n[idx] <- "xyz0"  
      idx <- which(n=="patient.orientation")  
      if (length(idx)!=0) n[idx] <- "orientation"  
      names(h) <- n
    }
    
    if (correction$rtplan) {
      n <- colnames(h$fraction.info )
      idx <- which(n=="planned.frac.nb")
      if (length(idx)!=0) n[idx] <- "nb.of.frac.planned"
      idx <- which(n=="beam.nb")  
      if (length(idx)!=0) n[idx] <- "nb.of.beam"  
      idx <- which(n=="brachy.app.nb")  
      if (length(idx)!=0) n[idx] <- "nb.of.brachy.app"  
      colnames(h$fraction.info ) <- n
    }
    
    if (correction$nopatient){
      n <- names(h)
      idx <-grep("^patient$",n)
      h <- c(h[1:idx], patient.name="", h[(idx+1):length(n)])
    }
  }
  
  if (from.dcm & update.needed & upgrade.to.latest.version) {
    
    colnames.df <- c("reference", 
                     "acquisition.date", "study.date","study.time","creation.date", 
                     "modality", 
                     "SOP.ID","transfer.syntax.UID", "implementation.ID",
                     "SOP.type", 
                     "scanning.sequence", "study.description", "serie.description", 
                     "study.ID", "study.UID","serie.UID", 
                     "PID",
                     "SOP.label.nb") 
    
    tag <- c("[(]0020,0052[)]$", 
             "^[(]0008,0023[)]$|^[(]3006,0008[)]$|^[(]300A,0006[)]$", "^[(]0008,0020[)]$", 
             "^[(]0008,0030[)]$", "^[(]0008,0012[)]$", 
             "^[(]0008,0060[)]$",
             "^[(]0008,0016[)]$", "^[(]0002,0010[)]$", "^[(]0002,0012[)]$", 
             "^[(]0008,0008[)]$",
             "^[(]0018,0020[)]$", "^[(]0008,1030[)]$", "^[(]0008,103E[)]$", 
             "^[(]0020,0010[)]$", "^[(]0020,000D[)]$", "^[(]0020,000E[)]$",
             "^[(]0010,0020[)]$",    
             "^[(]0020,0013[)]$")
    df <- data.frame (matrix(data=rep("",length(colnames.df)), 
                             ncol = length(colnames.df)))
    colnames(df)<-colnames.df
    
    name <- names(d[[1]])

    df[1,] <- sapply(tag, function(t) {
      value<-sapply(grep(t,name),function(i) d[[1]][[i]])
      if (length(value) == 0)
        return("")
      value <- value[!is.na(value)]
      if (length(value) == 0)
        return("")
      value <- sapply(value, function(v) {
        dum <-charToRaw(v)
        dum[dum > 125 | dum < 32] <- charToRaw(" ")
        return(rawToChar(dum))})
      value <- value[trimws(value) != ""]
      if (length(value) == 0)
        return("")
      return(sort(value)[1])
    })
    
    df$outfilename <- h$object.name
    df$ref.label <- gsub("^ref","",h$ref.pseudo)
    data.l <- lapply(1:length(d),function(i) list(address=a[[i]], data=d[[i]], 
                                                filename=h$object.info$dicom.file))
    
    L <- .obj.save.by.modality ( modality= castlow.str(df[1, ]$modality), object.info=df[1, ], data.l,
                                 only.header=TRUE, Rdcm.mode=TRUE)
    L[[1]]$header$object.alias <- h$object.alias
    L[[1]]$header$ref.pseudo <- h$ref.pseudo
    L[[1]]$header$frame.of.reference <- h$frame.of.reference
    if (!is.null(L[[1]]$header$roi.info)) L[[1]]$header$roi.info$roi.pseudo <- h$roi.info$roi.pseudo
    
    L[[1]]$header <- c(L[[1]]$header, h[is.na(match(names(h),names(L[[1]]$header)))])
    
    h<- L[[1]]$header
    espadon.version <- .espadon.version()
    
    
  }
  
  update.needed = update.needed|old_qs
  h$file.dirname <- dirname (Rdcm.filename)
  h$file.basename <- basename (Rdcm.filename)
  
  
  if (!address & !data) return(list(header=h, from.dcm=from.dcm,
                                    espadon.version=espadon.version,
                                    update.needed = update.needed))
  if (address & !data) return(list(header=h, address=a, from.dcm=from.dcm,
                                   espadon.version=espadon.version,
                                   update.needed = update.needed))
  if (!address & data) return(list(header=h, data=d, from.dcm=from.dcm,
                                   espadon.version=espadon.version,
                                   update.needed = update.needed))
  if (address & data) return(list(header=h, address=a, data=d, from.dcm=from.dcm,
                                  espadon.version=espadon.version,
                                  update.needed = update.needed))
}


