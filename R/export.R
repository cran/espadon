#' Export espadon objects in DICOM format
#' @description  The \code{export} function exports struct class objects and 
#' volume class objects with CT or RTDOSE modality in DICOM format.
#' @param obj espadon object of rtstruct, ct or rtdose modality to be exported.
#' @param format Format of the export. For next use.
#' @param ref.obj.list List of espadon objects which are referenced objects of \code{obj}. See Details.
#' @param use.original.UIs,use.original.ref.UIs Booleans. If \code{TRUE}, study instance UID, serie instance 
#' UID and image type attribute are those indicated in \code{$object.info} item.
#' Otherwise,They are regenerated. See Details.
#' @param file.prefix String. Prefix added to the generated filename, in case of \code{file.name} is \code{NULL}.
#' @param file.dirname String. Name of the directory in which files are generated.
#' @param file.name String. Base name of the generated files. in CT modality, a slice number is added as a suffix. 
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, whose structure it must keep. This 
#' dataframe is used to parse DICOM files.
#' @param ... Additional settings such as NAvalue (for "volume" data), '(0020,000D)', 
#' '(0020,000E)', '(0008,0008)'
#' @return Returns nothing, but generate DICOM files if conditions are required, 
#' and indicates the name or number of files created
#' @details The object you want to export may be known in a TPS thanks to these 
#' Unique Identifiers (UIs). If you want to create a DICOM object that is different 
#' and recognised by your TPS, it is important that the DICOM files you want to 
#' create have new UIs: in this case, set the \code{use.original.UIs} argument 
#' to \code{FALSE}. Otherwise, UIs in \code{$object.info} item of your object 
#' will be used. 
#' @details Your object may have been created from another DICOM object (i.e. a 
#' reference object). You can, for example, see these links with the function 
#' \link[espadon]{display.obj.links}.
#' @details If you want to keep this links, you must indicate which objects are 
#' references in the \code{ref.obj.list} argument, in the form of a list of espadon 
#' objects. If these reference objects have their own identifier and you wish to 
#' keep them, you must set the \code{use.original.ref.UIs} argument to \code{TRUE}.
#' @details Otherwise (\code{use.original.ref.UIs=FALSE}), the UIs of the reference 
#' objects will be regenerated. It is therefore important that the reference objects 
#' contain all their data such as \code{vol3D.data} or \code{roi.data}.
#' @details It may be useful to impose a study number (tag '(0020,000D)'), serial 
#' number (tag '(0020,000E)'), or your Image Type Attribute (tag '(0008,0008)'). 
#' In this case, you need to add the arguments#' 
#' \itemize{
#' \item  \code{'(0020,000D)' = your_study_UID}, 
#' \item  \code{'(0020,000E)' = your_serial_UID}, 
#' \item  \code{'(0008,0008)' = your_image_type_attribute}.
#' }
#' For MRI, add the argument 
#'\itemize{
#'\item  \code{'(0018,0020)' = your_scanning_sequence}. It will be set to 'SE' otherwise
#'\item  \code{'(0018,0021)' = your_sequence_variant}. It will be set to 'NONE' otherwise
#'}
#' You can change patient position by adding the argument :
#'\itemize{
#'\item  \code{'(0018,5100)' = your_patient_position}. It will be set to 'HFS' otherwise
#'}
#' @examples
#' # First, save toy patient objects to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "", 
#'                              dxyz = c (6, 6, 6))
#' dicom.dir <- file.path (tempdir(), "PM_dcm") 
#' export(patient$rtstruct[[1]], ref.obj.list = list (patient$ct[[1]]),
#'        file.dirname = dicom.dir,file.name="RS")
#' export(patient$ct[[1]], file.dirname = dicom.dir,file.name="CT")
#' list.files(dicom.dir)
#'
#' # check that the links have been preserved.
#' pat <- load.patient.from.dicom (dicom.dir, verbose = FALSE)
#' display.obj.links (pat)
#' 
#' # Cleaning  temporary directories
#' unlink (dicom.dir, recursive = TRUE)
#' 
#' @export
export <- function(obj,format="dcm", ref.obj.list = NULL, use.original.UIs = FALSE,
                   use.original.ref.UIs = TRUE,
                   file.prefix ="", file.dirname= '.', file.name = NULL, 
                   tag.dictionary = dicom.tag.dictionary(), ...){
  if (format[1] == "dcm"){
    if (!(is(obj,"volume") | is(obj,"struct"))) 
      stop("At present, only espadon objects with the ct, rtdose or rtstruct modality can be exported.")
    if  (!(obj$modality %in% c("ct","mr","rtdose","rtstruct"))) 
      stop("At present, only espadon objects with the ct, mr, rtdose or rtstruct modality can be exported.")
    if (tolower(obj$modality)=="ct" | tolower(obj$modality)=="mr") {
      export.img3Dplan.to.dicom (obj, ref.obj.list = ref.obj.list, use.original.UIs = use.original.UIs,
                                 use.original.ref.UIs=use.original.ref.UIs,
                                 file.prefix =file.prefix, file.dirname= file.dirname, file.name = file.name, 
                                 tag.dictionary = tag.dictionary,...)
    } else if (tolower(obj$modality)=="rtdose") {
      export.rtdose.to.dicom (obj, ref.obj.list = ref.obj.list, use.original.UIs = use.original.UIs,
                              use.original.ref.UIs=use.original.ref.UIs,
                              file.prefix =file.prefix, file.dirname= file.dirname, file.name = file.name, 
                              tag.dictionary = tag.dictionary,...)
    } else if (tolower(obj$modality)=="rtstruct"){ 
      export.rtstruct.to.dicom (obj, ref.obj.list = ref.obj.list, use.original.UIs = use.original.UIs,
                                use.original.ref.UIs=use.original.ref.UIs,
                                file.prefix =file.prefix, file.dirname= file.dirname, file.name = file.name, 
                                tag.dictionary = tag.dictionary,...)
    }
  }
}



################################################################################
export.img3Dplan.to.dicom <- function(obj, ref.obj.list = NULL, use.original.UIs = FALSE,
                                      use.original.ref.UIs = TRUE,
                                      file.prefix ="", file.dirname= '.', file.name = NULL, 
                                      tag.dictionary = dicom.tag.dictionary(),...){
  
  rownames(tag.dictionary) <- tag.dictionary$tag
  args <- tryCatch(list(...), error = function(e)list())
  study.UID <- .ascii.to.hex (args[["(0020,000D)"]], tag.dictionary["(0020,000D)","VR"])
  serie.UID <- .ascii.to.hex (args[["(0020,000E)"]], tag.dictionary["(0020,000E)","VR"])
  sop.type <- .ascii.to.hex (args[["(0008,0008)"]], tag.dictionary["(0008,0008)","VR"])
  patient.position <- .ascii.to.hex (args[["(0018,5100)"]], tag.dictionary["(0018,5100)","VR"])
  scanning.sequence <- .ascii.to.hex (args[["(0018,0020)"]], tag.dictionary["(0018,0020)","VR"])
  sequence.variant <- .ascii.to.hex (args[["(0018,0021)"]], tag.dictionary["(0018,0021)","VR"])
  NAvalue <- args[["NAvalue"]]
 
  study.vect <- c("patient","frame.of.reference")
  serie.vect <- c("patient","frame.of.reference","modality", "description", "study.date", "study.time")
  if (!is.null(obj$object.info) &  use.original.UIs){
    if (length(study.UID)==0) study.UID <- .ascii.to.hex (obj$object.info$study.UID, tag.dictionary["(0020,000D)","VR"])
    if (length(serie.UID)==0) serie.UID <- .ascii.to.hex (obj$object.info$serie.UID, tag.dictionary["(0020,000E)","VR"])
    if (length(sop.type)==0) sop.type <- .ascii.to.hex (obj$object.info$SOP.type, tag.dictionary["(0008,0008)","VR"])
  } else {
    if (length(study.UID)==0) study.UID <- create.UID(obj,study.vect,size = 50)
    if (length(serie.UID)==0) serie.UID <- create.UID(obj,serie.vect,size =50)
    if (length(sop.type)==0) sop.type <- .ascii.to.hex ('ORIGINAL\\PRIMARY', tag.dictionary["(0008,0008)","VR"])
  }
  if (length(patient.position)==0) patient.position <- .ascii.to.hex ("HFS", tag.dictionary["(0018,5100)","VR"])
  
  if (length(.ascii.to.hex(obj$frame.of.reference, "UI"))==0) 
    obj$frame.of.reference <-  rawToChar(create.UID (obj,c("patient", "ref.pseudo"), size =44))
  frame.of.reference <- .ascii.to.hex(obj$frame.of.reference, "UI")
  frame.of.reference <- c (unlist (.value.to.raw (length (frame.of.reference), 2, TRUE)),
                           frame.of.reference) 
  
  
  
  intercept <- NA
  slope <- NA
  vol3D <- round(obj$vol3D.data,6)
  
  if (length(vol3D)!=0){
    f <- is.na(vol3D)
    if (all(f)){
      warning(c("obj$vol3D.data contains only NA values. This data will not be exported."))
      vol3D <- NULL
    } else {
      obj.max.pixel  <- max(vol3D, na.rm =TRUE)
      obj.min.pixel  <- min(vol3D, na.rm =TRUE)
      intercept <- obj.min.pixel 
      slope <- sort(diff(sort(unique(as.vector(vol3D-intercept)))))[1]
      max.pixel <- round((obj.max.pixel-intercept)/slope)
      min.pixel <- round((obj.min.pixel-intercept)/slope)
      
      while((abs(max.pixel) > 2^15) |  (abs(min.pixel) > 2^15)){
        slope <- slope*2
        max.pixel <- round((obj.max.pixel-intercept)/slope)
        min.pixel <- round((obj.min.pixel-intercept)/slope)
      }
      
      if (any(f)) {
        if (is.null(NAvalue)){
          # NAvalue <-  - 2^ceiling(log(slope* max.pixel)/log(2))
          if (obj$modality == "ct") {NAvalue <-  - 2^15} else {NAvalue <- 0}
          warning(c("obj$vol3D.data contains NA values. This data will be replaced by the value ",NAvalue, 
                    ". If you wish to change this value, add the argument NAvalue=your_NA_value."))
        }
        vol3D[f] <- NAvalue
        min.pixel <- NAvalue
      }
      

      vol3D <- round((vol3D-intercept)/slope)
    }
  } else {
    warning(c("There is no obj$vol3D.data in ", obj$object.alias,"."))
  }
  
  L00081140 <- c()
  SOP <- create.SOPUID(obj)
  
  ######################################
  fn <- sapply(SOP$l1155, function(uid)  rev(unlist(strsplit(rawToChar(uid[uid!=as.raw(0)]),"[.]")))[1])
  
  if (length(file.name)==0){
    filen <- file.path(file.dirname, paste0(file.prefix, fn, ".dcm"))
  } else {
    le <- nchar(fn[1])
    suffix <- substr(fn,le + 1 -nchar(as.character(max(obj$k.idx))),le)
    if (length(suffix)==1) suffix <- ""
    file.name <- gsub("[.]dcm$","",file.name)
    filen <- file.path(file.dirname, paste0(file.name,"_", suffix,".dcm"))
  }
  dn <- unique(dirname(filen))
  for (fdn in dn) if (!dir.exists(fdn)) dir.create(fdn,recursive = T)
  
  #######################################
  
  SOP_UID.1150 <- lapply(SOP$l1150,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  SOP_UID.1155 <- lapply(SOP$l1155,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  
  version_name <- paste(c("R ",version$major,".",version$minor,
                          " espadon ", .espadon.version()),collapse ="")
  version_name <-  c(charToRaw(version_name),as.raw(0x00))
  version_name <- version_name[1:(2*floor(length(version_name)/2))]
  
  transfert_syntax <- c(charToRaw("1.2.840.10008.1.2.1"),as.raw(0x00))
  implementation_UID <- c( charToRaw(paste0(.espadon.version(),".",
                                            paste(charToRaw(obj$modality),collapse=""))),as.raw(0x00))
  
  ####################################
  preamble  <-  as.raw(c(rep(0,128),68,73,67,77))
  
  ####################################
  tag0002 <- c("(0002,0000)","(0002,0001)" ,"(0002,0002)", "(0002,0003)", 
               "(0002,0010)", "(0002,0012)", "(0002,0013)")
  
  L0002 <- lapply(tag0002,function(t) c())
  names(L0002) <- tag0002
  
  
  L0002[["(0002,0001)"]] <- c(.tag.to.hex ("(0002,0001)"), charToRaw (tag.dictionary["(0002,0001)","VR"]),
                              as.raw(c(0x00,0x00,0x02,0x00,0x00,0x00,0x00,0x01)))
  
  L0002[["(0002,0010)"]] <- c(.tag.to.hex ("(0002,0010)"), charToRaw (tag.dictionary["(0002,0010)","VR"]),
                              unlist(.value.to.raw (length(transfert_syntax),2,TRUE)),transfert_syntax)
  L0002[["(0002,0012)"]] <- c(.tag.to.hex ("(0002,0012)"), charToRaw (tag.dictionary["(0002,0012)","VR"]), 
                              unlist(.value.to.raw (length(implementation_UID),2,TRUE)),implementation_UID)
  L0002[["(0002,0013)"]] <- c(.tag.to.hex ("(0002,0013)"), charToRaw (tag.dictionary["(0002,0013)","VR"]), 
                              unlist(.value.to.raw (length(version_name),2,TRUE)), version_name)
  
  
  
  ####################################
  tag0008 <- c("(0008,0005)", "(0008,0008)", #(optional)
               "(0008,0012)", "(0008,0016)", "(0008,0018)", "(0008,0020)", "(0008,0022)",
               "(0008,0023)", "(0008,0030)", "(0008,0033)", "(0008,0050)", 
               "(0008,0060)", "(0008,0070)", "(0008,0090)",
               "(0008,1030)", "(0008,103E)", "(0008,1070)")
  VR0008 <- tag.dictionary[tag0008,"VR"]
  L0008 <- lapply(tag0008,function(t) c())
  names(L0008) <- tag0008
  names(VR0008) <- tag0008
  
  
  L0008[["(0008,0005)"]] <- c (.tag.to.hex ("(0008,0005)"), charToRaw (tag.dictionary["(0008,0005)","VR"]),
                               as.raw (c(0x0A,0x00)), charToRaw ("ISO_IR 100"))
  
  L0008[["(0008,0008)"]] <- c (.tag.to.hex ("(0008,0008)"), charToRaw (tag.dictionary["(0008,0008)","VR"]), 
                               unlist (.value.to.raw (length (sop.type), 2, TRUE)), sop.type)
  
  if (obj$creation.date!=""){
    new.raw <- .ascii.to.hex (obj$creation.date, tag.dictionary["(0008,0012)","VR"])
    m <- match(new.raw[1:8], as.raw(48:57))
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex ("(0008,0012)"), charToRaw (tag.dictionary["(0008,0012)","VR"]),
                                 as.raw (c (0x08,0x00)), new.raw[1:8])
  } else {
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex("(0008,0012)"), charToRaw (tag.dictionary["(0008,0012)","VR"]),
                                 as.raw (c(0x08,0x00)),charToRaw (format (Sys.Date (), "%Y%m%d"))) # devient la date d'assignation de SOP_UID
  }
  

  new.raw <- .ascii.to.hex (obj$study.date, tag.dictionary["(0008,0020)","VR"])
  L0008[["(0008,0020)"]]<-  c (.tag.to.hex ("(0008,0020)"), charToRaw (tag.dictionary["(0008,0020)","VR"]),
                               as.raw (c (0x08,0x00)), new.raw[1:8])
  
  
  new.raw <- .ascii.to.hex (trimws (obj$acq.date), VR0008["(0008,0022)"])
  L0008[["(0008,0022)"]] <- c (.tag.to.hex ("(0008,0022)"), charToRaw (tag.dictionary["(0008,0022)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  L0008[["(0008,0023)"]] <- c (.tag.to.hex ("(0008,0023)"), charToRaw (tag.dictionary["(0008,0023)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  
  new.raw <- .ascii.to.hex (trimws (obj$study.time), tag.dictionary["(0008,0030)","VR"])
  L0008[["(0008,0030)"]] <- c (.tag.to.hex ("(0008,0030)"), charToRaw (tag.dictionary["(0008,0030)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  L0008[["(0008,0033)"]] <- c (.tag.to.hex ("(0008,0033)"), charToRaw (tag.dictionary["(0008,0033)","VR"]), as.raw (c (0x00,0x00)))
  
  L0008[["(0008,0050)"]] <- c (.tag.to.hex ("(0008,0050)"), charToRaw (tag.dictionary["(0008,0050)","VR"]), as.raw (c (0x00,0x00)))
  
  if (obj$modality == "ct") {
    L0008[["(0008,0060)"]] <- c (.tag.to.hex ("(0008,0060)"), charToRaw (tag.dictionary["(0008,0060)","VR"]), 
                                 as.raw (c (0x02,0x00)), charToRaw ("CT"))
  } else {
    L0008[["(0008,0060)"]] <- c (.tag.to.hex ("(0008,0060)"), charToRaw (tag.dictionary["(0008,0060)","VR"]), 
                                 as.raw (c (0x02,0x00)), charToRaw ("MR"))
    
  }
  
  
  new.raw <- .ascii.to.hex ("R ESPADON MANUFACTURER", tag.dictionary["(0008,0070)","VR"])
  L0008[["(0008,0070)"]] <- c (.tag.to.hex ("(0008,0070)"), charToRaw (tag.dictionary["(0008,0070)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  L0008[["(0008,0090)"]] <- c (.tag.to.hex ("(0008,0090)"), charToRaw (tag.dictionary["(0008,0090)","VR"]), as.raw (c (0x00,0x00)))
  
  description <- unlist(strsplit(obj$description,"[|]")[[1]])
  new.raw <-.ascii.to.hex( iconv(description[1], "ISO-IR-100"), tag.dictionary["(0008,1030)","VR"])
  L0008[["(0008,1030)"]]<- c (.tag.to.hex ("(0008,1030)"), charToRaw (tag.dictionary["(0008,1030)","VR"]),  
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  if(!is.na(description[2])!=0){
    new.raw <-.ascii.to.hex( iconv(description[2], "ISO-IR-100"), tag.dictionary["(0008,103E)","VR"])
    L0008[["(0008,103E)"]]<- c (.tag.to.hex ("(0008,103E)"), charToRaw (tag.dictionary["(0008,103E)","VR"]),
                                unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  }
  
  #
  L0008[["(0008,1070)"]]<-  c (.tag.to.hex ("(0008,1070)"), charToRaw (tag.dictionary["(0008,1070)","VR"]), as.raw (c (0x00,0x00)))
  
  
  
  ####################################
  tag0010 <- c("(0010,0010)","(0010,0020)","(0010,0030)","(0010,0040)")
  L0010 <- lapply(tag0010,function(t) c())
  names(L0010) <- tag0010
  
  # if(nchar(obj$patient.name)!=0) obj$patient.name <- ""
  new.raw <- .ascii.to.hex( iconv(obj$patient.name, "ISO-IR-100"), tag.dictionary["(0010,0010)","VR"])
  L0010[["(0010,0010)"]]<-  c (.tag.to.hex("(0010,0010)"), charToRaw(tag.dictionary["(0010,0010)","VR"]), 
                               unlist(.value.to.raw (length(new.raw),2,T)),new.raw)
  
  # if(nchar(obj$patient)!=0) obj$patient <- ""
  new.raw <-.ascii.to.hex( iconv(obj$patient, "ISO-IR-100"), tag.dictionary["(0010,0020)","VR"])
  L0010[["(0010,0020)"]]<-  c (.tag.to.hex("(0010,0020)"), charToRaw(tag.dictionary["(0010,0020)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.bd), tag.dictionary["(0010,0030)","VR"])
  L0010[["(0010,0030)"]]<-  c (.tag.to.hex ("(0010,0030)"), charToRaw (tag.dictionary["(0010,0030)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.sex), tag.dictionary["(0010,0040)","VR"])
  L0010[["(0010,0040)"]]<-  c (.tag.to.hex ("(0010,0040)"), charToRaw (tag.dictionary["(0010,0040)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  ####################################
  tag0018 <- c("(0018,0050)", "(0018,0060)", "(0018,1020)","(0018,5100)")
  L0018 <- lapply(tag0018,function(t) c())
  names(L0018) <- tag0018
  if (tolower(obj$modality)=="mr" ) {
    L0018 <- c(list("(0018,0020)" =NULL, "(0018,0021)"= NULL,"(0018,0022)"= NULL,"(0018,0023)"= NULL), L0018)
    
    if (length(scanning.sequence)==0) {
      scanning.sequence <-.ascii.to.hex ('SE', tag.dictionary["(0018,0020)","VR"])
      warning(c("scanning.sequence is set to 'SE'. Add the argument '(0018,0020)' if not"))
    }
    L0018[["(0018,0020)"]] <- c (.tag.to.hex("(0018,0020)"), charToRaw(tag.dictionary["(0018,0020)","VR"]),
                                 unlist(.value.to.raw (length(scanning.sequence),2,TRUE)),scanning.sequence)
    
    if (length(sequence.variant)==0) {
      sequence.variant <-.ascii.to.hex ('NONE', tag.dictionary["(0018,0021)","VR"])
      warning(c("scanning.sequence is set to 'NONE'. Add the argument '(0018,0021)' if not"))
    }
    L0018[["(0018,0021)"]] <- c (.tag.to.hex("(0018,0021)"), charToRaw(tag.dictionary["(0018,0021)","VR"]),
                                 unlist(.value.to.raw (length(sequence.variant),2,TRUE)),sequence.variant)
    
    
    L0018[["(0018,0022)"]] <- c(.tag.to.hex("(0018,0022)"), charToRaw(tag.dictionary["(0018,0022)","VR"]),
                                as.raw(c(0x00, 0x00)))
    if (length(obj$k.idx)>1) {
      new.raw <-.ascii.to.hex("3D" ,  tag.dictionary["(0018,0023)","VR"])
    } else {
      new.raw <-.ascii.to.hex("2D" ,  tag.dictionary["(0018,0023)","VR"])
    }
    L0018[["(0018,0023)"]] <- c (.tag.to.hex("(0018,0023)"), charToRaw(tag.dictionary["(0018,0023)","VR"]),
                                 unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
      
    }
    
 
  
  
  slice.thickness <- obj$slice.thickness
  if (slice.thickness==0) slice.thickness<- c()
  new.raw <- .ascii.to.hex (as.character(slice.thickness), tag.dictionary["(0018,0050)","VR"])
  L0018[["(0018,0050)"]] <- c (.tag.to.hex("(0018,0050)"), charToRaw(tag.dictionary["(0018,0050)","VR"]),
                               unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
  
  L0018[["(0018,0060)"]] <- c (.tag.to.hex ("(0018,0060)"), charToRaw (tag.dictionary["(0018,0060)","VR"]), as.raw (c (0x00,0x00)))
  
  new.raw <- .ascii.to.hex (paste("espadon", .espadon.version()), tag.dictionary["(0018,1020)","VR"])
  L0018[["(0018,1020)"]] <- c(.tag.to.hex("(0018,1020)"), charToRaw(tag.dictionary["(0018,1020)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
  L0018[["(0018,5100)"]] <- c(.tag.to.hex("(0018,5100)"), charToRaw(tag.dictionary["(0018,5100)","VR"]),
                              unlist(.value.to.raw (length(patient.position),2,TRUE)),patient.position)
  
  ####################################
  
  tag0020 <- c("(0020,000D)", "(0020,000E)", "(0020,0010)", "(0020,0011)", "(0020,0013)",
               "(0020,0032)", "(0020,0037)", "(0020,0052)", "(0020,1040)","(0020,1041)")
  
  L0020 <- lapply(tag0020,function(t) c())
  names(L0020) <- tag0020
  
  L0020[["(0020,000D)"]] <- c (.tag.to.hex("(0020,000D)"), charToRaw(tag.dictionary["(0020,000D)","VR"]), 
                               unlist(.value.to.raw (length(study.UID), 2, T)), study.UID)
  
  L0020[["(0020,000E)"]] <- c (.tag.to.hex("(0020,000E)"), charToRaw (tag.dictionary["(0020,000E)","VR"]), 
                               unlist(.value.to.raw (length(serie.UID), 2, TRUE)), serie.UID)
  
  
  L0020[["(0020,0010)"]] <- c(.tag.to.hex("(0020,0010)"), charToRaw(tag.dictionary["(0020,0010)","VR"]),
                              as.raw(c(0x00, 0x00)))
  L0020[["(0020,0011)"]] <- c(.tag.to.hex("(0020,0011)"), charToRaw(tag.dictionary["(0020,0011)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  # L0020[["(0020,0013)"]] <- c(.tag.to.hex("(0020,0013)"), charToRaw(tag.dictionary["(0020,0013)","VR"]),
  #                             as.raw(c(0x00, 0x00)))
  
  
  new.raw <- .ascii.to.hex (paste(obj$orientation, collapse="\\"), tag.dictionary["(0020,0037)","VR"])
  L0020[["(0020,0037)"]] <- c(.tag.to.hex("(0020,0037)"), charToRaw(tag.dictionary["(0020,0037)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  L0020[["(0020,0052)"]] <- c(.tag.to.hex("(0020,0052)"),charToRaw(tag.dictionary["(0020,0052)","VR"]), frame.of.reference)
  
  L0020[["(0020,1040)"]] <- c(.tag.to.hex("(0020,1040)"), charToRaw(tag.dictionary["(0020,1040)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  
  ####################################
  tag0028 <- c("(0028,0002)", "(0028,0004)",  "(0028,0010)", "(0028,0011)","(0028,0030)", 
               "(0028,0100)", "(0028,0101)", "(0028,0102)","(0028,0103)", "(0028,0106)", 
               "(0028,0107)", "(0028,1052)","(0028,1053)","(0028,1054)")#, "(0028,0120)")
  
  L0028 <- lapply(tag0028,function(t) c())
  names(L0028) <- tag0028
  
  new.raw <-as.raw(c(1,0))
  L0028[["(0028,0002)"]] <- c(.tag.to.hex("(0028,0002)"), charToRaw(tag.dictionary["(0028,0002)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  new.raw <- .ascii.to.hex ("MONOCHROME2", tag.dictionary["(0028,0004)","VR"])
  L0028[["(0028,0004)"]] <- c(.tag.to.hex("(0028,0004)"), charToRaw(tag.dictionary["(0028,0004)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- c(as.raw(obj$n.ijk[2]%%256) ,as.raw((obj$n.ijk[2]/256)%%256))
  L0028[["(0028,0010)"]] <- c(.tag.to.hex("(0028,0010)"), charToRaw(tag.dictionary["(0028,0010)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  new.raw <-  c(as.raw(obj$n.ijk[1]%%256) ,as.raw((obj$n.ijk[1]/256)%%256))
  L0028[["(0028,0011)"]] <- c(.tag.to.hex("(0028,0011)"), charToRaw(tag.dictionary["(0028,0011)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- .ascii.to.hex (paste(obj$dxyz[1:2], collapse="\\"), tag.dictionary["(0028,0030)","VR"])
  L0028[["(0028,0030)"]] <- c(.tag.to.hex("(0028,0030)"), charToRaw(tag.dictionary["(0028,0030)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <-as.raw(c(16,0))
  L0028[["(0028,0100)"]] <- c(.tag.to.hex("(0028,0100)"), charToRaw(tag.dictionary["(0028,0100)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  L0028[["(0028,0101)"]] <- c(.tag.to.hex("(0028,0101)"), charToRaw(tag.dictionary["(0028,0101)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  new.raw <- as.raw(c(15,0))
  L0028[["(0028,0102)"]] <- c(.tag.to.hex("(0028,0102)"), charToRaw(tag.dictionary["(0028,0102)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  # new.raw <-  as.raw(c(0,0))
  # if (obj$min.pixel<0) 
  new.raw <-  as.raw(c(1,0))
  L0028[["(0028,0103)"]] <- c(.tag.to.hex("(0028,0103)"), charToRaw(tag.dictionary["(0028,0103)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- .ascii.to.hex (as.character(intercept), tag.dictionary["(0028,1052)","VR"])
  L0028[["(0028,1052)"]] <- c(.tag.to.hex("(0028,1052)"), charToRaw(tag.dictionary["(0028,1052)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- .ascii.to.hex (as.character(slope), tag.dictionary["(0028,1053)","VR"])
  L0028[["(0028,1053)"]] <- c(.tag.to.hex("(0028,1053)"), charToRaw(tag.dictionary["(0028,1053)","VR"]),
                              unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  if (nchar(trimws(obj$unit))!=0){
    new.raw <- .ascii.to.hex (as.character(obj$unit), tag.dictionary["(0028,1054)","VR"])
    L0028[["(0028,1054)"]] <- c(.tag.to.hex("(0028,1054)"), charToRaw(tag.dictionary["(0028,1054)","VR"]),
                                unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
  }
  
  ####################################
  
  for (sop.idx in 1:length(filen)){
    L0002[["(0002,0002)"]] <- c(.tag.to.hex ("(0002,0002)"), charToRaw (tag.dictionary["(0002,0002)","VR"]),
                                SOP_UID.1150[[sop.idx]])
    L0002[["(0002,0003)"]] <- c(.tag.to.hex ("(0002,0003)"), charToRaw (tag.dictionary["(0002,0003)","VR"]),
                                SOP_UID.1155[[sop.idx]])
    
    L0002[["(0002,0000)"]] <- c(.tag.to.hex ("(0002,0000)"), charToRaw (tag.dictionary["(0002,0000)","VR"]),
                                as.raw (c (0x04,0x00)),
                                unlist (  .value.to.raw (length (unlist (L0002[2:7])), 4, T)))
    
    L0008[["(0008,0016)"]] <- c (.tag.to.hex ("(0008,0016)"), charToRaw (tag.dictionary["(0008,0016)","VR"]), 
                                 SOP_UID.1150[[sop.idx]])
    L0008[["(0008,0018)"]] <- c (.tag.to.hex ("(0008,0018)"), charToRaw (tag.dictionary["(0008,0018)","VR"]), 
                                 SOP_UID.1155[[sop.idx]])
    
    new.raw <- .ascii.to.hex (as.character(sop.idx), tag.dictionary["(0020,0013)","VR"])
    L0020[["(0020,0013)"]] <- c(.tag.to.hex("(0020,0013)"), charToRaw(tag.dictionary["(0020,0013)","VR"]),
                                unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)
    
    new.raw <- .ascii.to.hex (paste(obj$xyz0[sop.idx,], collapse="\\"), tag.dictionary["(0020,0032)","VR"])
    L0020[["(0020,0032)"]] <- c(.tag.to.hex("(0020,0032)"), charToRaw(tag.dictionary["(0020,0032)","VR"]),
                                unlist(  .value.to.raw (length(new.raw),2,TRUE)), new.raw)

    
    new.raw <- .ascii.to.hex (paste(obj$xyz0[sop.idx,3], collapse="\\"), tag.dictionary["(0020,1041)","VR"])
    L0020[["(0020,1041)"]] <- c(.tag.to.hex("(0020,1041)"), charToRaw(tag.dictionary["(0020,1041)","VR"]),
                                unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
    
    if (length(vol3D)!=0){
      map <- as.vector(vol3D[,,sop.idx])
      min.pixel <- min(map)
      max.pixel <- max(map)
      new.raw <- c(as.raw(min.pixel%%256) ,as.raw((min.pixel/256)%%256))
      L0028[["(0028,0106)"]] <- c(.tag.to.hex("(0028,0106)"), charToRaw("US"),
                                  as.raw(c(0x02,00)), new.raw)
      
      new.raw <- c(as.raw(max.pixel%%256) ,as.raw((max.pixel/256)%%256))
      L0028[["(0028,0107)"]] <- c(.tag.to.hex("(0028,0107)"), charToRaw("US"),
                                  as.raw(c(0x02,00)), new.raw)
      
      
      map <- as.vector(rbind(as.raw(map%%256) ,as.raw((map/256)%%256)))
      L7FE00010 <- c(.tag.to.hex("(7FE0,0010)"),charToRaw("OW"), as.raw(c(0,0)),
                     unlist(.value.to.raw (length(map),4,TRUE)),map)
    }
    ####################################
    rawdata <- as.vector(c(preamble,unlist(L0002),unlist(L0008),unlist(L0010),unlist(L0018),
                           unlist(L0020),unlist(L0028),L7FE00010))
    zz <- file (filen[sop.idx], "wb")
    writeBin (rawdata, zz, size = 1)
    close (zz)
  }
  
  message ("number of files created : ", length(filen))
}

################################################################################

export.rtstruct.to.dicom <- function(obj, ref.obj.list = NULL, use.original.UIs = FALSE,
                                     use.original.ref.UIs = TRUE,
                                     file.prefix ="", file.dirname= '.', file.name = NULL, 
                                     tag.dictionary = dicom.tag.dictionary(),...){
  rownames(tag.dictionary) <- tag.dictionary$tag
  args <- tryCatch(list(...), error = function(e)list())
  study.UID <- .ascii.to.hex (args[["(0020,000D)"]], tag.dictionary["(0020,000D)","VR"])
  serie.UID <- .ascii.to.hex (args[["(0020,000E)"]], tag.dictionary["(0020,000E)","VR"])
  sop.type <- .ascii.to.hex (args[["(0008,0008)"]], tag.dictionary["(0008,0008)","VR"])
  # NAvalue <- args[["NAvalue"]]
  # if (is.null(NAvalue)) NAvalue <- 0
  study.vect <- c("patient","frame.of.reference")
  serie.vect <- c("patient","frame.of.reference","modality", "description", "study.date", "study.time")
  if (!is.null(obj$object.info) &  use.original.UIs){
    if (length(study.UID)==0) study.UID <- .ascii.to.hex (obj$object.info$study.UID, tag.dictionary["(0020,000D)","VR"])
    if (length(serie.UID)==0) serie.UID <- .ascii.to.hex (obj$object.info$serie.UID, tag.dictionary["(0020,000E)","VR"])
    if (length(sop.type)==0) sop.type <- .ascii.to.hex (obj$object.info$SOP.type, tag.dictionary["(0008,0008)","VR"])
  } else {
    if (length(study.UID)==0) study.UID <- create.UID(obj,study.vect,size = 50)
    if (length(serie.UID)==0) serie.UID <- create.UID(obj,serie.vect,size =50)
    if (length(sop.type)==0) sop.type <- .ascii.to.hex ('ORIGINAL\\PRIMARY', tag.dictionary["(0008,0008)","VR"])
    
  }
  
  if (length(.ascii.to.hex(obj$frame.of.reference, "UI"))==0) 
    obj$frame.of.reference <-  rawToChar(create.UID (obj,c("patient", "ref.pseudo"), size =44))
  frame.of.reference <- .ascii.to.hex(obj$frame.of.reference, "UI")
  frame.of.reference <- c (unlist (.value.to.raw (length (frame.of.reference), 2, TRUE)),
                           frame.of.reference) 
  
  ref.SOP_UID.1150 <- NULL
  ref.SOP_UID.1155 <- NULL
  z.ref <- NULL
  if (!is.null(ref.obj.list)){
    ref.serie<- do.call(rbind, lapply(ref.obj.list,function(obj_){
      if (!is.null(obj_$object.info) & use.original.ref.UIs) return(.ascii.to.hex (obj_$object.info$serie.UID, "UI"))
      create.UID(obj_,serie.vect,size =50)})) 

    f <- !duplicated.array(ref.serie)
    ref.obj.list <- ref.obj.list[f]
    ref.serie <- ref.serie[f,, drop = FALSE]
    
    SOP.L <- lapply(ref.obj.list, function(obj_) {
      if (length(.ascii.to.hex(obj_$frame.of.reference, "UI"))==0) 
        obj_$frame.of.reference <-  rawToChar(create.UID (obj_,c("patient", "ref.pseudo"), size =44))
      if (!is.null(obj_$object.info) & use.original.ref.UIs) {
        # l1155 <- obj_$object.info$SOP.label
        # ncl.f<- cumsum(nchar(l1155))
        # ncl.d <- c(1,ncl.f[-length(ncl.f)] -1)
        # l1155 <- charToRaw(paste(l1155,collapse=""))
        # l1155 <- lapply(1:length(ncl.f), function(i) {
        #   new.raw <- c(l1155[ncl.d[i]:ncl.f[i]] ,as.raw(0x00));  new.raw[1:(2 * floor(length(new.raw)/2))]})
        l1155 <- lapply(obj_$object.info$SOP.label, function(uid) .ascii.to.hex(uid,"UI"))
        l1150 <- .ascii.to.hex(obj_$object.info$SOP.ID,"UI")
        l1150 <- lapply(l1155,function(i) l1150)
        return(list(l1150=l1150,l1155=l1155))
        
      }
      return(create.SOPUID(obj_))
    }) 
    z.ref <- lapply(ref.obj.list, function(obj_) {
      if(obj_$modality!="ct" & obj_$modality!="mr") return(NULL)
      t.mat <- ref.cutplane.add(obj_)
      (as.matrix(cbind(obj_$xyz0,1)) %*% t(get.rigid.M(t.mat, obj_$ref.pseudo,paste0(obj_$ref.pseudo, "m"))))[,3]
    })
    
    ref.SOP_UID.1150 <- lapply(SOP.L, function(l )l$l1150)
    ref.SOP_UID.1155 <- lapply(SOP.L, function(l )l$l1155)
    ref.SOP_UID.1150 <- lapply(ref.SOP_UID.1150,function(l) lapply(l, function (uid)
      c(unlist(.value.to.raw (length(uid),2,TRUE)),uid)))
    ref.SOP_UID.1155 <- lapply(ref.SOP_UID.1155,function(l) lapply(l, function (uid)
      c(unlist(.value.to.raw (length(uid),2,TRUE)),uid)))
    rm(SOP.L)
    
  }
  
  
  # use.info <- Sys.info()[1:7]
  # use.info[6] <- paste(use.info[6],Sys.getenv("LOGNAME"),sep="")
  # use.info[7] <- format(Sys.time(), "%Y%m%d%H%M%s")
  # use.info <-c(hash(charToRaw(paste(use.info, collapse = " ")), key = NULL, size = 16),
  #              unlist(.value.to.raw(as.numeric(Sys.time()),2,F)))
  # fn <- paste(use.info,collapse="")
  SOP <- create.SOPUID(obj)
  
  
  ######################################
  fn <- sapply(SOP$l1155, function(uid)  rev(unlist(strsplit(rawToChar(uid[uid!=as.raw(0)]),"[.]")))[1])
  
  if (length(file.name)==0){
    filen <- file.path(file.dirname, paste0(file.prefix, fn, ".dcm"))
  } else {
    file.name <- gsub("[.]dcm$","",file.name)
    filen <- file.path(file.dirname, paste0(file.name,".dcm"))
  }
  dn <- unique(dirname(filen))
  for (fdn in dn) if (!dir.exists(fdn)) dir.create(fdn,recursive = T)
  
  #######################################
  
  SOP_UID.1150 <- lapply(SOP$l1150,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  SOP_UID.1155 <- lapply(SOP$l1155,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  
  version_name <- paste(c("R ",version$major,".",version$minor,
                          " espadon ", .espadon.version()),collapse ="")
  version_name <-  c(charToRaw(version_name),as.raw(0x00))
  version_name <- version_name[1:(2*floor(length(version_name)/2))]
  
  transfert_syntax <- c(charToRaw("1.2.840.10008.1.2.1"),as.raw(0x00))
  implementation_UID <- c( charToRaw(paste0(.espadon.version(),".",
                                            paste(charToRaw(obj$modality),collapse=""))),as.raw(0x00))
  
  sop.idx <- 1
  # for (sop.idx  in 1:length(SOP_UID.1155)) {
  ####################################
  preamble  <-  as.raw(c(rep(0,128),68,73,67,77))
  
  tag0002 <- c("(0002,0000)","(0002,0001)" ,"(0002,0002)", "(0002,0003)", 
               "(0002,0010)", "(0002,0012)", "(0002,0013)")
  L0002 <- lapply(tag0002,function(t) c())
  names(L0002) <- tag0002
  
  
  L0002[["(0002,0001)"]] <- c(.tag.to.hex ("(0002,0001)"), charToRaw (tag.dictionary["(0002,0001)","VR"]),
                              as.raw(c(0x00,0x00,0x02,0x00,0x00,0x00,0x00,0x01)))
  
  L0002[["(0002,0002)"]] <- c(.tag.to.hex ("(0002,0002)"), charToRaw (tag.dictionary["(0002,0002)","VR"]),
                              SOP_UID.1150[[sop.idx]])
  L0002[["(0002,0003)"]] <- c(.tag.to.hex ("(0002,0003)"), charToRaw (tag.dictionary["(0002,0003)","VR"]),
                              SOP_UID.1155[[sop.idx]])
  
  L0002[["(0002,0010)"]] <- c(.tag.to.hex ("(0002,0010)"), charToRaw (tag.dictionary["(0002,0010)","VR"]),
                              unlist(.value.to.raw (length(transfert_syntax),2,TRUE)),transfert_syntax)
  L0002[["(0002,0012)"]] <- c(.tag.to.hex ("(0002,0012)"), charToRaw (tag.dictionary["(0002,0012)","VR"]), 
                              unlist(.value.to.raw (length(implementation_UID),2,TRUE)),implementation_UID)
  L0002[["(0002,0013)"]] <- c(.tag.to.hex ("(0002,0013)"), charToRaw (tag.dictionary["(0002,0013)","VR"]), 
                              unlist(.value.to.raw (length(version_name),2,TRUE)), version_name)
  
  L0002[["(0002,0000)"]] <- c(.tag.to.hex ("(0002,0000)"), charToRaw (tag.dictionary["(0002,0000)","VR"]),
                              as.raw (c (0x04,0x00)),
                              unlist (.value.to.raw (length (unlist (L0002[2:7])), 4, T)))
  
  
  
  
  ####################################
  tag0008 <- c("(0008,0005)", "(0008,0012)", "(0008,0016)", "(0008,0018)", 
               "(0008,0020)", "(0008,0021)", "(0008,0023)", "(0008,0030)",
               "(0008,0050)", "(0008,0060)", "(0008,0070)", "(0008,0090)",
               "(0008,1030)", "(0008,103E)", "(0008,1070)")
  VR0008 <- tag.dictionary[tag0008,"VR"]
  L0008 <- lapply(tag0008,function(t) c())
  names(L0008) <- tag0008
  names(VR0008) <- tag0008
  
  
  L0008[["(0008,0005)"]] <- c (.tag.to.hex ("(0008,0005)"), charToRaw (VR0008["(0008,0005)"]),
                               as.raw (c(0x0A,0x00)), charToRaw ("ISO_IR 100"))
  
  if (obj$creation.date!=""){
    new.raw <- .ascii.to.hex (obj$creation.date, VR0008["(0008,0012)"])
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex ("(0008,0012)"), charToRaw (VR0008["(0008,0012)"]),
                                 as.raw (c (0x08,0x00)), new.raw[1:8])
  } else {
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex("(0008,0012)"), charToRaw (VR0008["(0008,0012)"]),
                                 as.raw (c(0x08,0x00)),charToRaw (format (Sys.Date (), "%Y%m%d"))) # devient la date d'assignation de SOP_UID
  }
  L0008[["(0008,0016)"]] <- c (.tag.to.hex ("(0008,0016)"), charToRaw (VR0008["(0008,0016)"]), SOP_UID.1150[sop.idx])
  L0008[["(0008,0018)"]] <- c (.tag.to.hex ("(0008,0018)"), charToRaw (VR0008["(0008,0018)"]), SOP_UID.1155[sop.idx])
  
  new.raw <- .ascii.to.hex (obj$study.date, VR0008["(0008,0020)"])
  m <- match(new.raw[1:8], as.raw(48:57))
  L0008[["(0008,0020)"]]<-  c (.tag.to.hex ("(0008,0020)"), charToRaw (VR0008["(0008,0020)"]),
                               as.raw (c (0x08,0x00)), new.raw[1:8])
  
  
  new.raw <- .ascii.to.hex (trimws (obj$creation.date), VR0008["(0008,0021)"])
  L0008[["(0008,0021)"]] <- c (.tag.to.hex ("(0008,0021)"), charToRaw (VR0008["(0008,0021)"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  L0008[["(0008,0023)"]] <- c (.tag.to.hex ("(0008,0023)"), charToRaw (VR0008["(0008,0023)"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  
  new.raw <- .ascii.to.hex (trimws (obj$study.time), VR0008["(0008,0030)"])
  L0008[["(0008,0030)"]] <- c (.tag.to.hex ("(0008,0030)"), charToRaw (VR0008["(0008,0030)"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  L0008[["(0008,0050)"]]<- c (.tag.to.hex("(0008,0050)"), charToRaw (VR0008["(0008,0050)"]), as.raw (c (0x00,0x00)))
  
  L0008[["(0008,0060)"]] <- c (.tag.to.hex ("(0008,0060)"), charToRaw (VR0008["(0008,0060)"]),
                               as.raw (c (0x08,0x00)), charToRaw ("RTSTRUCT"))
  
  new.raw <- .ascii.to.hex ("R ESPADON MANUFACTURER", tag.dictionary["(0008,0070)","VR"])
  L0008[["(0008,0070)"]] <- c (.tag.to.hex ("(0008,0070)"), charToRaw (tag.dictionary["(0008,0070)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  L0008[["(0008,0090)"]] <-  c (.tag.to.hex ("(0008,0090)"), charToRaw (VR0008["(0008,0090)"]),
                                as.raw(c(0x00, 0x00)))
  
  description <- unlist(strsplit(obj$description,"[|]")[[1]])
  if (length(description) == 0) description <-""
  new.raw <-.ascii.to.hex( iconv(description[1], "ISO-IR-100"), VR0008["(0008,1030)"])
  L0008[["(0008,1030)"]]<- c (.tag.to.hex ("(0008,1030)"), charToRaw (VR0008["(0008,1030)"]),  
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  if(!is.na(description[2])!=0){
    new.raw <-.ascii.to.hex( iconv(description[2], "ISO-IR-100"), VR0008["(0008,103E)"])
    L0008[["(0008,103E)"]]<- c (.tag.to.hex ("(0008,103E)"), charToRaw (VR0008["(0008,103E)"]),
                                unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  }
  
  #
  L0008[["(0008,1070)"]]<-  c (.tag.to.hex ("(0008,1070)"), charToRaw (VR0008["(0008,1070)"]),
                               as.raw(c(0x00, 0x00)))
  
  
  ####################################
  tag0010 <- c("(0010,0010)","(0010,0020)","(0010,0030)","(0010,0040)")
  L0010 <- lapply(tag0010,function(t) c())
  names(L0010) <- tag0010
  
  # if(nchar(obj$patient.name)!=0) obj$patient.name <- ""
  new.raw <- .ascii.to.hex( iconv(obj$patient.name, "ISO-IR-100"), tag.dictionary["(0010,0010)","VR"])
  L0010[["(0010,0010)"]]<-  c (.tag.to.hex("(0010,0010)"), charToRaw(tag.dictionary["(0010,0010)","VR"]), 
                               unlist(.value.to.raw (length(new.raw),2,T)),new.raw)
  
  # if(nchar(obj$patient)!=0) obj$patient <- ""
  new.raw <-.ascii.to.hex( iconv(obj$patient, "ISO-IR-100"), tag.dictionary["(0010,0020)","VR"])
  L0010[["(0010,0020)"]]<-  c (.tag.to.hex("(0010,0020)"), charToRaw(tag.dictionary["(0010,0020)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.bd), tag.dictionary["(0010,0030)","VR"])
  L0010[["(0010,0030)"]]<-  c (.tag.to.hex ("(0010,0030)"), charToRaw (tag.dictionary["(0010,0030)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.sex), tag.dictionary["(0010,0040)","VR"])
  L0010[["(0010,0040)"]]<-  c (.tag.to.hex ("(0010,0040)"), charToRaw (tag.dictionary["(0010,0040)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  ####################################
  # "(0012,0060) "Clinical Trial Coordinating Center Name Attribute
  
  # L0012[["(0012,0060)"]] <- c(.tag.to.hex("(0012,0060)"), charToRaw(tag.dictionary["(0012,0060)","VR"]),
  #                             as.raw(c(0x00, 0x00)))
  # 
  ####################################
  tag0018 <- "(0018,1020)"
  L0018 <- lapply(tag0018,function(t) c())
  names(L0018) <- tag0018
  new.raw <- .ascii.to.hex (paste("espadon", .espadon.version()), tag.dictionary["(0018,1020)","VR"])
  L0018[["(0018,1020)"]] <- c(.tag.to.hex("(0018,1020)"), charToRaw(tag.dictionary["(0018,1020)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
  
  ####################################
  tag0020 <- c( "(0020,000D)", "(0020,000E)","(0020,0010)","(0020,0011)",
                "(0020,0052)","(0020,1040)")
  L0020 <- lapply(tag0020,function(t) c())
  names(L0020) <- tag0020
  
  
  L0020[["(0020,000D)"]] <- c (.tag.to.hex("(0020,000D)"), charToRaw(tag.dictionary["(0020,000D)","VR"]), 
                               unlist(.value.to.raw (length(study.UID), 2, T)), study.UID)
  
  
  L0020[["(0020,000E)"]] <- c (.tag.to.hex("(0020,000E)"), charToRaw (tag.dictionary["(0020,000E)","VR"]), 
                               unlist(.value.to.raw (length(serie.UID), 2, TRUE)), serie.UID)
  
  L0020[["(0020,0010)"]] <- c(.tag.to.hex("(0020,0010)"), charToRaw(tag.dictionary["(0020,0010)","VR"]),
                              as.raw(c(0x00, 0x00)))
  L0020[["(0020,0011)"]] <- c(.tag.to.hex("(0020,0011)"), charToRaw(tag.dictionary["(0020,0011)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  
  L0020[["(0020,0052)"]] <- c(.tag.to.hex("(0020,0052)"), charToRaw(tag.dictionary["(0020,0052)","VR"]), frame.of.reference)
  
  L0020[["(0020,1040)"]] <- c(.tag.to.hex("(0020,1040)"), charToRaw(tag.dictionary["(0020,1040)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  ####################################
  tag3006 <- c("(3006,0002)","(3006,0008)")
  VR3006 <- tag.dictionary[tag3006,"VR"]
  L3006 <- lapply(tag3006,function(t) c())
  names(L3006) <- tag3006
  names(VR3006) <- tag3006


  if (is.null(obj$object.alias)) {
    label <- fn
  } else {
    label <- obj$object.alias
    le <- nchar(label)
    if (le ==0){
      label <-fn
    } else if (grepl('[_]ref[[:digit:]]+',label) & grepl('[_]do[[:digit:]]+',label)){
      label <- strsplit(obj$object.alias,"[_]")[[1]]
      label <- paste0(label[grep("^ref",label)],"_", label[grep("^do",label)])
    } 
  }
  le <- nchar(label)
  new.raw <- .ascii.to.hex (substr(label,max(1,le-15),le) , VR3006["(3006,0002)"])
  
  L3006[["(3006,0002)"]] <- c(.tag.to.hex("(3006,0002)"), charToRaw(VR3006["(3006,0002)"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  L3006[["(3006,0008)"]] <- c(.tag.to.hex("(3006,0008)"), charToRaw(VR3006["(3006,0008)"]),
                              as.raw(c(0x00, 0x00)))
  
  L30060010 <- c()
  if (!is.null( ref.SOP_UID.1150)){#(3006,0010)
    
    L30060016loop <- lapply(1:length(ref.SOP_UID.1150), function(ref.idx) 
      return(unlist(lapply(1:length(ref.SOP_UID.1150[[ref.idx]]), function(idx){
        tot <- c(.tag.to.hex("(0008,1150)"),charToRaw(tag.dictionary["(0008,1150)","VR"]),
                 ref.SOP_UID.1150[[ref.idx]][[idx]],
                 .tag.to.hex("(0008,1155)"),charToRaw(tag.dictionary["(0008,1155)","VR"]),
                 ref.SOP_UID.1155[[ref.idx]][[idx]])
        c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(tot),4,TRUE)), tot)}))
      ))
    
    
    L300060014loop <- unlist(lapply(1:length(L30060016loop),function(idx){
      tot <- c(.tag.to.hex("(0020,000E)"), charToRaw(tag.dictionary["(0020,000E)","VR"]),
               unlist(.value.to.raw (length(ref.serie[idx,]),2,TRUE)),ref.serie[idx,],
               .tag.to.hex("(3006,0016)"), charToRaw(tag.dictionary["(3006,0016)","VR"]), as.raw(c(0x00,0x00)), 
               unlist(.value.to.raw (length(L30060016loop[[idx]]),4,TRUE)), L30060016loop[[idx]])
      c(as.raw(c(0xfe, 0xff, 0x00, 0xe0)),
        unlist(.value.to.raw (length(tot),4,TRUE)),tot)
    }))
    
    
    
    tag30060012 <- c("(0008,1150)","(0008,1155)","(3006,0014)")
    VR30060012 <- tag.dictionary[tag30060012,"VR"]
    L30060012 <- list()
    L30060012[[1]] <-  c(.tag.to.hex("(0008,1150)"), charToRaw(tag.dictionary["(0008,1150)","VR"]),
                         SOP_UID.1150[[sop.idx]])
    L30060012[[2]] <-  c(.tag.to.hex("(0008,1155)"), charToRaw(tag.dictionary["(0008,1155)","VR"]),
                         SOP_UID.1155[[sop.idx]])
    L30060012[[3]] <-  c(.tag.to.hex("(3006,0014)"), charToRaw(tag.dictionary["(3006,0014)","VR"]) , 
                         as.raw(c(0x00,0x00)), 
                         unlist(.value.to.raw (length(L300060014loop),4,TRUE)), L300060014loop)
    L30060012 <- unlist(L30060012)
    L30060012 <- c(as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(L30060012),4,TRUE)),L30060012)
    
    
    tag30060010 <- c("(0020,0052)","(3006,0012)")
    L30060010 <- lapply(tag30060010,function(t) c())
    names(L30060010) <- tag30060010
    
    
    L30060010[["(0020,0052)"]] <- c(.tag.to.hex("(0020,0052)"),charToRaw(tag.dictionary["(0020,0052)","VR"]), frame.of.reference)
    L30060010[["(3006,0012)"]] <- c(.tag.to.hex("(3006,0012)"),charToRaw(tag.dictionary["(3006,0012)","VR"]), as.raw(c(0x00,0x00)), 
                                    unlist(.value.to.raw (length(L30060012),4,TRUE)), L30060012)
    L30060010 <- as.raw(unlist(L30060010))
    L30060010 <- c(as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(L30060010),4,TRUE)),L30060010)
    
    L30060010 <- c(.tag.to.hex("(3006,0010)"),charToRaw(tag.dictionary["(3006,0010)","VR"]), as.raw(c(0x00,0x00)), 
                   unlist(.value.to.raw (length(L30060010),4,TRUE)), L30060010)
  } 
  #####
  tag30060020loop <- c("(3006,0022)","(3006,0024)","(3006,0026)", "(3006,0028)","(3006,0036)")
  VR30060020loop <- tag.dictionary[tag30060020loop,"VR"]
  
  add.desc <- any(obj$roi.info$description!="")
  L30060020loop <-unlist(lapply(1:nrow(obj$roi.info), function(i) {
    v <-  .ascii.to.hex(as.character(as.integer(obj$roi.info$number[i])),VR30060020loop[1])
    tot <- c(.tag.to.hex("(3006,0022)"),charToRaw(tag.dictionary["(3006,0022)","VR"]), 
             unlist(.value.to.raw (length(v),2,TRUE)),v)
    
    tot <- c(tot, .tag.to.hex("(3006,0024)"),charToRaw(tag.dictionary["(3006,0024)","VR"]),  
             frame.of.reference)
    
    v <- .ascii.to.hex(iconv(obj$roi.info$name[i], "ISO-IR-100"),tag.dictionary["(3006,0026)","VR"])
    tot <- c(tot, c (.tag.to.hex("(3006,0026)"),charToRaw(tag.dictionary["(3006,0026)","VR"]),
                     unlist(.value.to.raw (length(v),2,TRUE)), v))
    if (add.desc){
      v <- .ascii.to.hex(iconv(obj$roi.info$description[i], "ISO-IR-100"),tag.dictionary["(3006,0028)","VR"])
      tot <- c(tot, c (.tag.to.hex("(3006,0028)"),charToRaw(tag.dictionary["(3006,0028)","VR"]),
                       unlist(.value.to.raw (length(v),2,TRUE)), v))
    }
    v <- .ascii.to.hex(toupper(obj$roi.info$generation.algorithm[i]),VR30060020loop[5])
    tot <- c(tot, c (.tag.to.hex(tag30060020loop[5]),charToRaw(VR30060020loop[5]),
                     unlist(.value.to.raw (length(v),2,TRUE)), v))
    c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(tot),4,TRUE)), tot)
  }))
  
  tag30060020 <- c( "(3006,0020)")
  VR30060020 <- tag.dictionary[tag30060020,"VR"]
  L30060020 <- list()
  L30060020[[1]] <- c(.tag.to.hex(tag30060020[1]),charToRaw(VR30060020[1]) , as.raw(c(0x00,0x00)), 
                      unlist(.value.to.raw (length(L30060020loop),4,TRUE)), L30060020loop)
  
  ####################################
  tag30060040 <- c( "(3006,002A)", "(3006,0040)", "(3006,0084)")
  VR30060040 <- tag.dictionary[tag30060040,"VR"]
  tag30060040loop <- c("(3006,0016)", "(3006,0042)","(3006,0046)","(3006,0048)", "(3006,0050)")
  VR30060040loop <- tag.dictionary[tag30060040loop,"VR"]
  
  L30060039loop <- unlist(lapply(1:nrow(obj$roi.info), function(i) {
    tot <- raw(0)
    if (!is.null(obj$roi.data[[i]])){
      L30060040loop <- unlist(lapply( 1:length(obj$roi.data[[i]]), function(j){
        pt <- round(as.matrix(obj$roi.data[[i]][[j]]$pt),10)# pour être sur que ça rentre sur 16octets
        le.pt <- nrow(pt)
        if (all(pt[1,]==pt[le.pt,])) {pt <- pt[-le.pt, ]; le.pt <- le.pt-1}
        tot_ <- NULL
        if (!is.null(ref.obj.list) & length(z.ref)>0){
          refz.L <- lapply (1:length(z.ref),function(ref.idx) {
            if (is.null(z.ref[[ref.idx]])) return(NULL)
            flag <- round(z.ref[[ref.idx]],6) == round(pt[1,3],6)
            if (all(!flag)){
              warning(paste(obj$roi.info$name[i],": Check Reference object at z=",
                            round(pt[1,3],6),"mm."),call. = FALSE)
              return(NULL)
            }
            tot <- c(.tag.to.hex("(0008,1150)"),charToRaw(tag.dictionary["(0008,1150)","VR"]),
                     ref.SOP_UID.1150[[ref.idx]][flag][[1]],
                     .tag.to.hex("(0008,1155)"),charToRaw(tag.dictionary["(0008,1155)","VR"]),
                     ref.SOP_UID.1155[[ref.idx]][flag][[1]])
            return(c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(tot),4,TRUE)), tot))
          })
          refz.L <- refz.L[ !sapply(refz.L,is.null)]
          if (length(refz.L)) 
            tot_ <- c(.tag.to.hex(tag30060040loop[1]),charToRaw(VR30060040loop[1]) , as.raw(c(0x00,0x00)), 
                                        unlist(.value.to.raw (length(refz.L[[1]]),4,TRUE)), refz.L[[1]])
            
        }
        
        

        v <- .ascii.to.hex(toupper(obj$roi.data[[i]][[j]]$type), VR30060040loop[2])
        tot_ <- c(tot_, .tag.to.hex(tag30060040loop[2]), charToRaw(VR30060040loop[2]),
                  unlist(.value.to.raw (length(v),2,TRUE)),v)

        v <- .ascii.to.hex(as.character(le.pt),VR30060040loop[3])
        tot_ <- c(tot_, .tag.to.hex(tag30060040loop[3]), charToRaw(VR30060040loop[3]),
                  unlist(.value.to.raw (length(v),2,TRUE)),v)
        
        v <- .ascii.to.hex(as.character(j-1),VR30060040loop[4])
        tot_ <- c(tot_, .tag.to.hex(tag30060040loop[4]), charToRaw(VR30060040loop[4]),
                  unlist(.value.to.raw (length(v),2,TRUE)),v)
        
        v <- .ascii.to.hex(paste(as.numeric(t(pt)),collapse="\\"),VR30060040loop[5])
        tot_ <- c(tot_, .tag.to.hex(tag30060040loop[5]), charToRaw(VR30060040loop[5]),
                  unlist(.value.to.raw (length(v),2,TRUE)),v)
        c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), 
           unlist(.value.to.raw (length(tot_),4,TRUE)), tot_)
      }))
      tot <- c(.tag.to.hex(tag30060040[2]),charToRaw(VR30060040[2]) , as.raw(c(0x00,0x00)), 
               unlist(.value.to.raw (length(L30060040loop),4,TRUE)), L30060040loop)
    }
    
    v <-.ascii.to.hex(paste(as.numeric(col2rgb(obj$roi.info$color[i])),collapse="\\"),VR30060040[1])
    tot <- c(.tag.to.hex(tag30060040[1]), charToRaw(VR30060040[1]),
             unlist(.value.to.raw (length(v),2,TRUE)),v, tot)
    
    v <- .ascii.to.hex(as.character(as.integer(obj$roi.info$number[i])),VR30060040[3])
    tot <- c(tot, .tag.to.hex(tag30060040[3]), charToRaw(VR30060040[3]),
             unlist(.value.to.raw (length(v),2,TRUE)),v)
    
    
    

    c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), 
       unlist(.value.to.raw (length(tot),4,TRUE)), tot)
  }))
  
  tag30060039 <- c("(3006,0039)")
  VR30060039 <- tag.dictionary[tag30060039,"VR"]
  L30060039 <- list()
  L30060039[[1]] <- c(.tag.to.hex(tag30060039[1]),charToRaw(VR30060039[1]) , as.raw(c(0x00,0x00)), 
                      unlist(.value.to.raw (length(L30060039loop),4,T)), L30060039loop)
  
  ####################################
  tag30060080loop <- c("(3006,0082)", "(3006,0084)", "(3006,0085)",
                       "(3006,00A4)", "(3006,00A6)")
  
  
  VR30060080loop <- tag.dictionary[tag30060080loop,"VR"]
  obj$roi.obs[is.na(obj$roi.obs)] <- ""
  
  L30060080loop <-unlist(lapply(1:nrow(obj$roi.obs), function(i) {
    v <-  .ascii.to.hex (as.character(as.integer(obj$roi.obs$nb[i])), VR30060080loop[1])
    tot <- c (.tag.to.hex (tag30060080loop[1]), charToRaw (VR30060080loop[1]), 
              unlist (.value.to.raw (length (v), 2, TRUE)), v)
    
    v <-  .ascii.to.hex (as.character(as.integer(obj$roi.obs$roi.nb[i])), VR30060080loop[2])
    tot <- c(tot, .tag.to.hex (tag30060080loop[2]), charToRaw (VR30060080loop[2]), 
             unlist (.value.to.raw (length (v), 2, TRUE)), v)
    
    v <-  .ascii.to.hex (iconv(obj$roi.obs$label[i], "ISO-IR-100"), VR30060080loop[3])
    tot <- c(tot, .tag.to.hex (tag30060080loop[3]), charToRaw (VR30060080loop[3]), 
             unlist (.value.to.raw (length (v), 2, TRUE)), v)
    
    v <- .ascii.to.hex(iconv(toupper(obj$roi.obs$type[i]), "ISO-IR-100"),VR30060080loop[4])
    tot <- c(tot, c (.tag.to.hex(tag30060080loop[4]),charToRaw(VR30060080loop[4]),
                     unlist(.value.to.raw (length(v),2,TRUE)), v))
    
    v <- .ascii.to.hex(iconv(obj$roi.obs$interpreter[i], "ISO-IR-100"),VR30060080loop[5])
    tot <- c(tot, c (.tag.to.hex(tag30060080loop[5]),charToRaw(VR30060080loop[5]),
                     unlist(.value.to.raw (length(v),2,TRUE)), v))
    
    c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(tot),4,T)), tot)
  }))
  
  
  tag30060080 <- c("(3006,0080)")
  VR30060080 <- tag.dictionary[tag30060080,"VR"]
  
  
  L30060080 <- list()
  L30060080[[1]] <- c(.tag.to.hex(tag30060080[1]),charToRaw(VR30060080[1]) , as.raw(c(0x00,0x00)), 
                      unlist(.value.to.raw (length(L30060080loop),4,TRUE)), L30060080loop)
  
  #############################################################################
  tag300E <- c("(300E,0002)","(300E,0004)","(300E,0005)","(300E,0008)")
  
  L300E <- lapply(tag300E,function(t) c())
  names(L300E) <- tag300E
  
  # Approval Status Attribute
  Approval.Status.Attribute <- "UNAPPROVED"
  Approval.Status.Attribute.vect <- c("APPROVED", "UNAPPROVED","REJECTED")
  if (!is.null(obj$approval.status)) {
    m <- match(toupper(trimws(obj$approval.status)), Approval.Status.Attribute.vect)
    if (!is.na(m)) Approval.Status.Attribute <- Approval.Status.Attribute.vect[m]
  }
  new.raw <- .ascii.to.hex(Approval.Status.Attribute, tag.dictionary["(300E,0002)","VR"])
  L300E[["(300E,0002)"]]<-  c (.tag.to.hex("(300E,0002)"), charToRaw( tag.dictionary["(300E,0002)","VR"]), 
                               unlist(.value.to.raw (length(new.raw),2,T)),new.raw)
  
  if (Approval.Status.Attribute %in%  c("APPROVED","REJECTED")){
    L300E[["(300E,0004)"]]<-  c (.tag.to.hex("(300E,0004)"), charToRaw( tag.dictionary["(300E,0004)","VR"]), as.raw(c(0,0)))
    L300E[["(300E,0005)"]]<-  c (.tag.to.hex("(300E,0005)"), charToRaw( tag.dictionary["(300E,0005)","VR"]), as.raw(c(0,0)))
    L300E[["(300E,0008)"]]<-  c (.tag.to.hex("(300E,0008)"), charToRaw( tag.dictionary["(300E,0008)","VR"]), as.raw(c(0,0)))
  }
  
  ####################################
  rawdata <- c(preamble,unlist(L0002),unlist(L0008),unlist(L0010),unlist(L0018),
               unlist(L0020),unlist(L3006), L30060010, unlist(L30060020), 
               unlist(L30060039), unlist(L30060080))
  zz <- file (filen, "wb")
  writeBin (rawdata, zz, size = 1)
  close (zz)
  # }
  
  message ("file created : ", basename(filen))
}

################################################################################
export.rtdose.to.dicom <- function(obj, ref.obj.list = NULL, use.original.UIs = FALSE,
                                   use.original.ref.UIs = TRUE,
                                   file.prefix ="", file.dirname= '.', file.name = NULL, 
                                   tag.dictionary = dicom.tag.dictionary(),...){
  
  rownames(tag.dictionary) <- tag.dictionary$tag
  args <- tryCatch(list(...), error = function(e)list())
  study.UID <- .ascii.to.hex (args[["(0020,000D)"]], tag.dictionary["(0020,000D)","VR"])
  serie.UID <- .ascii.to.hex (args[["(0020,000E)"]], tag.dictionary["(0020,000E)","VR"])
  sop.type <- .ascii.to.hex (args[["(0008,0008)"]], tag.dictionary["(0008,0008)","VR"])
  NAvalue <- args[["NAvalue"]]
  
  study.vect <- c("patient","frame.of.reference")
  serie.vect <- c("patient","frame.of.reference","modality", "description", "study.date", "study.time")
  
  if (!is.null(obj$object.info) &  use.original.UIs){
    if (length(study.UID)==0) study.UID <- .ascii.to.hex (obj$object.info$study.UID, tag.dictionary["(0020,000D)","VR"])
    if (length(serie.UID)==0) serie.UID <- .ascii.to.hex (obj$object.info$serie.UID, tag.dictionary["(0020,000E)","VR"])
    if (length(sop.type)==0) sop.type <- .ascii.to.hex (obj$object.info$SOP.type, tag.dictionary["(0008,0008)","VR"])
  } else {
    if (length(study.UID)==0) study.UID <- create.UID(obj,study.vect,size = 50)
    if (length(serie.UID)==0) serie.UID <- create.UID(obj,serie.vect,size =50)
    if (length(sop.type)==0) sop.type <- .ascii.to.hex ('ORIGINAL\\PRIMARY', tag.dictionary["(0008,0008)","VR"])
    
  }
  
  if (length(.ascii.to.hex(obj$frame.of.reference, "UI"))==0) 
    obj$frame.of.reference <-  rawToChar(create.UID (obj,c("patient", "ref.pseudo"), size =44))
  frame.of.reference <- .ascii.to.hex(obj$frame.of.reference, "UI")
  frame.of.reference <- c (unlist (.value.to.raw (length (frame.of.reference), 2, TRUE)),
                           frame.of.reference) 
  
  
  
  ref.modality<- sapply(ref.obj.list, function(l) l$modality)
  ref.obj.list <- ref.obj.list[ref.modality=="rtstruct"]
  
  L300C0060 <- c()# voir la référence pour chaque plan
  if (!is.null(ref.obj.list)){
    ref.serie<- do.call(rbind, lapply(ref.obj.list,function(obj_){
      if (!is.null(obj_$object.info) & use.original.ref.UIs) return(.ascii.to.hex (obj_$object.info$serie.UID, "UI"))
      create.UID(obj_,serie.vect,size =50)}))
    f <- !duplicated.array(ref.serie)
    ref.obj.list <- ref.obj.list[f]
    ref.serie <- ref.serie[f,, drop = FALSE]
    
    SOP.L <- lapply(ref.obj.list, function(obj_) {
      if (length(.ascii.to.hex(obj_$frame.of.reference, "UI"))==0)
        obj_$frame.of.reference <-  rawToChar(create.UID (obj_,c("patient", "ref.pseudo"), size =44))
      if (!is.null(obj_$object.info) & use.original.ref.UIs) {
        # l1155 <- obj_$object.info$SOP.label
        # ncl.f<- cumsum(nchar(l1155))
        # ncl.d <- c(1,ncl.f[-length(ncl.f)] -1)
        # l1155 <- charToRaw(paste(l1155,collapse=""))
        # l1155 <- lapply(1:length(ncl.f), function(i) {
        #   new.raw <- c(l1155[ncl.d[i]:ncl.f[i]] ,as.raw(0x00));  new.raw[1:(2 * floor(length(new.raw)/2))]})
        l1155 <- lapply(obj_$object.info$SOP.label, function(uid) .ascii.to.hex(uid,"UI"))
        l1150 <- .ascii.to.hex(obj_$object.info$SOP.ID,"UI")
        l1150 <- lapply(l1155,function(i) l1150)
        return(list(l1150=l1150,l1155=l1155))
        
      }
      return(create.SOPUID(obj_))
    })
    
    z.ref <- lapply(ref.obj.list, function(obj_) {
      if(obj_$modality!="ct" & obj_$modality!="mr") return(NULL)
      t.mat <- ref.cutplane.add(obj_)
      (as.matrix(cbind(obj_$xyz0,1)) %*% t(get.rigid.M(t.mat, obj_$ref.pseudo,paste0(obj_$ref.pseudo, "m"))))[,3]
    })
    
    ref.SOP_UID.1150 <- lapply(SOP.L, function(l )l$l1150)
    ref.SOP_UID.1155 <- lapply(SOP.L, function(l )l$l1155)
    ref.SOP_UID.1150 <- lapply(ref.SOP_UID.1150,function(l) lapply(l, function (uid)
      c(unlist(.value.to.raw (length(uid),2,TRUE)),uid)))
    ref.SOP_UID.1155 <- lapply(ref.SOP_UID.1155,function(l) lapply(l, function (uid)
      c(unlist(.value.to.raw (length(uid),2,TRUE)),uid)))
    rm(SOP.L)
    
    L3006060loop <- unlist(lapply(1:length(ref.SOP_UID.1150), function(ref.idx)
      return(unlist(lapply(1:length(ref.SOP_UID.1150[[ref.idx]]), function(idx){
        tot <- c(.tag.to.hex("(0008,1150)"),charToRaw(tag.dictionary["(0008,1150)","VR"]),
                 ref.SOP_UID.1150[[ref.idx]][[idx]],
                 .tag.to.hex("(0008,1155)"),charToRaw(tag.dictionary["(0008,1155)","VR"]),
                 ref.SOP_UID.1155[[ref.idx]][[idx]])
        c (as.raw(c(0xfe, 0xff, 0x00, 0xe0)), unlist(.value.to.raw (length(tot),4,TRUE)), tot)}))
      )))
    L300C0060 <- c(.tag.to.hex("(300C,0060)"), charToRaw(tag.dictionary["(300C,0060)","VR"]), as.raw(c(0x00,0x00)),
                   unlist(.value.to.raw (length(L3006060loop),4,TRUE)), L3006060loop)
    
  }
  
  
  
  SOP <- create.SOPUID(obj)
  
  ######################################
  fn <- sapply(SOP$l1155, function(uid)  rev(unlist(strsplit(rawToChar(uid[uid!=as.raw(0)]),"[.]")))[1])
  
  if (length(file.name)==0){
    filen <- file.path(file.dirname, paste0(file.prefix, fn, ".dcm"))
  } else {
    
    file.name <- gsub("[.]dcm$","",file.name)
    filen <- file.path(file.dirname, paste0(file.name,".dcm"))
  }
  dn <- unique(dirname(filen))
  for (fdn in dn) if (!dir.exists(fdn)) dir.create(fdn,recursive = T)
  
  #######################################
  
  SOP_UID.1150 <- lapply(SOP$l1150,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  SOP_UID.1155 <- lapply(SOP$l1155,function(uid) c(unlist(.value.to.raw (length(uid),2,TRUE)),uid))
  
  
  version_name <- paste(c("R ",version$major,".",version$minor,
                          " espadon ", .espadon.version()),collapse ="")
  version_name <-  c(charToRaw(version_name),as.raw(0x00))
  version_name <- version_name[1:(2*floor(length(version_name)/2))]
  
  transfert_syntax <- c(charToRaw("1.2.840.10008.1.2.1"),as.raw(0x00))
  implementation_UID <- c( charToRaw(paste0(.espadon.version(),".",
                                            paste(charToRaw(obj$modality),collapse=""))),as.raw(0x00))
  ####################################
  sop.idx <- 1
  
  ####################################
  preamble  <-  as.raw(c(rep(0,128),68,73,67,77))
  
  ####################################
  tag0002 <- c("(0002,0000)","(0002,0001)" ,"(0002,0002)", "(0002,0003)", 
               "(0002,0010)", "(0002,0012)", "(0002,0013)")
  L0002 <- lapply(tag0002,function(t) c())
  names(L0002) <- tag0002
  
  
  L0002[["(0002,0001)"]] <- c(.tag.to.hex ("(0002,0001)"), charToRaw (tag.dictionary["(0002,0001)","VR"]),
                              as.raw(c(0x00,0x00,0x02,0x00,0x00,0x00,0x00,0x01)))
  L0002[["(0002,0002)"]] <- c(.tag.to.hex ("(0002,0002)"), charToRaw (tag.dictionary["(0002,0002)","VR"]),
                              SOP_UID.1150[[sop.idx]])
  L0002[["(0002,0003)"]] <- c(.tag.to.hex ("(0002,0003)"), charToRaw (tag.dictionary["(0002,0003)","VR"]),
                              SOP_UID.1155[[sop.idx]])
  L0002[["(0002,0010)"]] <- c(.tag.to.hex ("(0002,0010)"), charToRaw (tag.dictionary["(0002,0010)","VR"]),
                              unlist(.value.to.raw (length(transfert_syntax),2,TRUE)),transfert_syntax)
  
  L0002[["(0002,0012)"]] <- c(.tag.to.hex ("(0002,0012)"), charToRaw (tag.dictionary["(0002,0012)","VR"]), 
                              unlist(.value.to.raw (length(implementation_UID),2,TRUE)),implementation_UID)
  
  L0002[["(0002,0013)"]] <- c(.tag.to.hex ("(0002,0013)"), charToRaw (tag.dictionary["(0002,0013)","VR"]), 
                              unlist(.value.to.raw (length(version_name),2,TRUE)), version_name) 
  
  L0002[["(0002,0000)"]] <- c(.tag.to.hex ("(0002,0000)"), charToRaw (tag.dictionary["(0002,0000)","VR"]),
                              as.raw (c (0x04,0x00)),
                              unlist (.value.to.raw (length (unlist (L0002[2:7])), 4, T)))
  # rawdata <- as.vector(c(preamble, unlist(L0002)))
  # dicom.browser(rawdata)
  ####################################
  tag0008 <- c("(0008,0005)", "(0008,0012)", "(0008,0016)", "(0008,0018)", "(0008,0020)", 
               "(0008,0023)", "(0008,0030)", "(0008,0050)", "(0008,0060)", "(0008,0070)", "(0008,0090)",
               "(0008,1030)", "(0008,103E)", "(0008,1070)")
  VR0008 <- tag.dictionary[tag0008,"VR"]
  L0008 <- lapply(tag0008,function(t) c())
  names(L0008) <- tag0008
  names(VR0008) <- tag0008
  
  
  L0008[["(0008,0005)"]] <- c (.tag.to.hex ("(0008,0005)"), charToRaw (tag.dictionary["(0008,0005)","VR"]),
                               as.raw (c(0x0A,0x00)), charToRaw ("ISO_IR 100"))
  
  if (obj$creation.date!=""){
    new.raw <- .ascii.to.hex (obj$creation.date, tag.dictionary["(0008,0012)","VR"])
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex ("(0008,0012)"), charToRaw (tag.dictionary["(0008,0012)","VR"]),
                                 as.raw (c (0x08,0x00)), new.raw[1:8])
  } else {
    L0008[["(0008,0012)"]]<-  c (.tag.to.hex("(0008,0012)"), charToRaw (tag.dictionary["(0008,0012)","VR"]),
                                 as.raw (c(0x08,0x00)),charToRaw (format (Sys.Date (), "%Y%m%d"))) # devient la date d'assignation de SOP_UID
  }
  L0008[["(0008,0016)"]] <- c (.tag.to.hex ("(0008,0016)"), charToRaw (tag.dictionary["(0008,0016)","VR"]), SOP_UID.1150[[sop.idx]])
  L0008[["(0008,0018)"]] <- c (.tag.to.hex ("(0008,0018)"), charToRaw (tag.dictionary["(0008,0018)","VR"]), SOP_UID.1155[[sop.idx]])
  
  new.raw <- .ascii.to.hex (obj$study.date, tag.dictionary["(0008,0020)","VR"])
  m <- match(new.raw[1:8], as.raw(48:57))
  L0008[["(0008,0020)"]]<-  c (.tag.to.hex ("(0008,0020)"), charToRaw (tag.dictionary["(0008,0020)","VR"]),
                               as.raw (c (0x08,0x00)), new.raw[1:8])
  
  
  new.raw <- .ascii.to.hex (trimws (obj$acq.date), VR0008["(0008,0023)"])
  L0008[["(0008,0023)"]] <- c (.tag.to.hex ("(0008,0023)"), charToRaw (tag.dictionary["(0008,0023)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  
  new.raw <- .ascii.to.hex (trimws (obj$study.time), tag.dictionary["(0008,0030)","VR"])
  L0008[["(0008,0030)"]] <- c (.tag.to.hex ("(0008,0030)"), charToRaw (tag.dictionary["(0008,0030)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  L0008[["(0008,0050)"]]<- c (.tag.to.hex("(0008,0050)"), charToRaw (tag.dictionary["(0008,0050)","VR"]), as.raw (c (0x00,0x00)))
  
  L0008[["(0008,0060)"]] <- c (.tag.to.hex ("(0008,0060)"), charToRaw (tag.dictionary["(0008,0060)","VR"]),
                               as.raw (c (0x06,0x00)), charToRaw ("RTDOSE"))
  
  new.raw <- .ascii.to.hex ("R ESPADON MANUFACTURER", tag.dictionary["(0008,0070)","VR"])
  L0008[["(0008,0070)"]] <- c (.tag.to.hex ("(0008,0070)"), charToRaw (tag.dictionary["(0008,0070)","VR"]), 
                               unlist (.value.to.raw (length (new.raw), 2, TRUE)), new.raw)
  
  L0008[["(0008,0090)"]] <- c (.tag.to.hex ("(0008,0090)"),charToRaw (tag.dictionary["(0008,0090)","VR"]),
                               as.raw(c(0x00, 0x00)))
  
  description <- unlist(strsplit(obj$description,"[|]")[[1]])
  new.raw <-.ascii.to.hex( iconv(description[1], "ISO-IR-100"), tag.dictionary["(0008,1030)","VR"])
  L0008[["(0008,1030)"]]<- c (.tag.to.hex ("(0008,1030)"), charToRaw (tag.dictionary["(0008,1030)","VR"]),  
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  if(!is.na(description[2])!=0){
    new.raw <-.ascii.to.hex( iconv(description[2], "ISO-IR-100"), tag.dictionary["(0008,103E)","VR"])
    L0008[["(0008,103E)"]]<- c (.tag.to.hex ("(0008,103E)"), charToRaw (tag.dictionary["(0008,103E)","VR"]),
                                unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  }
  
  #
  L0008[["(0008,1070)"]]<-  c (.tag.to.hex ("(0008,1070)"), charToRaw (tag.dictionary["(0008,1070)","VR"]),
                               as.raw(c(0x00, 0x00)))
  
  
  ####################################
  tag0010 <- c("(0010,0010)","(0010,0020)","(0010,0030)","(0010,0040)")
  L0010 <- lapply(tag0010,function(t) c())
  names(L0010) <- tag0010
  
  # if(nchar(obj$patient.name)!=0) obj$patient.name <- ""
  new.raw <- .ascii.to.hex( iconv(obj$patient.name, "ISO-IR-100"), tag.dictionary["(0010,0010)","VR"])
  L0010[["(0010,0010)"]]<-  c (.tag.to.hex("(0010,0010)"), charToRaw(tag.dictionary["(0010,0010)","VR"]), 
                               unlist(.value.to.raw (length(new.raw),2,T)),new.raw)
  
  # if(nchar(obj$patient)!=0) obj$patient <- ""
  new.raw <-.ascii.to.hex( iconv(obj$patient, "ISO-IR-100"), tag.dictionary["(0010,0020)","VR"])
  L0010[["(0010,0020)"]]<-  c (.tag.to.hex("(0010,0020)"), charToRaw(tag.dictionary["(0010,0020)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.bd), tag.dictionary["(0010,0030)","VR"])
  L0010[["(0010,0030)"]]<-  c (.tag.to.hex ("(0010,0030)"), charToRaw (tag.dictionary["(0010,0030)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  new.raw <- .ascii.to.hex (trimws(obj$patient.sex), tag.dictionary["(0010,0040)","VR"])
  L0010[["(0010,0040)"]]<-  c (.tag.to.hex ("(0010,0040)"), charToRaw (tag.dictionary["(0010,0040)","VR"]),  
                               unlist (.value.to.raw (length (new.raw), 2, T)),new.raw)
  
  
  ####################################
  
  tag0018 <- c("(0018,0050)", "(0018,1020)")
  L0018 <- lapply(tag0018,function(t) c())
  names(L0018) <- tag0018
  slice.thickness <- obj$slice.thickness
  if (slice.thickness==0) slice.thickness<- c()
  new.raw <- .ascii.to.hex (as.character(slice.thickness), tag.dictionary["(0018,0050)","VR"])
  new.raw <- .ascii.to.hex (as.character(obj$slice.thickness), tag.dictionary["(0018,0050)","VR"])
  L0018[["(0018,0050)"]] <- c(.tag.to.hex("(0018,0050)"), charToRaw(tag.dictionary["(0018,0050)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
  new.raw <- .ascii.to.hex (paste("espadon", .espadon.version()), tag.dictionary["(0018,1020)","VR"])
  L0018[["(0018,1020)"]] <- c(.tag.to.hex("(0018,1020)"), charToRaw(tag.dictionary["(0018,1020)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)),new.raw)
  
  ####################################
  
  tag0020 <- c("(0020,000D)", "(0020,000E)", "(0020,0010)", "(0020,0011)", "(0020,0013)",
               "(0020,0032)", "(0020,0037)", "(0020,0052)", "(0020,1040)")
  
  L0020 <- lapply(tag0020,function(t) c())
  names(L0020) <- tag0020
  
  L0020[["(0020,000D)"]] <- c (.tag.to.hex("(0020,000D)"), charToRaw(tag.dictionary["(0020,000D)","VR"]), 
                               unlist(.value.to.raw (length(study.UID), 2, T)), study.UID)
  
  
  L0020[["(0020,000E)"]] <- c (.tag.to.hex("(0020,000E)"), charToRaw (tag.dictionary["(0020,000E)","VR"]), 
                               unlist(.value.to.raw (length(serie.UID), 2, TRUE)), serie.UID)
  
  
  L0020[["(0020,0010)"]] <- c(.tag.to.hex("(0020,0010)"), charToRaw(tag.dictionary["(0020,0010)","VR"]),
                              as.raw(c(0x00, 0x00)))
  L0020[["(0020,0011)"]] <- c(.tag.to.hex("(0020,0011)"), charToRaw(tag.dictionary["(0020,0011)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  L0020[["(0020,0013)"]] <- c(.tag.to.hex("(0020,0013)"), charToRaw(tag.dictionary["(0020,0013)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  new.raw <- .ascii.to.hex (paste(obj$xyz0[sop.idx,], collapse="\\"), tag.dictionary["(0020,0032)","VR"])
  L0020[["(0020,0032)"]] <- c(.tag.to.hex("(0020,0032)"), charToRaw(tag.dictionary["(0020,0032)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- .ascii.to.hex (paste(obj$orientation, collapse="\\"), tag.dictionary["(0020,0037)","VR"])
  L0020[["(0020,0037)"]] <- c(.tag.to.hex("(0020,0037)"), charToRaw(tag.dictionary["(0020,0037)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  L0020[["(0020,0052)"]] <- c(.tag.to.hex("(0020,0052)"),charToRaw(tag.dictionary["(0020,0052)","VR"]), frame.of.reference)
  
  L0020[["(0020,1040)"]] <- c(.tag.to.hex("(0020,1040)"), charToRaw(tag.dictionary["(0020,1040)","VR"]),
                              as.raw(c(0x00, 0x00)))
  
  
  ####################################
  tag0028 <- c("(0028,0002)", "(0028,0004)",  "(0028,0008)","(0028,0009)", "(0028,0010)", "(0028,0011)","(0028,0030)", 
               "(0028,0100)", "(0028,0101)", "(0028,0102)","(0028,0103)", "(0028,0106)", 
               "(0028,0107)")#, "(0028,0120)")
  
  L0028 <- lapply(tag0028,function(t) c())
  names(L0028) <- tag0028
  
  new.raw <-as.raw(c(1,0))
  L0028[["(0028,0002)"]] <- c(.tag.to.hex("(0028,0002)"), charToRaw(tag.dictionary["(0028,0002)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  new.raw <- .ascii.to.hex ("MONOCHROME2", tag.dictionary["(0028,0004)","VR"])
  L0028[["(0028,0004)"]] <- c(.tag.to.hex("(0028,0004)"), charToRaw(tag.dictionary["(0028,0004)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <-   .ascii.to.hex (as.character(obj$n.ijk[3]), tag.dictionary["(0028,0008)","VR"])
  L0028[["(0028,0008)"]] <- c(.tag.to.hex("(0028,0008)"), charToRaw(tag.dictionary["(0028,0008)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <-  .tag.to.hex("(3004,000C)")
  L0028[["(0028,0009)"]] <- c(.tag.to.hex("(0028,0009)"), charToRaw(tag.dictionary["(0028,0009)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- c(as.raw(obj$n.ijk[2]%%256) ,as.raw((obj$n.ijk[2]/256)%%256))
  L0028[["(0028,0010)"]] <- c(.tag.to.hex("(0028,0010)"), charToRaw(tag.dictionary["(0028,0010)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  new.raw <-  c(as.raw(obj$n.ijk[1]%%256) ,as.raw((obj$n.ijk[1]/256)%%256))
  L0028[["(0028,0011)"]] <- c(.tag.to.hex("(0028,0011)"), charToRaw(tag.dictionary["(0028,0011)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <- .ascii.to.hex (paste(obj$dxyz[1:2], collapse="\\"), tag.dictionary["(0028,0030)","VR"])
  L0028[["(0028,0030)"]] <- c(.tag.to.hex("(0028,0030)"), charToRaw(tag.dictionary["(0028,0030)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  new.raw <-as.raw(c(16,0))
  L0028[["(0028,0100)"]] <- c(.tag.to.hex("(0028,0100)"), charToRaw(tag.dictionary["(0028,0100)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  L0028[["(0028,0101)"]] <- c(.tag.to.hex("(0028,0101)"), charToRaw(tag.dictionary["(0028,0101)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  new.raw <- as.raw(c(15,0))
  L0028[["(0028,0102)"]] <- c(.tag.to.hex("(0028,0102)"), charToRaw(tag.dictionary["(0028,0102)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  new.raw <-  as.raw(c(0,0))
  if (obj$min.pixel<0) new.raw <-  as.raw(c(1,0))
  L0028[["(0028,0103)"]] <- c(.tag.to.hex("(0028,0103)"), charToRaw(tag.dictionary["(0028,0103)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  #codage de 0.1mGy
  
  vol3D <- as.vector(round(obj$vol3D.data,6))
  slope <- NA
  if (length(vol3D)!=0){
    f <- is.na(vol3D)
    if (all(f)){
      warning(c("obj$vol3D.data contains only NA values. This data will not be exported."))
      vol3D <- NULL
    } else {
      
      if (any(f)) {
        if (is.null(NAvalue)){
          NAvalue <- 0
          warning(c("obj$vol3D.data contains NA values. This data will be replaced by the value 0. ", 
                    "If you wish to change this value, add the argument NAvalue=your_NA_value."))
        }
        vol3D[f] <- NAvalue
      }
      
      slope <- sort(diff(sort(unique(vol3D))))[1]

      obj.max.pixel  <- max(vol3D, na.rm =TRUE)
      obj.min.pixel  <- min(vol3D, na.rm =TRUE)
      max.pixel <- round(obj.max.pixel/slope)
      min.pixel <- round(obj.min.pixel/slope)
      while((abs(max.pixel) > 2^15) |  (abs(min.pixel) > 2^15)){
        slope <- slope*2
        max.pixel <- round(obj.max.pixel/slope)
        min.pixel <- round(obj.min.pixel/slope)
      }
      

      
      new.raw <- c(as.raw(min.pixel%%256) ,as.raw((min.pixel/256)%%256))
      L0028[["(0028,0106)"]] <- c(.tag.to.hex("(0028,0106)"), charToRaw("US"),
                                  as.raw(c(0x02,00)), new.raw)
      
      new.raw <- c(as.raw(max.pixel%%256) ,as.raw((max.pixel/256)%%256))
      L0028[["(0028,0107)"]] <- c(.tag.to.hex("(0028,0107)"), charToRaw("US"),
                                  as.raw(c(0x02,00)), new.raw)
      
      vol3D <- vol3D/slope
      vol3D <- as.vector(rbind(as.raw(vol3D%%256) ,as.raw((vol3D/256)%%256)))
    }
  } else {
    warning(c("There is no obj$vol3D.data in ", obj$object.alias,"."))
  }
  
  
  
  
  ####################################
  
  
  tag3004 <- c("(3004,0002)", "(3004,0004)", "(3004,000A)", "(3004,000C)","(3004,000E)")
  L3004 <- lapply(tag3004,function(t) c())
  names(L3004) <- tag3004
  
  new.raw <- .ascii.to.hex (toupper(obj$unit), tag.dictionary["(3004,0002)","VR"])
  L3004[["(3004,0002)"]] <- c(.tag.to.hex("(3004,0002)"), charToRaw(tag.dictionary["(3004,0002)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  if (is.null(obj$rtdose.info)){
    if (obj$min.pixel<0) {
      type.raw <- .ascii.to.hex ("ERROR", tag.dictionary["(3004,0004)","VR"])
    } else {
      type.raw <- .ascii.to.hex ("PHYSICAL", tag.dictionary["(3004,0004)","VR"])}
    summation.raw <- .ascii.to.hex ("PLAN", tag.dictionary["(3004,000A)","VR"])
    
  } else {
    type.raw <- .ascii.to.hex (obj$rtdose.info$type, tag.dictionary["(3004,0004)","VR"])
    summation.raw <- .ascii.to.hex (obj$rtdose.info$Dose.Summation.Type, tag.dictionary["(3004,000A)","VR"]) 
  }
  
  L3004[["(3004,0004)"]] <- c(.tag.to.hex("(3004,0004)"), charToRaw(tag.dictionary["(3004,0004)","VR"]),
                              unlist(.value.to.raw (length(type.raw),2,TRUE)), type.raw)
  
  L3004[["(3004,000A)"]] <- c(.tag.to.hex("(3004,000A)"), charToRaw(tag.dictionary["(3004,000A)","VR"]),
                              unlist(.value.to.raw (length(summation.raw),2,TRUE)), summation.raw)
  
  new.raw <- .ascii.to.hex (paste(obj$k.idx * obj$dxyz[3], collapse="\\"), tag.dictionary["(3004,000C)","VR"])
  L3004[["(3004,000C)"]] <- c(.tag.to.hex("(3004,000C)"), charToRaw(tag.dictionary["(3004,000C)","VR"]),
                              unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
  
  
  if (length(vol3D)!=0){
    new.raw <- .ascii.to.hex (paste(slope, collapse="\\"), tag.dictionary["(3004,000E)","VR"])
    L3004[["(3004,000E)"]] <- c(.tag.to.hex("(3004,000E)"), charToRaw(tag.dictionary["(3004,000E)","VR"]),
                                unlist(.value.to.raw (length(new.raw),2,TRUE)), new.raw)
    L7FE00010 <- c(.tag.to.hex("(7FE0,0010)"),charToRaw("OW"), as.raw(c(0,0)),
                   unlist(.value.to.raw (length(vol3D),4,TRUE)),vol3D)
  }
  
  
  
  ####################################
  rawdata <- c(preamble,unlist(L0002),unlist(L0008),unlist(L0010),unlist(L0018),
               unlist(L0020),unlist(L0028), unlist(L3004),L300C0060,L7FE00010)
  zz <- file (filen, "wb")
  writeBin (rawdata, zz, size = 1)
  close (zz)
  
  message ("file created : ", basename(filen))
}



################################################################################
.tag.to.hex <- function(tag, Low_Order_First=TRUE){
  hexcode =c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F")
  rtag <- match(unlist(strsplit(tag,""))[c(2:5,7:10)],hexcode)-1
  if (Low_Order_First) return(as.raw(c(16*rtag[3] + rtag[4], 16*rtag[1] + rtag[2],
                                       16*rtag[7] + rtag[8], 16*rtag[5] + rtag[6])))
  return(as.raw(c(16*rtag[1] + rtag[2], 16*rtag[3] + rtag[4],
                  16*rtag[5] + rtag[6], 16*rtag[7] + rtag[8])))
  
}
#' @importFrom sodium hash
.ascii.to.hex <- function(str,VR){
  if (length(str)==0) return(raw(0))
  if (is.na(str)) return(raw(0))
  if (nchar(str)==0) return(raw(0))

  
  switch(VR,
         "AE" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw, as.raw(c(0x0a, 0x0c, 0x0d, 0x1b, 0x5c))))
           if (length(new.raw)>16 | any(m)) new.raw <- raw(0)
         },
         
         "AS" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           if (length(new.raw)!=4 | any(is.na(m)) | !(new.raw[4] %in% charToRaw ("DWMY"))) new.raw <- raw(0)
         },
         
         "CS" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- match(new.raw, as.raw(c(65:90,48:57,92,95,32)))
           # if (length(new.raw)>16 | any(is.na(m))) new.raw <- raw(0)
           if (any(is.na(m))) new.raw <- raw(0)
         },
         "DA" = {
           new.raw <- charToRaw (str)
           m <- match(new.raw[1:8], as.raw(48:57))
           if (length(new.raw)!=8 | any(is.na(m))) new.raw <- raw(0)
         },
         
         "DS" = {
           new.raw <- charToRaw (str)
           m <- match(new.raw,as.raw(0x20))
           new.raw <- c(new.raw[is.na(m)],as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- match(new.raw, charToRaw("0123456789+-Ee. \\"))
           if (any(is.na(m))) new.raw <- raw(0)
         },
         
         "DT" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- match(new.raw, charToRaw("0123456789+-. "))
           if (length(new.raw)>26 | any(is.na(m))) new.raw <- raw(0)
         },
         
         "IS" = {
           mem <-unlist(strsplit(str,"[\\]"))
           mem <- suppressWarnings(as.numeric(mem))
           if (any(is.na(mem))){
             new.raw <- raw(0)
           } else if (any(mem < -2^31) | any(mem > +2^31 - 1 )){
             new.raw <- raw(0)
           } else {
             new.raw <- charToRaw (str)
             m <- match(new.raw,as.raw(0x20))
             new.raw <- c(new.raw[is.na(m)],as.raw(0x20))
             new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
             m <- match(new.raw, charToRaw("\\0123456789+- "))
             if (length(new.raw)>12 | any(is.na(m))) new.raw <- raw(0)
           }
         },
         
         "LO" = {
           new.raw <- c(charToRaw (trimws(str, which="right")),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw,as.raw(c(0:26,28:31,92))))
           if (length(new.raw)>64 | any(m)) new.raw <- raw(0)
         },
         
         "LT" = {
           new.raw <- c(charToRaw (trimws(str, which="right")),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
           
           if (length(new.raw)>64 | any(m)) new.raw <- raw(0)
         }, 
         
         "PN" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           str_ <- rawToChar(new.raw)
           le <- 0
           if (nchar(str_)>0) le <- max(nchar(unlist(strsplit(str_,"^", fixed=TRUE))))
           m <- !is.na(match(new.raw, as.raw(c(10, 12, 13, 92))))
           if (le>64 | any(m)) new.raw <- raw(0)
         },
         
         "SH" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw, as.raw(c(0:26,28:31,92))))
           if (length(new.raw)>16 | any(m)) new.raw <- raw(0)
         },
         
         "ST" = {
           new.raw <- c(charToRaw (trimws(str, which="right")),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
           if (length(new.raw)>1024 | any(m)) new.raw <- raw(0)
         },
         
         "TM" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- match(new.raw, charToRaw("0123456789. "))
           if (length(new.raw)>16 | any(is.na(m))) new.raw <- raw(0)
         }, 
         
         "UI" = {
           new.raw <- c(charToRaw (trimws(str)),as.raw(0x00))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- match(new.raw, as.raw(c(0,46,48:57)))
           if (length(new.raw)>64 | any(is.na(m))) new.raw <- raw(0)
         }, 
         
         "UN" = {
           new.raw <- c(charToRaw (trimws(str, which="right")),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
         }, 
         
         "UT"= {
           new.raw <- c(charToRaw (trimws(str, which="right")),as.raw(0x20))
           new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
           m <- !is.na(match(new.raw,as.raw(c(0:9,11,14:26,28:31))))
           if (length(new.raw)>2^32-2 | any(m))new.raw <- raw(0)
         })
  if (length(new.raw)==0) warning(str," is not compliant with ", VR," DICOM representation" )
  return(new.raw)
}

get.obj.from.alias <- function (obj.alias, patient){
  if (is.null(obj.alias))  return(list())
  obj.alias <- obj.alias[obj.alias!=""]
  if (length(obj.alias)==0) return(list())
  
  start.idx <- which(names(patient)=="T.MAT") + 1
  end.idx <- length(patient)
  n <-names(patient)
  tab <- do.call(rbind,lapply(start.idx:end.idx, function(idx) 
    data.frame(alias = names(patient[[idx]]), modality = n[idx])))
  L <-lapply(obj.alias, function(n){
    f <- tab$alias==n
    if (any(f)) return(patient[[tab$modality[f]]][[n]])
    return(NULL)
  })
  names(L) <- obj.alias
  L
}


create.SOPUID <- function(obj,size=64){
  
  obj$ref.object.alias <- NULL
  obj$object.info <- NULL
  obj$ref.object.info <- NULL
  
  if (!(is(obj,"volume") | is(obj,"struct"))) return(NULL)
  if  (!(obj$modality %in% c("ct","mr","rtdose","rtstruct"))) return(NULL)
  
  obj$ref.object.info <- NULL
  obj$object.info <- NULL
  le.esp <- nchar(.espadon.UID()) + 1
  espadon.UID <- charToRaw(paste0(.espadon.UID(),"."))
  
  nb_cut <- 1
  
  nc_nbcut <- 0
  raw.mat <- NULL
  if (obj$modality %in% c("ct","mr")) {
    nb_cut <- max(obj$k.idx)
    nc_nbcut <- nchar(as.character(nb_cut))
    prefix <-  rawToChar(rep(as.raw(0x30),nc_nbcut-1))
    nb.str <- paste0(prefix,obj$k.idx)
    le.vect <- cumsum(nchar(nb.str))
    raw.mat  <- charToRaw(paste(nb.str,collapse=""))
    raw.mat <- do.call(cbind,lapply(nc_nbcut:1, function(i) raw.mat[le.vect-i+1]))
    nb_cut <- nrow(raw.mat)
  }
  
  UID.base <- matrix(create.UID(obj, size = floor((size-nc_nbcut)/2)*2),nrow=1, byrow = TRUE)
  
  l1155 <- cbind( UID.base [rep(1,nb_cut),,drop=FALSE], raw.mat, as.raw(0))
  
  l1155 <- l1155[,1:(floor(ncol(l1155)/2)*2), drop=FALSE]
  l1155 <- lapply(1:nrow(l1155),function(i) l1155[i,])
  
  l1150 <- switch (obj$modality,
                   "ct"= "1.2.840.10008.5.1.4.1.1.2",
                   "mr"= "1.2.840.10008.5.1.4.1.1.4",
                   "rtdose" = "1.2.840.10008.5.1.4.1.1.481.2",
                   "rtstruct" = "1.2.840.10008.5.1.4.1.1.481.3")
  
  # l1150 <-  matrix(rep(.ascii.to.hex(l1150,"UI") ,nb_cut), nrow= nb_cut, byrow = TRUE)
  l1150 <- .ascii.to.hex(l1150,"UI")
  l1150 <- lapply(l1155, function(l)l1150)
  return(list(l1150=l1150,l1155=l1155))     
  
}
#' @importFrom qs2 qs_serialize
create.UID <- function(obj, el.name=NULL,supp.info.l=NULL, size=64){ #>44
  
  obj$ref.object.alias <- NULL
  obj$object.info <- NULL
  obj$ref.object.info <- NULL
  
  
  if (is.null(el.name)){
    UID <- qs_serialize(list(obj,supp.info.l))
  }else {
    UID <- qs_serialize(list(obj[el.name],supp.info.l))
  }
  le.esp <- nchar(.espadon.UID()) + 1
  espadon.UID <- charToRaw(paste0(.espadon.UID(),"."))
  UID <- paste(.espadon.UID(),
               paste(as.numeric(hash(UID, key = NULL,  size = max(size,16+le.esp)-le.esp)) %% 10, 
                     collapse=""), sep=".")
  new.raw <- c(charToRaw (UID),as.raw(0x00))
  new.raw <- new.raw[1:(2 * floor(length(new.raw)/2))]
  return(new.raw)
}


get.obj.from.alias <- function (obj.alias, patient){
  if (is.null(obj.alias))  return(list())
  obj.alias <- obj.alias[obj.alias!=""]
  if (length(obj.alias)==0) return(list())
  
  start.idx <- which(names(patient)=="T.MAT") + 1
  end.idx <- length(patient)
  n <-names(patient)
  tab <- do.call(rbind,lapply(start.idx:end.idx, function(idx) 
    data.frame(alias = names(patient[[idx]]), modality = n[idx])))
  L <-lapply(obj.alias, function(n){
    f <- tab$alias==n
    if (any(f)) return(patient[[tab$modality[f]]][[n]])
    return(NULL)
  })
  names(L) <- obj.alias
  L
}


