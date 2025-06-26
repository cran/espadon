################################################################################
#' Save a T.MAT class object
#' @description The \code{save.T.MAT} function saves the data required by 
#' \link[espadon]{load.T.MAT}, \link[espadon]{load.patient.from.dicom} or 
#' \link[espadon]{load.patient.from.Rdcm} to generate \code{T.MAT}, as 
#'  pre-formatted Rdcm files.
#' @param T.MAT "t.mat" class object to save.
#' @param dirname Directory where new reg .Rdcm files will be saved.
#' @return Returns \code{TRUE}, if all reg files generating T.MAT are saved.
#' @details Reg files from DICOM files cannot be updated with the \code{save.T.MAT}
#' function. Only transfer matrices added with \link[espadon]{ref.add} or 
#' \link[espadon]{ref.cutplane.add} will be saved.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#'                              
#' # Save T.MAT to a temporary file pat.dir
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' save.T.MAT (patient$T.MAT, dirname = pat.dir)
#' list.files(pat.dir)
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @importFrom qs2 qs_serialize
#' @importFrom methods is
#' @export
save.T.MAT <- function (T.MAT, dirname) {
  if (!is (T.MAT, "t.mat")) stop("T.MAT should be a t.mat class object.")
  if (dirname==""){dirname="."}
  if (!dir.exists(dirname)) {
    if (!dir.create(dirname,recursive = T))
      stop(paste(dirname,"doesn't exist and can't be created"))
  }

  if (nrow(T.MAT$reg.info$file)==0) {
    message("no T.MAT reg to save")
    return(TRUE)
  }
  ok <- rep(FALSE, nrow(T.MAT$reg.info$file))
  for (idx in 1:nrow(T.MAT$reg.info$file)){
    m <- T.MAT$matrix.list [[T.MAT$reg.info$file$t[idx]]]
    tab <- T.MAT$matrix.description[T.MAT$matrix.description$t==
                                      T.MAT$reg.info$file$t[idx],]
    m.name <- as.character (T.MAT$ref.info[T.MAT$ref.info$ref.pseudo == as.character(tab$src),"ref"])

    if (T.MAT$reg.info$file$path[idx]!="local") {
      obj <- load.Rdcm.raw.data(T.MAT$reg.info$file$path[idx], data=FALSE, address=FALSE)
      if ((!obj$from.dcm) & !all(obj$header$ref.data[[m.name]]$matrix==m)){
        obj$header$ref.data[[m.name]]$matrix <- m
      } else {
        obj <- NULL
      }
      
    } else {
      reg <- list ()
      reg$patient <- T.MAT$reg.info$patient[1,1]
      reg$patient.name <- T.MAT$reg.info$patient[1,2]
      reg$patient.bd <- T.MAT$reg.info$patient[1,3]
      reg$patient.sex <- T.MAT$reg.info$patient[1,4]
      reg$file.basename <- ""
      reg$file.dirname <- ""
      reg$object.name <- ""
      reg$object.alias <- ""
      reg$frame.of.reference <-  as.character (T.MAT$ref.info[T.MAT$ref.info$ref.pseudo == as.character(tab$dest),"ref"])
      reg$ref.pseudo <- as.character(tab$dest)
      reg$modality <- "reg"
      reg$description <- gsub("<-"," from ", tab$t)
      reg$creation.date <- format(Sys.Date(), "%Y%m%d")
      reg$object.name <- gsub("<-","_from_", tab$t)
      reg$object.alias <- reg$object.name
      reg$file.basename <- paste0(reg$object.name,".Rdcm")
      reg$nb.of.ref <- 2
      reg$ref.data <- list()
      
      reg$ref.data[[reg$frame.of.reference]] <- list (
        src =reg$frame.of.reference,
        type = "RIGID ",
        matrix = diag (4))
      
      reg$ref.data[[m.name]] <- list (
        src =m.name,
        type = "RIGID ",
        matrix = m)
      
      
      if (dirname=="") {dirname="."}
      T.MAT$reg.info$file$path[idx] <- file.path(dirname,reg$file.basename )
      obj<- list()
      obj$header <- reg
      
    }
    
    if (is.null(obj)) {
      ok[idx] <- TRUE
    }else {
      if (file.exists(T.MAT$reg.info$file$path[idx])) warning (paste(reg$file.basename, "already exists."))
      h <- qs_serialize(c(obj$header,list(espadon.version=.espadon.version())))
      a <- numeric(0)
      d <- qs_serialize(obj$data)
      zz <-  file(T.MAT$reg.info$file$path[idx], "wb")
      writeBin(length(h),zz,size=4, endian="little")
      writeBin(length(a),zz,size=4, endian="little")
      writeBin(length(d),zz,size=4, endian="little")
      writeBin(h,zz,endian="little")
      writeBin(a,zz,endian="little")
      writeBin(d,zz,endian="little")
      close(zz)
      ok[idx] <- file.exists(T.MAT$reg.info$file$path[idx])
    }
  }
  
  return(all(ok))
}