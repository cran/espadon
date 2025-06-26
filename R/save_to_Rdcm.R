#' Save a \pkg{espadon} object in a pre-formatted *.Rdcm file
#' @description The function \code{save.to.Rdcm} allows you to save an object 
#' created by \pkg{espadon} in a pre-formatted *.Rdcm file. This object will also be 
#' accessible by the \code{load.patient.from.Rdcm} function.
#' @param obj \pkg{espadon} object of class \code{"volume"}, \code{"struct"}, \code{"mesh"},
#' \code{"histo"}, \code{"dvh"}, \code{"histo2D"}.
#' @param object.name Character string, representing the name of the object, 
#' default to \code{obj$object.alias}.
#' @param dirname Directory where new files from \code{obj} will be saved.
#' @note \code{save.to.Rdcm} can not replace an *.Rdcm file created by 
#' \link[espadon]{dicom.to.Rdcm.converter}.
#' @return Returns \code{TRUE}, if \code{paste0(object.name,".Rdcm")} exists in 
#' \code{dirname}.
#' @return Returns \code{FALSE}, if \code{object.name} is not a valid file name, 
#' or if the file that is created would replace a *.Rdcm file created by 
#' \link[espadon]{dicom.to.Rdcm.converter}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#'                              
#' # Save T.MAT to a temporary file pat.dir
#' pat.dir <- file.path (tempdir(), "PM_Rdcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
#' save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
#' list.files(pat.dir)
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @importFrom qs2 qs_serialize
#' @export
save.to.Rdcm <- function (obj, object.name = obj$object.alias, dirname = obj$file.dirname) {
  
  
  obj$object.name <- gsub( '[[:blank:]]','_',paste(gsub('[[:punct:]]+', '', unlist(strsplit(object.name,"_"))),collapse="_"))
  if (obj$object.name=="") {
    warning ("specify a valid object.name")
    return (FALSE)
  }
  
  obj$file.basename <- paste0(obj$object.name,".Rdcm")
  
  if (dirname==""){dirname="."}
  if (!dir.exists(dirname)) {
    if (!dir.create(dirname,recursive = T))
      stop(paste(dirname,"doesn't exist and can't be created"))
    }

  Rdcm.filename <- file.path(dirname,obj$file.basename )
  
  if (file.exists(Rdcm.filename)){
    old_obj <- load.Rdcm.raw.data (Rdcm.filename, data = FALSE)
    if (old_obj$from.dcm){
      warning ("this object already exists and comes from dicom files. It cannot be saved under this name.")
      return (FALSE)
    } else {warning ("this object already exists.")
    }
  }

  if (!class(obj) %in% c("volume", "struct", "mesh", "histo", "dvh", "histo2D")){
    warning ("not supported class.")
    return(FALSE)
  }
  obj_ <- list()
  
  switch (class(obj),
          "volume" = {
            obj_$header <- obj
            obj_$data <- obj$vol3D.data
            obj_$header$vol3D.data <- NULL},
          
          "struct" = {
            obj_$header <- obj
            obj_$header$roi.info <- obj$roi.info[,1:6];
            obj_$data <- list(roi.info = obj$roi.info,  roi.data = obj$roi.data)
            obj_$header$roi.data <- NULL},
          
          "mesh"= {
            obj_$header <- obj
            obj_$data <- obj$mesh
            obj_$header$mesh <- NULL},
          
          "histo"= {
            obj_$header <- obj
            idx <- grep("^MC[.]", names(obj))
            if (length(idx)>0){
              obj_$data <- obj[idx]
            }
            obj_$header[idx] <- NULL},
          
          "dvh"= {
            obj_$header <- obj
            idx <- grep("^MC[.]", names(obj))
            if (length(idx)>0){
              obj_$data <- obj[idx]
            }
            obj_$header[idx] <- NULL},
          
          "histo2D"= {
            obj_$header <- obj
            obj_$data <- obj$density.map
            obj_$header$density.map <- NULL},
          
          obj_$header <- obj
  )
  h <- qs_serialize(c(obj_$header,list(espadon.version=.espadon.version())))
  a <- numeric(0) #qserialize(obj_$address)
  d <- qs_serialize(obj_$data)
  zz <-  file(Rdcm.filename, "wb")
  writeBin(length(h),zz,size=4, endian="little")
  writeBin(length(a),zz,size=4, endian="little")
  writeBin(length(d),zz,size=4, endian="little")
  writeBin(h,zz,endian="little")
  writeBin(a,zz,endian="little")
  writeBin(d,zz,endian="little")
  close(zz)
  return (file.exists(Rdcm.filename))
}