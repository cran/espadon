#' Loading an \pkg{espadon} object from DICOM files or folder
#' @description Loading an \pkg{espadon} object from DICOM files or folder.
#' @param dcm.files String vector, representing the list of the full names of the  
#' DICOM files of the same DICOM object, or its directory.
#' @param data Boolean. Only valid for objects usable by the \pkg{espadon} package, 
#' namely ct, mr, rtdose, rtstruct, pt... If \code{data = TRUE}, either the values 
#' of the voxels when modality is (ct, mr, rtdose, pt), or the coordinates of the 
#' RoI when modality is rtstruct, are loaded into memory.
#' @param ref.pseudo String, \code{$ref.pseudo} (i.e. pseudonym of the frame of 
#' reference) 
#' to assign to the loaded object.
#' @param tag.dictionary Dataframe, by default equal to 
#' \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.
#' @param verbose Boolean. If \code{TRUE}, a progress bar indicates the progress 
#' of the conversion.
#' @return Returns an \pkg{espadon} object of class "dvh","histo","histo2D","mesh", 
#' "rtplan","struct", "undef" or "volume" depending on the object modality. See 
#' \link[espadon]{espadon.class} for class definitions.
#' @seealso \link[espadon]{load.obj.data} and \link[espadon]{load.obj.from.Rdcm}


#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file pat.dir for testing.
#' pat.dir <- file.path (tempdir(), "PM_dcm") 
#' dir.create (pat.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = pat.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#'
#' # loading of rt-plan object
#' RTplan <- load.obj.from.dicom (dcm.filename)
#' str (RTplan)
                     
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)

#' @export
load.obj.from.dicom <- function (dcm.files, data = TRUE, ref.pseudo = "ref1",
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
  
  dicomlist <- .load.dcm(dcm.filenames, data=data, verbose = verbose, 
                         tag.dictionary = tag.dictionary)
  if (length(dicomlist)==0) return (NULL)
  if (length(dicomlist)>1) {
    message ("More than one object is found. Only first one is loaded.")
  }
  obj <- .load.object(dicomlist[[1]], data=data)
  obj$ref.pseudo <- ref.pseudo
 return (obj)
}
