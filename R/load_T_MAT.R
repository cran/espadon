####################################################################################
#' Loading of information about transfer matrices between frames of reference of 
#' patient objects.
#' @description The \code{load.T.MAT} function lists all the frames of reference
#' of the objects included in the patient directory. It concatenates all the 
#' information of the reg matrices of a directory, creating, among  other things, 
#' a list of 4x4 transfer matrices between frames of reference.
#' @param dirname Directory containing a patient's DICOM objects.
#' @return Returns a "t.mat" class object. It is a list that includes :
#' \itemize{
#' \item \code{$ref.info}: dataframe giving the correspondence between the frame of
#' reference (column \code{$ref}) of the DICOM object (TAG (0020,0052) ) and its
#' pseudonym (column \code{$ref_pseudo}).
#' \item \code{$reg.info}:list of dataframes : the first one gives the PID, 
#' birthday, and sex of the patient, the second one gives the name of the source 
#' file of transfer matrices.
#' \item \code{$matrix.description}: dataframe giving the transfer matrix names
#' (column \code{$t}), its source frame of reference (column \code{$src}), the 
#' destination frame of reference (column \code{$dest}), and its type (\code{$type}). 
#' Note: only the RIGID type is supported.
#' \item \code{$matrix.list}: list of 4X4 transfer matrices. This list contains
#' at least as many Identity matrices as there are \code{ref.pseudo}.
#' }
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
#' T.MAT <- load.T.MAT (pat.dir)
#' T.MAT
#' 
#' # Cleaning  temporary directory
#' unlink (pat.dir, recursive = TRUE)
 
#' @export
load.T.MAT <- function (dirname) {
  lf <- list.files (dirname, pattern ="[.]Rdcm", full.names = T)
  if (length(lf)==0) return(NULL)
  reg.list <-lapply (lf, function (fn){
    h <- load.Rdcm.raw.data (fn, data=FALSE, address=FALSE)$header
    if (h$modality!="reg") return(list(ref.pseudo =  h$ref.pseudo,
                                       ref = h$frame.of.reference, 
                                       patient =h$patient,
                                       patient.bd =h$patient.bd,
                                       patient.sex=h$patient.sex))
    return(list(ref.pseudo =  h$ref.pseudo,
                ref = h$frame.of.reference, 
                patient =h$patient,
                patient.bd =h$patient.bd,
                patient.sex=h$patient.sex, reg.list=h))
  })
  .load.T.MAT.by.reglist (reg.list)
}

