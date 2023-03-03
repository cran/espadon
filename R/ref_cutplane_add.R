#' Adding volume's cutting planes frame of reference in T.MAT
#' @description The \code{ref.cutplane.add} function adds in T.MAT the transfer 
#' matrices from or to volume's cutting planes frame of reference.
#' @param vol "volume" class object.
#' @param origin Vector of the x, y, z coordinates of the origin point of 
#' the cut planes frame of reference.
#' @param ref.cutplane Name of the volume's cutting planes frame of reference. 
#' By default \code{ref.cutplane = paste0 (vol$ref.pseudo,".m")}.
#' @param T.MAT "t.mat" class object created by 
#' \link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, then only the link 
#' between \code{vol$ref.pseudo} and \code{ref.name} is established.
#' @return Returns a "t.mat" class object, which contains the transfer matrices 
#' from or to volume's cutting planes frame of reference. If the \code{T.MAT} is \code{NULL}, 
#' then the returned object will contain only 4 matrices: "src.ref<-src.ref",
#' "src.ref<-ref.cutplane", "ref.cutplane<-ref.cutplane", "ref.cutplane<-src.ref".
#' @seealso \link[espadon]{ref.add}, \link[espadon]{ref.srctodest.add}, 
#' \link[espadon]{ref.remove}.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = "mr", roi.name = "", dxyz = c (4, 4, 4))
#' MR <- patient$mr[[1]]
#' MR$xyz.from.ijk
#'
#' # creation of t.mat, containing the transfer matrix to the frame of reference 
#' # of the MR cutting planes
#' t.mat <- ref.cutplane.add (MR)
#' 
#' # Change of frame of reference
#' MR.m <- vol.in.new.ref (MR, paste0 (MR$ref.pseudo, "m"), t.mat)
#' 
#' MR.m$xyz.from.ijk

#' @export
#' @importFrom methods is
ref.cutplane.add <- function (vol, origin = c(0,0,0),
                              ref.cutplane = paste0 (vol$ref.pseudo, "m"),
                              T.MAT = NULL){
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  t.mat <- ( ref.add( src.ref= vol$ref.pseudo, 
                    orientation = vol$orientation, 
                    origin = origin, 
                    new.ref.pseudo = ref.cutplane,
                    T.MAT = T.MAT))
  t.mat$ref.info[t.mat$ref.info$ref.pseudo==vol$ref.pseudo,"ref"] <- vol$frame.of.reference
  t.mat$reg.info$patient <- unique(rbind(t.mat$reg.info$patient,
                                         data.frame(patient=vol$patient,
                                                    patient.name=vol$patient.name,
                                                    patient.bd=vol$patient.bd,
                                                    patient.sex=vol$patient.sex)
                                   ))
  
  return(t.mat)
}