#' Volume selected by binary volume
#' \loadmathjax
#' @description The \code{get.volume.from.bin} function calculates the volume in
#' \mjeqn{cm^3}{ascii} of the selection specified by a "volume" class object 
#' of \code{"binary"} modality.
#' @param bin "volume" class object, of "binary" modality.
#' @return Returns the volume of the binary selection, in \mjeqn{cm^3}{ascii}.
#' @seealso \link[espadon]{get.volume.from.roi}
#' @examples
#' # loading of toy-patient objects
#' step <- 4
#' patient <- toy.load.patient (modality = c ("ct", "rtstruct"), roi.name = "brain",
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]] 
#' 
#' # creation of a binary object
#' bin.brain <- bin.from.roi (vol = CT, struct = S, roi.sname = "bra", 
#'                            verbose = FALSE)
#' # Volume calculation
#' get.volume.from.bin (bin.brain)

#' @export
#' @importFrom methods is
get.volume.from.bin <- function (bin) {
  if (!is (bin, "volume")) 
    stop ("bin should be a volume class object.")
  if ((bin$modality!="binary")) 
    stop ("bin must be modality binary.")

  
  return (sum (bin$vol3D.data, na.rm = TRUE) * abs(prod (bin$dxyz))/1000)
}
