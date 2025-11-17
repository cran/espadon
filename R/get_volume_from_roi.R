#' Volume of a region of interest (RoI)
#' @description The \code{get.volume.from.roi} function extracts the volume
#' in\mjeqn{cm^3}{ascii} of one or more RoI, from the \code{$roi.info} of the 
#' "struct" class object.
#' @param struct "struct" class object.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} 
#' object. By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @return Returns a vector of the volumes in \mjeqn{cm^3}{ascii} of the requested 
#' RoI.
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all set 
#' to NULL, all RoI are selected.
#' @seealso \link[espadon]{get.volume.from.bin}, \link[espadon]{select.names}
#' @examples
#' # loading of toy-patient objects
#' step <- 4
#' patient <- toy.load.patient (modality = c ("rtstruct"), 
#'                              dxyz = rep (step, 3))
#' S <- patient$rtstruct[[1]] 
#'
#' # Volume extraction
#' vol <- get.volume.from.roi (S, roi.sname = "bra", roi.idx = c (1, 3))
#' names (vol)
#' vol

#' @export
#' @importFrom methods is
get.volume.from.roi <- function (struct, roi.name = NULL, roi.sname = NULL, 
                                 roi.idx = NULL) {
  
  # if (length (roi.idx) != 1) {
  #   warning ("multiple names or no names forbidden.")
  #   return (NULL)
  # }
  if (!is (struct, "struct")) 
    stop ("struct should be a struct class object.")
 
  if(is.null(struct$roi.data))
    stop ("empty roi.data\n")

  
  roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  
  v <- struct$roi.info$vol[roi.idx]
  names(v) <- struct$roi.info$roi.pseudo[roi.idx]
  return (v)
}