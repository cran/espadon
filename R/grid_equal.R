#' Comparison of the grids of two volume objects
#' @description The \code{grid.equal} function checks that two volumes share the 
#' same grid, i.e. the same frame of reference, the same origin point, and the 
#' same dx, dy, dz steps.
#' @param vol1,vol2 "volume" class objects
#' @return Returns \code{TRUE} if the 2 volumes share the same grid.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (modality = c ("ct","mr","rtdose"), roi.name = "", 
#'                              dxyz = c (4, 4, 4), beam.nb = 1)
#'
#' # Comparison of the grids
#' grid.equal (patient$rtdose[[1]], patient$ct[[1]])
#' grid.equal (patient$mr[[1]], patient$ct[[1]])

#' @export
#' @importFrom methods is
grid.equal <- function (vol1, vol2)  {
  if (!is (vol1, "volume") | !is (vol2, "volume")) return (FALSE)
  if (vol1$ref.pseudo != vol2$ref.pseudo) return (FALSE)
  if (!all (abs(vol1$xyz.from.ijk - vol2$xyz.from.ijk) < 1e-6)) return (FALSE)
  if (!all (vol1$n.ijk == vol2$n.ijk)) return (FALSE)
  return (TRUE)
}