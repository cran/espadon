#' @keywords internal
"_PACKAGE"
#' @author Cathy & Jean-Marc Fontbonne
#' @docType package
#' @name espadon-package
#' @aliases espadon
#' @title Easy Study of Patient DICOM Data in Oncology
#' @description \pkg{espadon} works in a native way (user friendly):
#' * on images (CT, MR, PT) and fully manages changes of referential (REG)
#' * on dosimetry (rt-dose)
#' * on structures (rt-struct)
#' @description It is also able to use any DICOM format file, as long as the user knows where to look for the information.
#' In addition to the simplified use of the above-mentioned formats, espadon contains many functions that allow the user to rework documents to produce new features that can be integrated into machine learning.
#' @description \strong{\pkg{espadon} integrates functionalities:}
#' - file loading, patient-centered information fusion
#' - handling of 3D images
#' \itemize{
#' \item changes of referential
#' \item resampling, filtering
#' }
#' - Contour manipulation, allowing the creation of new contours
#' \itemize{
#' \item fusion, intersection, inversion, erosion, dilation, opening, closing
#' \item segmentation based on the image
#' }
#' - 2D representation
#' \itemize{
#' \item with application of masks
#' \item and transport of structures in the new reference frames
#' }
#' - ... and 3D
#' \itemize{
#' \item production of mesh for the representation
#' }
#' - for measurements
#' \itemize{
#' \item calculation of histograms (1D, 2D), DVH (integrating Monte-Carlo to simulate organ movements)
#' \item Measurement of surfaces, radii of curvature (on the mesh) and volumes
#' }
#' - measurement of standard dosimetry indicators
#'
#' @description In the calculation flow, various objects are created and accessible to the developer for his own use.
#'

#' @useDynLib espadon
#' @importFrom Rcpp evalCpp
NULL