
#' Vector product of two vectors
#' @param v1  Vector of x, y, z coordinates
#' @param v2  Vector of x, y, z coordinates
#'
#' @return Returns the x, y, z coordinates of the vector product of v1 and v2

#' @examples
#' vector.product(c (1, 0, 0), c (0, 1, 0))
#'
#' @export
vector.product <- function (v1,v2){
  v1_<- c(v1[2:length(v2)],v1[1])
  v2_<- c(v2[2:length(v2)],v2[1])
  v <-  v1 *v2_ - v2*v1_
  as.numeric(c(v[2:length(v)],v[1]))
}