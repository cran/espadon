#' Creation of orientation
#' @description The \code{orientation.create} function creates the orientation 
#' vectors of a plane:
#' \itemize{
#' \item from 3 points A, B and C (see details),
#' \item or from 2 vectors B and C, resp. defining x and y-axis (see details),
#' \item or from 2 points A, B defining x-axis, and the normal vector to the plane (see details), 
#' \item or from a vector B defining x-axis, and the normal vector to the plane (see details). 
#' } 
#' @param A Vector of the x, y and z coordinates of point \code{A}, by default 
#' equal to \code{c(0,0,0)} in the case where \code{B} and \code{C} are vectors.
#' @param B Vector of x, y and z coordinates of point or vector \code{B}.
#' @param C Vector of x, y and z coordinates of point or vector \code{C}.
#' @param normal Vector of x, y and z coordinates of normal vector.
#' @return Returns the orientation of the plane. 
#' That means the concatenation of 2 vectors, defining an orthonormal basis of 
#' the plane.

#' @details  When using \code{B} and \code{C}, \code{B-A} define the x-axis 
#' unit vector. The unit vector of the y-axis is orthonormal to the x-axis, coplanar 
#' with \code{A}, \code{B} and \code{C}, and in the direction of \code{A} to \code{C}.
#' @details  When using \code{B} and \code{normal}, the unit vector of the x-axis is 
#' orthonormal to the \code{normal} vector, in the direction of \code{A} to \code{B}.
#' The unit vector of the y-axis is defined so as to constitute a direct orthonormal 
#' basis with the unit vector of the x-axis and the normal vector of the plane.

#' @export
#'
#' @examples
#' A <- c (-29.93, 18.85, 4.34)
#' B <- c (28.73, 15.36, 4.46)
#' C <- c (1.53, 75.21, 13.51)
#' orientation.create (A, B, C)

orientation.create  <- function(A = c (0, 0, 0), B= NULL, C = NULL, 
                                normal = NULL){
  

  if (is.null(B) | (is.null(C) & is.null (normal))) stop ("Allocate B and C or normal.")
  u <- as.numeric(B - A)
  
  if (!is.null(C)) {
    u <- u / sqrt(as.numeric(u%*%u))
    v <- as.numeric (C - A)
    k <- as.numeric (v %*% u)
    v <- v- k*u
    n <- sqrt (as.numeric (v%*%v))
    if(round(n,6)==0) stop ("The 3 points are aligned.")
    v <- v / n
   
  } else if (!is.null (normal)){
    normal <-  normal / sqrt (as.numeric (normal%*%normal))
    k <- as.numeric (u %*% normal)
    u <- u - k * normal
    n <- sqrt(as.numeric(u%*%u))
    if(round(n,6)==0) stop ("B-A and normal are collinear.")
    u <- u / n
    v <-  vector.product(normal,u)
    n <- sqrt (as.numeric (v%*%v))
    v <- v / n
  }
  
  return (c(u,v))
 
}