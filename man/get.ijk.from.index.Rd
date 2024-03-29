% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ijk_from_index.R
\name{get.ijk.from.index}
\alias{get.ijk.from.index}
\title{Conversion of the indices of a point into ijk vector}
\usage{
get.ijk.from.index(idx, vol)
}
\arguments{
\item{idx}{Index, or matrix of voxel indices of the array \code{vol$vol3D.data}.}

\item{vol}{"volume" class object.}
}
\value{
Returns an i, j, k column matrix of the DICOM indices of the points
of \code{vol$vol3D.data}.
}
\description{
The \code{get.ijk.from.index} function converts the voxel indices
of \code{vol$vol3D.data} (for example, obtained with the function \code{which})
into a vector or matrix of DICOM indices i, j, k.
}
\examples{
# loading of toy-patient objects (decrease dxyz and increase beam.nb for 
# better result)
step <- 4
patient <- toy.load.patient (modality = "rtdose", roi.name = "", 
                             dxyz = rep (step, 3), beam.nb = 3)
D <- patient$rtdose[[1]]

# voxels location where the dose is greater than 99.9\% Dmax
Dmax <- max (D$vol3D.data, na.rm = TRUE) # D$max.pixel
get.ijk.from.index (which (D$vol3D.data >= 0.999 * Dmax), D)
# or
get.ijk.from.index (which (D$vol3D.data >= 0.999 * Dmax, arr.ind = TRUE), D)

ijk <- as.numeric (get.ijk.from.index (which.max (D$vol3D.data), D))
display.kplane (D, k = ijk[3]) 
}
\seealso{
\link[espadon]{get.value.from.ijk}, \link[espadon]{display.kplane}
}
