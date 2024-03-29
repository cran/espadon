% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vol_sum.R
\name{vol.sum}
\alias{vol.sum}
\title{Sum of 2 volumes}
\usage{
vol.sum(vol1, vol2, alias = "", description = NULL)
}
\arguments{
\item{vol1, vol2}{"volume" class objects. The
2 volumes must have the same modality, and the same grid (i.e. share the same
position of the voxels).}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL} (default value), it will be set to
\code{paste (vol1$object.alias, "+", vol2$object.alias)}.}
}
\value{
Returns a "volume" class object (see \link[espadon]{espadon.class}
for class definitions), with the same grid and modality
as \code{vol1} and \code{vol2}, sum of \code{vol1} and \code{vol2}.
}
\description{
The \code{vol.sum} function adds two "volume" class objects
of the same grid and of the same modality.
}
\examples{
# loading of toy-patient objects (decrease dxyz and increase beam.nb for 
# better result)
step <- 5
pat<- toy.load.patient (modality = c ( "rtdose"), dxyz = rep (step, 3),
                        beam.nb = 3)

# Double dose
D <- vol.sum (pat$rtdose[[1]], pat$rtdose[[1]])
pat$rtdose[[1]]$max.pixel
D$max.pixel
}
