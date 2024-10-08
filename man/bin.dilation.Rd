% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_dilation.R
\name{bin.dilation}
\alias{bin.dilation}
\title{Binary volume dilation}
\usage{
bin.dilation(vol, radius = 10, alias = "", description = NULL)
}
\arguments{
\item{vol}{"volume" class object, of \code{"binary"} modality}

\item{radius}{Positive number, or xyz-vector of 3 positive numbers.By default,
radius = 10.}

\item{alias}{Character string, \code{$object.alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL}
(default value), it will be set to \code{paste (vol$object.alias, "dilataion r =", radius)}.}
}
\value{
Returns a "volume" class object of \code{"binary"} modality (see
\link[espadon]{espadon.class} for class definitions), with
the same grid as \code{vol}, in which the selected volume has been enlarged
by the \code{radius}.
}
\description{
The \code{bin.dilation} function enlarges a
"volume" class object, of \code{"binary"} modality, by means of
convolution with a sphere.
Dilation is useful for :
\itemize{
\item filling holes that are smaller than the \code{radius},
\item enlarging capes,
\item filling narrow channels,
\item merging two shapes close to each other.
}
}
\examples{
# loading of toy-patient objects (decrease dxyz for better result)
step <- 4
patient <- toy.load.patient (modality = "mr", roi.name = "", 
                             dxyz = rep (step, 3))
MR <- patient$mr[[1]]

# generation of a binary volume
b <- bin.from.vol(MR, min = 15,max = 30)

b.dilation <- bin.dilation (b, radius = step)
display.plane (bottom = MR, top = b, main = "Before dilation", 
               view.coord = -20, interpolate = FALSE)
display.plane (bottom = MR, top = b.dilation, main = "After dilation", 
               view.coord = -20,interpolate = FALSE)
}
\seealso{
\link[espadon]{bin.erosion}, \link[espadon]{bin.opening},
\link[espadon]{bin.closing}, \link[espadon]{add.margin}, \link[espadon]{nesting.cube}.
}
