% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_volume_from_roi.R
\name{get.volume.from.roi}
\alias{get.volume.from.roi}
\title{Volume of a region of interest (RoI)
\loadmathjax}
\usage{
get.volume.from.roi(struct, roi.name = NULL, roi.sname = NULL, roi.idx = NULL)
}
\arguments{
\item{struct}{"struct" class object.}

\item{roi.name}{Vector of exact names of the RoI in the \code{struct} object.
By default \code{roi.name = NULL}. See Details.}

\item{roi.sname}{Vector of names or parts of names of the RoI in the \code{struct}
object. By default \code{roi.sname = NULL}. See Details.}

\item{roi.idx}{Vector of indices of the RoI that belong to the \code{struct}
object. By default \code{roi.idx = NULL}. See Details.}
}
\value{
Returns a vector of the volumes in \mjeqn{cm^3}{ascii} of the requested
RoI.
}
\description{
The \code{get.volume.from.roi} function extracts the volume
in\mjeqn{cm^3}{ascii} of one or more RoI, from the \code{$roi.info} of the
"struct" class object.
}
\details{
If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all set
to NULL, all RoI are selected.
}
\examples{
# loading of toy-patient objects
step <- 4
patient <- toy.load.patient (modality = c ("rtstruct"), 
                             dxyz = rep (step, 3))
S <- patient$rtstruct[[1]] 

# Volume extraction
vol <- get.volume.from.roi (S, roi.sname = "bra", roi.idx = c (1, 3))
names (vol)
vol
}
\seealso{
\link[espadon]{get.volume.from.bin}, \link[espadon]{select.names}
}
