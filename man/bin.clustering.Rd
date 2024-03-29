% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_clustering.R
\name{bin.clustering}
\alias{bin.clustering}
\title{Binary volume clustering
\loadmathjax}
\usage{
bin.clustering(vol, alias = "", description = NULL)
}
\arguments{
\item{vol}{"volume" class object, of \code{"binary"} modality}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL} (default value), it will be set to
\code{paste (vol$object.alias,"clustering")}}
}
\value{
Returns a "volume" class object (see \link[espadon]{espadon.class}
for class definitions), of \code{"cluster"} modality. This object contains the
\code{$cluster.info} field, detailing the label and volumes in \mjeqn{cm^{3}}{ascii}
of the different clusters. Note that the label "0" is used for the background.
}
\description{
The \code{bin.clustering} function groups and labels TRUE voxels
that have a 6-connectivity (i.e. sharing a common side).
}
\examples{
# loading of toy-patient objects (decrease dxyz for better result)
step <- 4
patient <- toy.load.patient (modality = "ct", 
                             dxyz = rep (step, 3))
CT <- patient$ct[[1]]

# generation of a binary volume
b <- bin.from.vol(CT, min = -80, max = 20)

# Display of the n = 3 largest volumes
n <- 3
cluster.b <- bin.clustering (b)

col <- c ("#00000000", rainbow (n))
breaks <- seq (-0.5, n + 0.5, length.out = n+2)
display.plane (CT, top = b, main = "Before clustering",
               view.coord = 50, top.col = col, top.breaks = breaks, 
               interpolate = FALSE)
display.plane (CT, top = cluster.b, main = "After clustering", 
               view.coord = 50, top.col = col, top.breaks = breaks, 
               interpolate = FALSE)
}
