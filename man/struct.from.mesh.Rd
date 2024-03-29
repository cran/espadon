% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/struct_from_mesh.R
\name{struct.from.mesh}
\alias{struct.from.mesh}
\title{Creation of struct class object from an espadon mesh}
\usage{
struct.from.mesh(
  mesh,
  z,
  thickness = NULL,
  roi.name = mesh$object.alias,
  roi.nb = 1,
  roi.color = "#ff0000",
  roi.type = "",
  alias = "",
  description = NULL,
  force.closed = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{mesh}{espadon mesh class object.}

\item{z}{z-coordinate vector where mesh contours are computed.}

\item{thickness}{struct thickness between 2 adjacent contours. If NULL (default),
it is deduced from \code{z}.}

\item{roi.name}{Character string, representing the name of created RoI.}

\item{roi.nb}{Positive integer, representing the number of created RoI.}

\item{roi.color}{Color of the created RoI, in hex code format ("#RRGGBB").}

\item{roi.type}{Type of RoI, from among "", "EXTERNAL", "PTV", "CTV", "GTV",
"TREATED_VOLUME", "IRRAD_VOLUME", "OAR", "BOLUS", "AVOIDANCE", "ORGAN", "MARKER",
"REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", "CAVITY", "BRACHY_CHANNEL",
"BRACHY_ACCESSORY", "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", "SUPPORT", "FIXATION",
"DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT".}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the the created object.}

\item{force.closed}{Boolean. Set to \code{TRUE} when the mesh represents the
surface of a closed volume.}

\item{verbose}{Boolean. If \code{TRUE} (default), a progress bar indicates
the state of calculation.}

\item{...}{Others parameters}
}
\value{
Returns a "struct" class object (see \link[espadon]{espadon.class}
for class definition), including the unique \code{roi.name} as region of interest.
}
\description{
The \code{struct.from.mesh} function creates a struct object with
a unique RoI, defined by the contours of a mesh.
}
\examples{
# Creation of an espadon mesh of a cube
M <- obj.create (class = "mesh")
M$mesh <- Rvcg::vcgIsotropicRemeshing (Rvcg::vcgBox(),0.5) 
M$nb.faces <- ncol (M$mesh$it)

S <- struct.from.mesh (M, z = seq(-1,1,0.5))
display.3D.contour(S)
}
