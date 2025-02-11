% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/struct_create.R
\name{struct.create}
\alias{struct.create}
\title{Struct creating from contours list}
\usage{
struct.create(
  contours.list,
  roi.name,
  roi.nb = 1,
  roi.color = "#ff0000",
  roi.type = "",
  ref.pseudo = "ref1",
  frame.of.reference = "",
  alias = "",
  description = NULL
)
}
\arguments{
\item{contours.list}{list of data frames or 3-column matrices, representing the
xyz coordinates of polygon vertices.  Each polygon must have a unique
z-coordinate. All coordinates must be uniformly distributed according to a
unique inter-slice distance.}

\item{roi.name}{Character string, representing the name of created RoI.}

\item{roi.nb}{Positive integer, representing the number of created RoI.}

\item{roi.color}{Color of the created RoI, in hex code format ("#RRGGBB").}

\item{roi.type}{Type of RoI, from among "", "EXTERNAL", "PTV", "CTV", "GTV",
"TREATED_VOLUME", "IRRAD_VOLUME", "OAR", "BOLUS", "AVOIDANCE", "ORGAN", "MARKER",
"REGISTRATION", "ISOCENTER", "CONTRAST_AGENT", "CAVITY", "BRACHY_CHANNEL",
"BRACHY_ACCESSORY", "BRACHY_SRC_APP", "BRACHY_CHNL_SHLD", "SUPPORT", "FIXATION",
"DOSE_REGION","CONTROL" and "DOSE_MEASUREMENT".}

\item{ref.pseudo}{Character string, frame of reference pseudonym of the
created object.By defaukt equal to "ref1"}

\item{frame.of.reference}{Character string, frame of reference of the
created object.}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the the created object.}
}
\value{
Returns a "struct" class object (see \link[espadon]{espadon.class}
for class definition), including the unique \code{roi.name} as region of interest.
}
\description{
The \code{struct.create} function creates a struct object from a
list of polygons, representing the contours of a shape.
}
\examples{
contours.z <- -50:50
theta <- seq(0,2*pi, length.out = 100)
contours <- lapply(contours.z,function(z){
  if (z<(-25)) return(data.frame(x = (50 + z) * cos(theta),
                                 y = (50 + z) * sin(theta),
                                 z = z))
  if (z>25) return(data.frame(x = (50 - z) * cos(theta),
                              y = (50 - z) * sin(theta),
                              z = z))
  return(data.frame(x = 25 * cos(theta),
                    y = 25 * sin(theta),
                    z = z))
})

contours <- contours[!sapply(contours, is.null)]
S <- struct.create(contours, roi.name="myshape",
                   roi.nb = 1,
                   roi.color = "#ff0000",
                   roi.type = "",
                   ref.pseudo = "ref1", 
                  alias="", description = NULL)
display.3D.contour(S)
}
\seealso{
\link[espadon]{struct.from.mesh}.
}
