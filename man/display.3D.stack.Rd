% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display_3D_stack.R
\name{display.3D.stack}
\alias{display.3D.stack}
\title{Display in 3D the selected planes of an \pkg{espadon} class volume}
\usage{
display.3D.stack(
  vol,
  k.idx = unique(vol$k.idx[seq(1, vol$n.ijk[3], length.out = 10)]),
  display.ref = vol$ref.pseudo,
  T.MAT = NULL,
  col = grey.colors(10, start = 0, end = 1, alpha = c(rep(0, 1), rep(1, 9))),
  breaks = NULL,
  cube = TRUE,
  border = TRUE,
  ktext = TRUE,
  line.col = "#379DA2",
  line.lwd = 1,
  cex = 1
)
}
\arguments{
\item{vol}{"volume" class object to display.}

\item{k.idx}{vector of cutting plane numbers to be displayed, to be chosen in
\code{vol$k.idx}. By default \code{k.idx} is a vector of 10 uniformly
distributed cutting planes in the volume.}

\item{display.ref}{Character string. Pseudonym of the frame of reference used
for display.}

\item{T.MAT}{"t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm}
or \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{vol} must
be displayed in \code{display.ref = vol$ref.pseudo}.}

\item{col}{Vector, representing the color palette of the image. Transparent
colors are not represented.}

\item{breaks}{One of :
\itemize{
\item \code{NULL} : The minimum and the maximum value of the \code{vol} define
the range.
\item Vector giving the breakpoints of each color.
}}

\item{cube}{Boolean. If \code{TRUE} (default), the "volume" edges are displayed.}

\item{border}{Boolean. If \code{TRUE} (default), the borders of the planes defined
in \code{k.idx} are displayed.}

\item{ktext}{Boolean. If \code{TRUE} (default), the  selected cutting plane numbers
are displayed.}

\item{line.col}{Color of cube, planes and texts displayed.}

\item{line.lwd}{Line width of the border and cube, by default at 1.}

\item{cex}{Numeric character expansion factor of displayed plan numbers.}
}
\value{
Returns a display of the \code{k.idx} cutting planes of \code{vol},
in the current \pkg{RGL} window if it exists, in a new window otherwise. The
colors of the palettes are managed by \code{col} and \code{breaks}.
}
\description{
The \code{display.3D.stack} function displays in 3D the requested
cutting planes of a "volume" class object.
}
\examples{
# loading of toy-patient objects (decrease dxyz for better result)
step <- 4
patient <- toy.load.patient (modality = "ct", dxyz = rep (step, 3))

# display o 3 planes
library (rgl)
open3d()
display.3D.stack (patient$ct[[1]],
                  col = pal.RVV (200, alpha = c(rep(0,90), rep (1, 110))))
}
