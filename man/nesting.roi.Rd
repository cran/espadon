% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nesting_roi.R
\name{nesting.roi}
\alias{nesting.roi}
\title{Restrict volume to RoI}
\usage{
nesting.roi(
  obj,
  struct,
  roi.name = NULL,
  roi.sname = NULL,
  roi.idx = NULL,
  xyz.margin = c(0, 0, 0),
  obj.restrict = FALSE,
  T.MAT = NULL,
  alias = "",
  description = NULL,
  ...
)
}
\arguments{
\item{obj}{object of class volume or mesh.}

\item{struct}{"struct" class object.}

\item{roi.name}{Vector of exact names of the RoI in the \code{struct} object.
By default \code{roi.name = NULL}. See Details.}

\item{roi.sname}{Names or parts of names of the RoI in the \code{struct} object.
By default \code{roi.sname = NULL}. See Details.}

\item{roi.idx}{Index of the RoI that belong to the \code{struct} object.
By default \code{roi.idx = NULL}. See Details.}

\item{xyz.margin}{Vector of length 3, representing the distances in mm to be added
to the x, y and z directions of the rectangular parallelepiped circumscribed
to the chosen RoI, in the cutting planes frame of reference. By default
xyz.margin = c (0, 0, 0).}

\item{obj.restrict}{Boolean. Used if \code{obj} is of class"volume". If
\code{obj.restrict = TRUE}, the rectangular parallelepiped circumscribed to
the selected voxels, enlarged by xyz.margin cannot exceed the initial volume.}

\item{T.MAT}{"t.mat" class object, created by
\link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm}
or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, \code{struct$ref.pseudo}
must be equal to \code{obj$ref.pseudo}.}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the the created object.
If \code{description = NULL}, it will be that of the \code{obj}, plus
"restricted to" the selected RoI.}

\item{...}{Additional arguments such as \code{vol} (depracated), replaced by \code{obj}.}
}
\value{
Returns a "volume" class object, in which 3D volume is limited to the
rectangular parallelepiped circumscribed to the chosen RoI, increased by the
requested margins.
}
\description{
The \code{nesting.roi} function restricts a "volume" class
object to the rectangular parallelepiped circumscribed to the chosen RoI.
}
\details{
If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all set
to NULL, all RoI are selected.
}
\examples{
# loading of toy-patient objects (decrease dxyz for  better result)
step <- 4
patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
                             roi.name = "brain", dxyz = rep (step, 3))
CT <- patient$ct[[1]]
S <- patient$rtstruct[[1]]

CT.brain <- nesting.roi (CT, S, roi.sname = "brain")
CT.brain.with.margin <- nesting.roi (CT, S, roi.sname = "brain",
                                         xyz.margin = c (10,10,10))

# display at the center of gravity of the cerebellum Gz
Gz <- S$roi.info [grep("^brain",S$roi.info$roi.pseudo),]$Gz
display.plane (bottom = CT.brain, view.coord = Gz,
               struct = S, bottom.col = pal.RVV (1000),
               bottom.breaks = seq (-1000, 1000, length.out = 1001),
               bg = "#00ff00",  interpolate  = FALSE, legend.shift = -20)
display.plane (bottom = CT.brain.with.margin,view.coord = Gz,
               struct = S,  bottom.col = pal.RVV (1000),
               bottom.breaks = seq(-1000, 1000, length.out = 1001),
               bg = "#00ff00", interpolate  = FALSE, legend.shift = -20)
}
\seealso{
\link[espadon]{add.margin}, \link[espadon]{nesting.cube} and
\link[espadon]{nesting.bin}.
}
