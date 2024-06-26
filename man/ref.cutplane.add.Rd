% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ref_cutplane_add.R
\name{ref.cutplane.add}
\alias{ref.cutplane.add}
\title{Adding volume's cutting planes frame of reference in T.MAT}
\usage{
ref.cutplane.add(
  vol,
  origin = c(0, 0, 0),
  ref.cutplane = paste0(vol$ref.pseudo, "m"),
  T.MAT = NULL
)
}
\arguments{
\item{vol}{"volume" class object.}

\item{origin}{Vector of the x, y, z coordinates of the origin point of
the cutting planes frame of reference.}

\item{ref.cutplane}{Name of the volume's cutting planes frame of reference.
By default \code{ref.cutplane = paste0 (vol$ref.pseudo,".m")}.}

\item{T.MAT}{"t.mat" class object created by
\link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm}
or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, then only the link
between \code{vol$ref.pseudo} and \code{ref.name} is established.}
}
\value{
Returns a "t.mat" class object, which contains the transfer matrices
from or to volume's cutting planes frame of reference. If the \code{T.MAT} is \code{NULL},
then the returned object will contain only 4 matrices: "src.ref<-src.ref",
"src.ref<-ref.cutplane", "ref.cutplane<-ref.cutplane", "ref.cutplane<-src.ref".
}
\description{
The \code{ref.cutplane.add} function adds in T.MAT the transfer
matrices from or to volume's cutting planes frame of reference.
}
\examples{
# loading of toy-patient objects
patient <- toy.load.patient (modality = "mr", roi.name = "", dxyz = c (4, 4, 4))
MR <- patient$mr[[1]]
MR$xyz.from.ijk

# creation of t.mat, containing the transfer matrix to the frame of reference 
# of the MR cutting planes
t.mat <- ref.cutplane.add (MR)

# Change of frame of reference
MR.m <- vol.in.new.ref (MR, paste0 (MR$ref.pseudo, "m"), t.mat)

MR.m$xyz.from.ijk
}
\seealso{
\link[espadon]{ref.add}, \link[espadon]{ref.srctodest.add},
\link[espadon]{ref.remove}.
}
