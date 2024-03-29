% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ref_add.R
\name{ref.add}
\alias{ref.add}
\title{Adding a frame of reference in T.MAT}
\usage{
ref.add(
  src.ref,
  orientation = c(1, 0, 0, 0, 1, 0),
  origin = c(0, 0, 0),
  new.ref.pseudo = "newref",
  T.MAT = NULL
)
}
\arguments{
\item{src.ref}{Character string, pseudonym of the frame of reference in which
the \code{orientation} vector and the origin point \code{origin} are defined.}

\item{orientation}{Vector of 6 or 9 elements, composed of the coordinates of the 2 orthonormal vectors (i, j),
or of the 3 orthonormal vectors (i, j, k) of the new coordinate system,
in the \code{src.ref} frame of reference.}

\item{origin}{Vector of the x, y, z coordinates of the origin point of the
new frame of reference in the \code{src.ref} frame of reference. Default to
c (0, 0, 0).}

\item{new.ref.pseudo}{Character string, pseudonym of the new frame of
reference to add.}

\item{T.MAT}{"t.mat" class object created by
\link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm}
or \link[espadon]{load.T.MAT}. If \code{T.MAT = NULL}, then only the link
between \code{src.ref} and \code{new.ref.pseudo} is computed.}
}
\value{
Returns a "t.mat" class object, which contains the transfer
matrices from or to \code{new.ref.pseudo} pseudonym of the new frame
of reference. If the \code{T.MAT} is \code{NULL}, then the returned object
will contain only 4 matrices: "src.ref<-src.ref", "src.ref<-new.ref.pseudo",
"new.ref.pseudo<- new.ref.pseudo", "new.ref.pseudo<-src.ref".

Returns a \code{NULL} if  \code{orientation} is not well defined.
}
\description{
The \code{ref.add} function adds the transfer matrices from or
to a new frame of reference defined from 2 unit vectors and an origin point.
}
\examples{
# Adding of the reference frame "ref1_60", which is a 60 degree rotation of 
# reference frame "ref1".
orientation <- c (cos (pi / 3), sin (pi / 3), 0, 
                  -sin (pi / 3), cos (pi / 3), 0)

local.Tmat <- ref.add (src.ref = "ref1", orientation = orientation,
                  new.ref.pseudo = "ref1_60")

str(local.Tmat)
}
\seealso{
\link[espadon]{ref.cutplane.add}, \link[espadon]{ref.remove},
\link[espadon]{ref.srctodest.add}.
}
