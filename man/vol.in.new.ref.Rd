% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vol_in_new_ref.R
\name{vol.in.new.ref}
\alias{vol.in.new.ref}
\title{Change of frame of reference of a volume}
\usage{
vol.in.new.ref(vol, new.ref.pseudo, T.MAT, alias = "", description = NULL)
}
\arguments{
\item{vol}{"volume" class object.}

\item{new.ref.pseudo}{pseudonym of the frame of reference in which the volume
should be located. This \code{new.ref.pseudo} must exist in the \code{T.MAT}
list.}

\item{T.MAT}{"t.mat" class object, created by
\link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
\link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL} (default value), it will be that of the \code{vol} volume.}
}
\value{
Returns "volume" class object in the new frame of reference
\code{new.ref.pseudo}.
}
\description{
The \code{vol.in.new.ref} function allows you to change the
frame of reference of a volume.
}
\examples{
# loading of toy-patient objects (decrease dxyz for  better result)
patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
                             dxyz = c (4, 4, 4))
CT <- patient$ct[[1]]
CT.in.new.ref <- vol.in.new.ref (CT, patient$mr[[1]]$ref.pseudo, patient$T.MAT)
}
\seealso{
\link[espadon]{struct.in.new.ref}
}
