% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/struct_in_new_ref.R
\name{struct.in.new.ref}
\alias{struct.in.new.ref}
\title{Change of frame of reference of a "struct" class object.}
\usage{
struct.in.new.ref(struct, new.ref.pseudo, T.MAT, alias = "")
}
\arguments{
\item{struct}{"struct" class object.}

\item{new.ref.pseudo}{pseudonym of the frame of reference in which the struct
should be located. This \code{new.ref.pseudo} must exist in the \code{T.MAT} list.}

\item{T.MAT}{"t.mat" class object, created by
\link[espadon]{load.patient.from.dicom}, \link[espadon]{load.patient.from.Rdcm},
\link[espadon]{load.T.MAT} or \link[espadon]{ref.add}.}

\item{alias}{Character string, \code{$alias} of the created object.}
}
\value{
Returns "struct" class object in the new frame of reference
\code{new.ref.pseudo}.
}
\description{
The \code{struct.in.new.ref} function allows you to change the
frame of reference of a struct.
}
\examples{
# loading of toy-patient objects (decrease dxyz for better result)
step <- 4
patient <- toy.load.patient (modality = c("mr", "rtstruct"), roi.name = "", 
                             dxyz = rep (step, 3))
S <- patient$rtstruct[[1]]
S.in.new.ref <- struct.in.new.ref (S, patient$mr[[1]]$ref.pseudo, patient$T.MAT)
}
\seealso{
\link[espadon]{vol.in.new.ref}
}
