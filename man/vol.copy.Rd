% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vol_copy.R
\name{vol.copy}
\alias{vol.copy}
\title{Creating a volume from another one}
\usage{
vol.copy(vol, alias = "", modality = NULL, description = NULL, number = NULL)
}
\arguments{
\item{vol}{"volume" class object, template of the created object.}

\item{alias}{Character string, \code{$object.alias} of the created object.}

\item{modality}{Character string, modality of the created volume. If
\code{modality = NULL}, then the created object will have the modality of
\code{vol}.}

\item{description}{Character string, description of the returned object. If
\code{descritption = NULL}, then the created object will have the description
of \code{vol}.}

\item{number}{number of the returned volume. If \code{number = NULL}, then
the returned object will have the number of \code{vol}.}
}
\value{
Returns a "volume" class object (see \link[espadon]{espadon.class}
for class definitions), with the same grid as \code{vol}, in which \code{$vol3D.data}
is initialized to \code{NA}.
}
\description{
The \code{vol.copy} function creates a "volume" class object,
with the same grid as the \code{vol} volume object.
}
\examples{
# loading of toy-patient objects (decrease dxyz for  better result)
patient <- toy.load.patient (modality = "ct", roi.name = "",dxyz = c (4, 4, 4))
CT <- patient$ct[[1]]

# creating a volume
vol.from.CT <- vol.copy (CT, alias = "ct reference")
str (vol.from.CT)
}
