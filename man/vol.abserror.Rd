% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vol_abserror.R
\name{vol.abserror}
\alias{vol.abserror}
\title{Absolute error volume}
\usage{
vol.abserror(vol, vol.ref, T.MAT = NULL, alias = "", description = NULL)
}
\arguments{
\item{vol, vol.ref}{"volume" class objects.}

\item{T.MAT}{"t.mat" class object to link the reference frames of \code{vol}
and \code{vol.ref}. \code{T.MAT} can be created by
\link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If
\code{T.MAT = NULL}, \code{vol$ref.pseudo} must be equal to \code{vol.ref$ref.pseudo}.}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL} (default value), it will be set to "absolute error".}
}
\value{
Returns a "volume" class object (see \link[espadon]{espadon.class}
for class definitions), with the same grid as \code{vol.ref}, and representing
the absolute error between \code{vol} and \code{vol.ref}.
}
\description{
The function \code{vol.abserror} provides the absolute error volume between 2 volumes.
}
\examples{
# loading of toy-patient objects (decrease dxyz)
step <- 5
patient <- toy.load.patient (modality = c("ct", "sct","rtstruct"), 
                             roi.name = c("eye", "brain","gizzard"),
                             dxyz = rep (step, 3))

patient$ct[[1]]$description
patient$ct[[2]]$description
# Creation of the absolute error volume between ct and synthetic ct
vAE <- vol.abserror (patient$ct[[2]], patient$ct[[1]], T.MAT = patient$T.MAT)

# Display
palette_vAE <- colorRampPalette(c("#00005F", "#0000FF", "#00FFFF", "#00FF00",
                                  "#FFFF00", "#FF7F00", "#FF0000", "#7F0000",
                                  "#5F0000")) (100)
breaks_vAE <- seq(floor (vAE$min.pixel), ceiling (vAE$max.pixel), 
                  length.out = 101)
                  
layout (mat = matrix(c(rep(1,6),2,2), ncol=4))                                
plot (vAE, view.coord =61, view.type = "trans", 
      col = palette_vAE, breaks = breaks_vAE)
display.palette(palette_vAE, breaks = breaks_vAE, 
                cex.axis = 1.2, main = vAE$unit)
par(mfrow=c(1,1))
}
\seealso{
\link[espadon]{vol.error}
}
