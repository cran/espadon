% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vol_repair.R
\name{vol.repair}
\alias{vol.repair}
\title{repairing missing planes of volumes}
\usage{
vol.repair(vol, alias = "", description = NULL)
}
\arguments{
\item{vol}{"volume" class object.}

\item{alias}{Character string, \code{$alias} of the created object.}

\item{description}{Character string, describing the created object. If
\code{description = NULL} (default value), it will be set to
\code{paste (vol$object.alias, "repair")}.}
}
\value{
Returns a "volume" class object (see \link[espadon]{espadon.class}
for class definitions), with no missing plane, if \code{vol} is to be repaired.
Returns \code{vol} otherwise.
}
\description{
The \code{vol.repair} function repairs missing planes in volumes.
}
\details{
Missing planes at download can generate errors or unpredictible results
in espadon processing. The \code{vol.repair} function detects such missing
planes and recreates their value by interpolation.
}
\examples{
step <- 4
patient <- toy.load.patient (modality = c("ct", "mr", "rtstruct", "rtdose"),
                             roi.name  = "",
                             dxyz = rep (step, 3), beam.nb = 3)
CT <- patient$ct[[1]]

# this function removes a plane in a volume in order to simulate
# a dicom transfer issue
remove.plane <- function (vol, k) {
  idx <- which (vol$k.idx == k)
  vol$n.ijk[3] <- vol$n.ijk[3] - 1
  vol$xyz0 <- vol$xyz0[-idx, ]
  vol$k.idx <- vol$k.idx[-idx]
  vol$missing.k.idx <- TRUE
  vol$vol3D.data <- vol$vol3D.data[, , -idx]
  return (vol)
}

# Creation of CT.damaged without the 29th slice.
CT.damaged<- remove.plane (CT, 29)
CT.fix <- vol.repair (CT.damaged)

# Display 
par (mfrow=c(3, 3))
for (k in 28:30) {
display.kplane (CT, k, main = paste("CT @ k =",k),interpolate = FALSE)
display.kplane (CT.damaged, k, main = "damaged CT",interpolate = FALSE)
display.kplane (CT.fix, k, main = "fixed CT", interpolate = FALSE)
}
}
