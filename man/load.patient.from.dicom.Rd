% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_patient_from_dicom.R
\name{load.patient.from.dicom}
\alias{load.patient.from.dicom}
\title{Loading patient data from DICOM files}
\usage{
load.patient.from.dicom(
  dcm.files,
  data = FALSE,
  dvh = FALSE,
  ignore.duplicates = FALSE,
  tag.dictionary = dicom.tag.dictionary(),
  verbose = TRUE
)
}
\arguments{
\item{dcm.files}{String vector, representing the list of the full names of the
DICOM files of the same patient, or its directories.}

\item{data}{Boolean. If \code{data = TRUE}, the voxels value of the "volume"
class objects, or the coordinates of the RoI (region of interest)
of the \code{struct} class objects, are loaded into memory.}

\item{dvh}{Boolean. if \code{dvh = TRUE} and if they exist, patient DVH are
loaded, for convenience. They are not used as is in \pkg{espadon} package.}

\item{ignore.duplicates}{Boolean. If \code{TRUE}, the function ignores duplicated objects.}

\item{tag.dictionary}{Dataframe, by default equal to
\link[espadon]{dicom.tag.dictionary},
whose structure it must keep. This dataframe is used to parse DICOM files.}

\item{verbose}{Boolean. If \code{TRUE}, a progress bar indicates loading progress.}
}
\value{
Returns an \pkg{espadon} object of class "patient", describing the
information from \code{dcm.files}. See \link[espadon]{espadon.class} for a
description of the "patient" class.
}
\description{
The \code{load.patient.from.dicom} function is used to load or
pre-load in memory all patient objects from DICOM files.
}
\examples{
# First, save toy.dicom.raw () raw data to a temporary file pat.dir for testing.
pat.dir <- file.path (tempdir(), "toy_dcm") 
dir.create (pat.dir, recursive = TRUE) 
dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = pat.dir,
                          fileext = ".dcm")
zz <- file (dcm.filename, "wb")
writeBin (toy.dicom.raw (), zz, size = 1)
close (zz)

# loading patient. Here the toy patient ha only a unique rt-plan object
patient <- load.patient.from.dicom (pat.dir, data = FALSE)
str (patient, max = 2)
# description of object
patient$description
# transfer matrices :
patient$T.MAT
# rt-plan object
str (patient$rtplan[[1]])
# Cleaning  temporary directory
unlink (pat.dir, recursive = TRUE)
}
\seealso{
\link[espadon]{dicom.to.Rdcm.converter}, \link[espadon]{load.patient.from.Rdcm},
\link[espadon]{load.obj.data},  \link[espadon]{load.obj.from.dicom},
\link[espadon]{load.obj.from.Rdcm} and \link[espadon]{load.T.MAT}.
}
