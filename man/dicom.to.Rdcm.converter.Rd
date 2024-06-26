% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dicom_to_Rdcm_converter.R
\name{dicom.to.Rdcm.converter}
\alias{dicom.to.Rdcm.converter}
\title{Conversion of DICOM object into files that can be interpreted by the \pkg{espadon}
package}
\usage{
dicom.to.Rdcm.converter(
  dcm.files,
  pat.dest.dir,
  update = TRUE,
  ignore.duplicates = FALSE,
  tag.dictionary = dicom.tag.dictionary(),
  verbose = TRUE
)
}
\arguments{
\item{dcm.files}{String vector, representing the list of the full names of
the DICOM files of the same patient, or its directory.}

\item{pat.dest.dir}{Character string representing the full name of the patient's
directory, which will contain files converted into Rdcm objects.}

\item{update}{Boolean. If set to \code{TRUE}, and if \code{pat.dest.dir}
contains previously converted files, these files are updated,even if they are
duplicated. They retain the same \pkg{espadon} reference frame assignment.}

\item{ignore.duplicates}{Boolean. If \code{TRUE}, the function ignores duplicated objects.}

\item{tag.dictionary}{Dataframe, by default equal to
\link[espadon]{dicom.tag.dictionary},
whose structure it must keep. This dataframe is used to parse DICOM files.}

\item{verbose}{Boolean. If \code{TRUE}, a progress bar indicates the progress
of the conversion.}
}
\value{
Returns the list of basenames of the created files.

Returns \code{NULL} if there are no DICOM files in \code{dcm.files}
}
\description{
The \code{dicom.to.Rdcm.converter} function creates, for each DICOM
object, a *.Rdcm file usefull for using \pkg{espadon} package. Each Rdcm file created is
referenced by the date of acquisition of the object (if it is not available,
its creation date), the patient's PIN, a reference number, an object number
in this reference system, and the object modality (mr, ct, rtstruct...).
}
\note{
For each DICOM object, \code{dicom.to.Rdcm.converter}
creates a *.Rdcm file whose basename is made up of the date of the acquisition
(or creation date if previous not found), the patient's PIN, the pseudonym of
the frame of reference ("ref1", "ref2"...), the number of the volume object
in the directory in this frame of reference ("do1", "do2"...), and the object
modality ("mr", "ct", "rtdose", "rtstruct"...).

For example: \code{BASE = "20160514_a008e9ac_ref2_do1_mr"}
}
\examples{
# First, save toy.dicom.raw () raw data to a temporary file for testing.
pat.src.dir <- file.path (tempdir(), "PM_dcm") 
dir.create (pat.src.dir, recursive = TRUE) 
dcm.filename <- tempfile (pattern = "PM_rtplan", tmpdir = pat.src.dir,
                          fileext = ".dcm")
zz <- file (dcm.filename, "wb")
writeBin (toy.dicom.raw (), zz, size = 1)
close (zz)

# Create a temporary destination directory where the *.Rdcm file will be saved
pat.dest.dir <- file.path (tempdir(), "PM_Rdcm")

dicom.to.Rdcm.converter (pat.src.dir, pat.dest.dir, update = TRUE)
# or
dicom.to.Rdcm.converter (dcm.filename, pat.dest.dir, update = TRUE)

list.files (pat.dest.dir)

# Cleaning  temporary directories
unlink (pat.src.dir, recursive = TRUE)
unlink (pat.dest.dir, recursive = TRUE)
}
