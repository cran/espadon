% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_deployment.R
\name{study.deployment}
\alias{study.deployment}
\title{Deployment of DICOM files from multiple patients}
\usage{
study.deployment(
  pats.dir,
  deploy.dir,
  design.matrix = matrix(TRUE, nrow = length(dir(pats.dir)), ncol = 1, dimnames =
    list(basename(dir(pats.dir)), "expert_1")),
  pid.prefix = "",
  white.list = c("instance", "reference"),
  black.list = c("frame of reference", "class"),
  tag.dictionary = dicom.tag.dictionary()
)
}
\arguments{
\item{pats.dir}{Name of the directory in which all patient directories are
stored, each containing the DICOM files to be duplicated.}

\item{deploy.dir}{Name of the directory where all patient files will be
duplicated.}

\item{design.matrix}{Boolean matrix. See Details.}

\item{pid.prefix}{string vector of length 1 or string vector of length
\code{ncol(design.matrix)}, representing the prefix added to the new
unique identifier of the deployed patient (tag (0010,0020)).}

\item{white.list}{Names vector, representing a part of the DICOM tag name
UI value representation, other than those defined by the DICOM
standard, which will be modified. By default, the UID name containing 'instance'
or 'reference' will be modified.}

\item{black.list}{Names vector, representing a part of the DICOM tag name
UI value representation, other than those defined by the DICOM
standard, which will not be modified. By default, the frame of reference UID
will not be modified.}

\item{tag.dictionary}{Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary},
whose structure it must keep. This dataframe is used to parse DICOM files.}
}
\value{
Creates the \code{deploy.dir} directory, containing the expert
directories defined by the \code{design.matrix} column names. Each expert
directory contains as many patient directories as defined by the
\code{design.matrix} row names. All patients will be independent of each other.
The new created patients have the pats.dir as name, and expert name as first
name, and an independent patient ID, with prefix \code{pid.prefix}.
}
\description{
The \code{study.deployment} function duplicates DICOM data from
multiple patients, so that it becomes data independent of the original data.
This function simplifies the analysis of multi-center or multi-expert studies
in dosimetry challenges, contouring consensus searches, etc.
}
\details{
The \code{design.matrix} argument defines how patients DICOM files will be deployed.
The names of the lines must match the names of the directories contained in \code{pats.dir}.
The names of the columns are for example the different experts or hospitals
who will study the patient files. These experts will only review the patients
files defined by \code{rownames(design.matrix)[design.matrix[,"expert"]]}.
}
\examples{
# First, save toy.dicom.raw () raw data to a temporary file/pats.dir/toy_PM 
# for testing.
toy_PM.dir <- file.path (tempdir(), "pats.dir","toy_PM") 
dir.create (toy_PM.dir, recursive = TRUE) 
dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = toy_PM.dir,
                          fileext = ".dcm")
zz <- file (dcm.filename, "wb")
writeBin (toy.dicom.raw (), zz, size = 1)
close (zz)

# function test:
pats.dir <- dirname (toy_PM.dir)
deploy.dir <- file.path (tempdir(), "deploy.dir") 
design.matrix <- matrix(TRUE, nrow = length (dir (pats.dir)), ncol=3,
                        dimnames = list (basename (dir (pats.dir)),
                                         c("Dr Quinn","Dr Who","Dr House")))
design.matrix
study.deployment (pats.dir, deploy.dir, design.matrix, 
                 pid.prefix = c("zz_", "yy_", "xx_"))

# check result
list.files(deploy.dir, recursive = TRUE)
load.patient.from.dicom(deploy.dir)$patient
# Cleaning  temporary directory
unlink (pats.dir, recursive = TRUE)
unlink (deploy.dir, recursive = TRUE)
}
