% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rdcm_inventory.R
\name{Rdcm.inventory}
\alias{Rdcm.inventory}
\title{Inventory of \pkg{espadon} objects from Rdcm files}
\usage{
Rdcm.inventory(dirname, upgrade.to.latest.version = FALSE)
}
\arguments{
\item{dirname}{Character string, representing the full name of patient
directory, including Rdcm files.}

\item{upgrade.to.latest.version}{Boolean. If \code{TRUE}, the function attempts
to upgrade to the latest version, parsing the DICOM data. It may take longer
to load the data. Consider using the \link[espadon]{Rdcm.upgrade} function.}
}
\value{
Returns a dataframe, providing information of DICOM objects.
}
\description{
The \code{Rdcm.inventory} function creates, from Rdcm files in a
patient's directory, a dataframe describing objects.
}
\examples{
# First, save toy patient objects to a temporary file pat.dir for testing.
pat.dir <- file.path (tempdir(), "PM_Rdcm") 
dir.create (pat.dir, recursive = TRUE) 
patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
                             dxyz = c (4, 4, 4))
save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
save.T.MAT (patient$T.MAT, dirname = pat.dir)


Rdcm.inventory (pat.dir)

# Cleaning  temporary directory
unlink (pat.dir, recursive = TRUE)
}
