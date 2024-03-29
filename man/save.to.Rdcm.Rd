% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_to_Rdcm.R
\name{save.to.Rdcm}
\alias{save.to.Rdcm}
\title{Save a \pkg{espadon} object in a pre-formatted *.Rdcm file}
\usage{
save.to.Rdcm(obj, object.name = obj$object.alias, dirname = obj$file.dirname)
}
\arguments{
\item{obj}{\pkg{espadon} object of class \code{"volume"}, \code{"struct"}, \code{"mesh"},
\code{"histo"}, \code{"dvh"}, \code{"histo2D"}.}

\item{object.name}{Character string, representing the name of the object,
default to \code{obj$object.alias}.}

\item{dirname}{Directory where new files from \code{obj} will be saved.}
}
\value{
Returns \code{TRUE}, if \code{paste0(object.name,".Rdcm")} exists in
\code{dirname}.

Returns \code{FALSE}, if \code{object.name} is not a valid file name,
or if the file that is created would replace a *.Rdcm file created by
\link[espadon]{dicom.to.Rdcm.converter}.
}
\description{
The function \code{save.to.Rdcm} allows you to save an object
created by \pkg{espadon} in a pre-formatted *.Rdcm file. This object will also be
accessible by the \code{load.patient.from.Rdcm} function.
}
\note{
\code{save.to.Rdcm} can not replace an *.Rdcm file created by
\link[espadon]{dicom.to.Rdcm.converter}.
}
\examples{
# loading of toy-patient objects (decrease dxyz for better result)
step <- 5
patient <- toy.load.patient (modality = c("ct", "mr"), roi.name = "", 
                             dxyz = rep (step, 3))
                             
# Save T.MAT to a temporary file pat.dir
pat.dir <- file.path (tempdir(), "PM_Rdcm") 
dir.create (pat.dir, recursive = TRUE) 
save.to.Rdcm (patient$ct[[1]], dirname = pat.dir)
save.to.Rdcm (patient$mr[[1]], dirname = pat.dir)
list.files(pat.dir)

# Cleaning  temporary directory
unlink (pat.dir, recursive = TRUE)
}
