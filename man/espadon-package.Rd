% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/espadon-package.R
\docType{package}
\name{espadon-package}
\alias{espadon-package}
\alias{espadon}
\title{espadon: Easy Study of Patient DICOM Data in Oncology}
\description{
Exploitation, processing and 2D-3D visualization of DICOM-RT files (structures, dosimetry, imagery) for medical physics and clinical research, in a patient-oriented perspective.

\pkg{espadon} works in a native way (user friendly):
\itemize{
\item on images (CT, MR, PT) and fully manages changes of referential (REG)
\item on dosimetry (rt-dose)
\item on structures (rt-struct)
}

It is also able to use any DICOM format file, as long as the user knows where to look for the information.
In addition to the simplified use of the above-mentioned formats, espadon contains many functions that allow the user to rework documents to produce new features that can be integrated into machine learning.

\strong{\pkg{espadon} integrates functionalities:}
\itemize{
\item file loading, patient-centered information fusion
\item handling of 3D images
\itemize{
\item changes of referential
\item resampling, filtering
}
\item Contour manipulation, allowing the creation of new contours
\itemize{
\item fusion, intersection, inversion, erosion, dilation, opening, closing
\item segmentation based on the image
}
\item 2D representation
\itemize{
\item with application of masks
\item and transport of structures in the new reference frames
}
\item ... and 3D
\itemize{
\item production of mesh for the representation
}
\item for measurements
\itemize{
\item calculation of histograms (1D, 2D), DVH (integrating Monte-Carlo to simulate organ movements)
\item Measurement of surfaces, radii of curvature (on the mesh) and volumes
}
\item measurement of standard dosimetry indicators
}

In the calculation flow, various objects are created and accessible to the developer for his own use.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://espadon.cnrs.fr}
}

}
\author{
\strong{Maintainer}: Cathy Fontbonne \email{contact.espadon@lpccaen.in2p3.fr} (LPC-CAEN, France)

Authors:
\itemize{
  \item Jean-Marc Fontbonne (LPC-CAEN, France)
}

Other contributors:
\itemize{
  \item Nathan Azemar (LPC-CAEN, France) [contributor]
}


Cathy & Jean-Marc Fontbonne
}
\keyword{internal}
