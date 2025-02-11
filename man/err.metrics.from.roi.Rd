% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/err_metrics_from_roi.R
\name{err.metrics.from.roi}
\alias{err.metrics.from.roi}
\title{Error evaluation metrics for 2-volume comparisons}
\usage{
err.metrics.from.roi(
  obj,
  obj.ref,
  struct = NULL,
  roi.name = NULL,
  roi.sname = NULL,
  roi.idx = NULL,
  T.MAT = NULL
)
}
\arguments{
\item{obj}{"volume" class object to be compared.}

\item{obj.ref}{"volume" class reference object.}

\item{struct}{"struct" class object or NULL.}

\item{roi.name}{Vector of exact names of the RoI in the \code{struct} object.
By default \code{roi.name = NULL}. See Details.}

\item{roi.sname}{Vector of names or parts of names of the RoI in the \code{struct} object.
By default \code{roi.sname = NULL}. See Details.}

\item{roi.idx}{Vector of indices of the RoI that belong to the \code{struct}
object.
By default \code{roi.idx = NULL}. See Details.}

\item{T.MAT}{"t.mat" class object to link the reference frames of \code{obj}
\code{obj.ref} and \code{struct}. \code{T.MAT} can be created by
\link[espadon]{load.patient.from.Rdcm} or \link[espadon]{load.T.MAT}. If
\code{T.MAT = NULL}, \code{struct$ref.pseudo} must be equal to
\code{obj$ref.pseudo} and \code{obj.ref$ref.pseudo}.}
}
\value{
Returns, in the zones delimited by the requested RoI, the following metrics:
\itemize{
\item ME: Mean Error
\item MAE: Mean Absolute Error
\item MSE: Mean Squared Error
\item RMSE: Root Mean Squared Error
}
}
\description{
The \code{err.metrics.from.roi} calculates various metrics
(ME, MAE, MSE, RMSE) to compare 2 “volume” class objects in the zones delimited
by the requested RoI.
}
\details{
If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are
all set to \code{NULL}, all RoI whose volume is greater than 0 are selected.
}
\examples{
# loading of toy-patient objects (decrease dxyz)
step <- 5
patient <- toy.load.patient (modality = c("ct", "sct","rtstruct"), 
                             roi.name = c("eye", "brain","gizzard"),
                             dxyz = rep (step, 3))

patient$ct[[1]]$description
patient$ct[[2]]$description
# Calculation of eye zone and ptv metrics
err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
                     struct = patient$rtstruct[[1]], roi.sname = c("eye","ptv"),
                     T.MAT= patient$T.MAT)
# Calculation of metrics on delimited zones on all RoIs
err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
                     struct = patient$rtstruct[[1]],
                     T.MAT= patient$T.MAT)
                     
# Calculation on all volume
err.metrics.from.roi(obj = patient$ct[[2]], obj.ref = patient$ct[[1]], 
                     T.MAT= patient$T.MAT)
                     
}
