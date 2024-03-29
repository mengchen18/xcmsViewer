% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_consensusSpectrumLite.R
\name{consensusSpectrumLite}
\alias{consensusSpectrumLite}
\title{calculate consensus spectrum}
\usage{
consensusSpectrumLite(x, ppmtol = 20, minProp = 0.3, df = TRUE)
}
\arguments{
\item{x}{a data.frame with at least three columns named as "mz", 'intensity' and "ID". 
ID tells which peak is from which spectra.}

\item{ppmtol}{ppm tolerance}

\item{minProp}{minimum proportion that a peak should be included}

\item{df}{whether to return a data.frame, if false, paste into charecter vector}
}
\description{
calculate consensus spectrum
}
