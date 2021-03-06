% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_spectraComparison.R
\name{cospec}
\alias{cospec}
\title{Measure the cosine distance (spectra angle) of two MS2 spectra}
\usage{
cospec(measured, standard, ppmtol = 10, alpha = 1)
}
\arguments{
\item{measured}{the measured spectra}

\item{standard}{the standard spectra compare against}

\item{ppmtol}{ppm tolerance}

\item{alpha}{a numerical value. the power applied to the peak intensity, usually between 0-1 to give less weight 
on super high intensity peak (therefore, more weights on the number of matched peaks)}
}
\description{
Measure the cosine distance (spectra angle) of two MS2 spectra
}
