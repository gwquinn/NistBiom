% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/det.R
\name{get.FMR}
\alias{get.FMR}
\title{Get FMR at the specified threshold(s).}
\usage{
get.FMR(DET, ...)
}
\arguments{
\item{DET}{A data.table object or path to a DET file.}

\item{t}{A single decision threshold (or vector of thresholds).}
}
\value{
A vector of false match rates the same length as t.
}
\description{
Get FMR at the specified threshold(s).
}
