% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegressionModels.R
\name{linmod}
\alias{linmod}
\title{linmod}
\usage{
linmod(x, y, xlab = "X", ylab = "Y")
}
\arguments{
\item{x}{Numeric vector}

\item{y}{Numeric vector}

\item{xlab}{Character string. X variable label.}

\item{ylab}{Character string. Y variable label.}
}
\description{
Linear Regression between x and y. Returns data.frame to be displayed in UI
containing relationship (i.e linear), model (ie WTonCL), estimate and pval.
}
