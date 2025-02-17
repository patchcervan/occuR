% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occu_tmb.R
\name{get_predicted_values}
\alias{get_predicted_values}
\title{Get predicted values for psi and p from parameters and matrices}
\usage{
get_predicted_values(fix, ran, mats)
}
\arguments{
\item{fix}{fixed parameter vector}

\item{ran}{random effects parameter vector}

\item{mats}{matrices from make_matrices}
}
\value{
list of two vectors: predicted values for psi and p
}
\description{
Get predicted values for psi and p from parameters and matrices
}
