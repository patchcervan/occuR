% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occu_tmb.R
\name{check_inp}
\alias{check_inp}
\title{Simple input checks for fit_occu}
\usage{
check_inp(forms, visit_data, site_data, start, print)
}
\arguments{
\item{forms}{list of formulae for psi (occupancy) and p (detection)}

\item{visit_data}{data.table with row for each visit with a "obs" column for detection record}

\item{site_data}{data.table with row for each site x occasion}

\item{start}{named list of starting values for beta_psi and beta_p - fixed effects parameters}

\item{print}{if TRUE then possibly useful info is printed out}
}
\value{
stops execution if input is invalid
}
\description{
Simple input checks for fit_occu
}
