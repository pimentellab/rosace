% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runROSACE.R
\name{MCMCDiagnostics}
\alias{MCMCDiagnostics}
\title{Extract MCMC diagnostics from CmdStanMCMC object}
\usage{
MCMCDiagnostics(fit, sampler = FALSE)
}
\arguments{
\item{fit}{CmdStanMCMC object}

\item{sampler}{whether to get diagnostics summary or
extract detailed sampler diagnostics}
}
\value{
list of diagnostics summary or data.frame of sampler diagnostics
}
\description{
Extract MCMC diagnostics from CmdStanMCMC object
}
