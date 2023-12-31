% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{ImputeData}
\alias{ImputeData}
\alias{ImputeData.AssayGrowth}
\alias{ImputeData.Rosace}
\title{Impute raw count data and get imputed count}
\usage{
ImputeData(object, impute.method, ...)

\method{ImputeData}{AssayGrowth}(object, impute.method, ...)

\method{ImputeData}{Rosace}(object, impute.method, name, key, ...)
}
\arguments{
\item{object}{An object: matrix (AssayGrowth), AssayGrowth, or Rosace}

\item{impute.method}{The imputation method to use: "knn" or "zero"}

\item{...}{Additional arguments to be passed to the imputation method.}

\item{name}{character string of AssayGrowth name}

\item{key}{key of AssayGrowth}
}
\value{
An object with imputed count
}
\description{
Impute the NA data either by K-Nearest Neighbor method or fill the NA with 0.
}
