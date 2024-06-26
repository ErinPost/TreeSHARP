% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtm_normal_stan.R
\name{dtm_normal_stan}
\alias{dtm_normal_stan}
\title{Bayesian Dirichlet Tree Multinomial regression with Stan}
\usage{
dtm_normal_stan(X, Y, tree, sigma, ...)
}
\arguments{
\item{X}{Numeric matrix of input values.}

\item{Y}{Numeric matrix of output values.}

\item{tree}{An object of class \code{ape} specificying the DTM structure of Y.}

\item{sigma}{A positive value indicating the SD of the normal distribution}

\item{...}{Arguments passed to \code{rstan::sampling} (e.g. iter, chains).}
}
\value{
An object of class \code{stanfit} returned by \code{rstan::sampling}
}
\description{
Bayesian Dirichlet Tree Multinomial regression with Stan
}
