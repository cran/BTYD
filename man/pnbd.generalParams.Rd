% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pnbd.R
\name{pnbd.generalParams}
\alias{pnbd.generalParams}
\title{Define general parameters}
\usage{
pnbd.generalParams(params, x, t.x, T.cal, func, hardie = TRUE)
}
\arguments{
\item{params}{Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
that order. r and alpha are unobserved parameters for the NBD transaction
process. s and beta are unobserved parameters for the Pareto (exponential
gamma) dropout process.}

\item{x}{number of repeat transactions in the calibration period T.cal, or a
vector of transaction frequencies.}

\item{t.x}{time of most recent repeat transaction, or a vector of recencies.}

\item{T.cal}{length of calibration period, or a vector of calibration period
lengths.}

\item{func}{name of the function calling dc.InputCheck; either \code{pnbd.LL}
or \code{pnbd.PAlive}.}

\item{hardie}{if TRUE, use \code{\link{h2f1}} instead of
\code{\link[hypergeo]{hypergeo}}.}
}
\value{
A vector of log likelihood values if \code{func} is \code{pnbd.LL},
or a vector of probabilities that a customer is still alive if \code{func}
is \code{pnbd.PAlive}.
}
\description{
This is to ensure consistency across all functions that require the
likelihood function, or the log of it, and to make sure that the same
implementation of the hypergeometric function is used everywhere for building
\code{A0}.
}
\details{
This function is only ever called by either \code{\link{pnbd.LL}} or
\code{\link{pnbd.PAlive}} so it returns directly the output that is expected
from those calling functions: either the log likelihood for a set of
customers, or the probability that a set of customers with characteristics
given by \code{x}, \code{t.x} and \code{T.cal}, having estimated a set of
\code{params}, is still alive. Either set of customers can be of size 1.
}
\seealso{
\code{\link{pnbd.LL}}

\code{\link{pnbd.PAlive}}

\code{\link{pnbd.DERT}}
}
