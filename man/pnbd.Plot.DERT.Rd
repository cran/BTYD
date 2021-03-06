% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pnbd.R
\name{pnbd.Plot.DERT}
\alias{pnbd.Plot.DERT}
\title{Pareto/NBD Plot Discounted Expected Residual Transactions}
\usage{
pnbd.Plot.DERT(params, x, t.x, T.cal, d, hardie = TRUE, type = "persp")
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

\item{d}{the discount rate to be used. Make sure that it matches up with your
chosen time period (do not use an annual rate for monthly data, for
example).}

\item{hardie}{if TRUE, use \code{\link{h2f1}} instead of
\code{\link[hypergeo]{hypergeo}}.}

\item{type}{must be either "persp" (perspective - 3 dimensional) or
"contour". Determines the type of plot produced by this function.}
}
\value{
A matrix with discounted expected residual transaction values for
every combination of calibration period frequency \code{x} and calibration
period recency \code{t.x}.
}
\description{
Plots discounted expected residual transactions for different combinations of
calibration period frequency and recency.
}
\details{
The length of the calibration period \code{T.cal} must be a single value, not
a vector.
}
\examples{
# The RFM and CLV paper uses all 78 weeks of the cdnow data to
# estimate parameters. These parameters can be estimated as follows:
# elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD2"),2,3)
# elog[, 'date'] <- as.Date(elog[, 'date'], format = '\%Y\%m\%d')
# cal.cbs <- dc.ElogToCbsCbt(elog)$cal$cbs
# pnbd.EstimateParameters(cal.cbs, hardie = TRUE)

# (The final function was run several times with its own output as
# input for starting parameters, to ensure that the result converged).

params <- c(0.5629966, 12.5590370, 0.4081095, 10.5148048)

# 15\% compounded annually has been converted to 0.0027 compounded continously,
# as we are dealing with weekly data and not annual data.
d <- 0.0027

pnbd.Plot.DERT(params = params, 
               x = 0:14, 
               t.x = 0:77, 
               T.cal = 77.86, 
               d = d, 
               hardie = TRUE, 
               type = "persp")
pnbd.Plot.DERT(params = params, 
               x = 0:14, 
               t.x = 0:77, 
               T.cal = 77.86, 
               d = d, 
               hardie = TRUE, 
               type="contour")
}
\references{
Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. "RFM and CLV:
Using Iso-Value Curves for Customer Base Analysis." Journal of Marketing
Research Vol.42, pp.415-430. November. 2005.
\url{http://www.brucehardie.com/papers.html}

Note that this paper refers to what this package is calling
discounted expected residual transactions (DERT) simply as discounted
expected transactions (DET).
}
