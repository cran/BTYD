\name{cdnowSummary}
\docType{data}
\alias{cdnowSummary}
\title{CDNOW repeat transaction data summary}
\description{

  Data representing the purchasing behavior of 2357 CDNOW customers
  between January 1997 and June 1998, summarized as a customer-by-time
  matrix and a vector of cumulative weekly transactions.
  
}
\usage{data(cdnowSummary)}
\format{
 \describe{
 \item{$cbs}{A customer-by-time matrix with four columns: frequency
   (\dQuote{x}), recency (\dQuote{t.x}), length of observation in the calibration
   period (\dQuote{T.cal}), and number of transactions in the holdout period
   (\dQuote{x.star}). Each row represents a customer.}
 \item{$cu.tracking}{A vector containing cumulative transactions for
   every week in both the calibration and estimating periods (78 weeks
   total). This vector contains the sum of transactions across all
   customers.}
 \item{$est.params}{A vector containing estimated values for the four
   Pareto/NBD parameters: r, alpha, s, and beta, in that order. This
   estimation was made using \code{\link{pnbd.EstimateParameters}}, and
   is included here to avoid having to run the relatively time-consuming
   parameter estimation function in examples.}
  \item{$m.x}{A vector containing the average value of each customer's
    repeat transactions. Used in examples for spend functions.}
 }
}

\source{
  The data was put together using
  data conversion functions included in this package. The original event
  log is included (see \code{\link{cdnowElog}}).
  }

\details{

  The customers in this data represent 1/10th of the cohort of customers who
  made their first transactions with CDNOW in the first quarter of
  1997. CDNOW was an online retailer, selling music and related products
  on the web since 1994.
  
 }
\keyword{datasets}
