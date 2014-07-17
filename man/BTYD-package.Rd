\name{BTYD-package}
\alias{BTYD-package}
\alias{BTYD}
\docType{package}
\title{
Implementing Buy 'Til You Die Models
}
\description{
  
This package contains functions for data preparation, parameter estimation, scoring,
and plotting for the BG/BB, BG/NBD and Pareto/NBD models.

}
\details{
  \tabular{ll}{
    Package: \tab BTYD\cr
    Type: \tab Package\cr
    Version: \tab 2.3\cr
    Date: \tab 2014-07-17\cr
    License: \tab GPL-3\cr
    }

  This project was funded and sponsored by the Wharton Customer Analytics
  Initiative (wcai.wharton.upenn.edu).

  This package implements the BG/BB, BG/NBD and Pareto/NBD models, which
  capture/project customer purchase patterns in a typical
  non-contractual setting.

  While these models are developed on a customer-by-customer basis, they
  do not necessarily require data at such a granular level. The
  Pareto/NBD requires a \dQuote{customer-by-sufficient-statistic} matrix
  (CBS), which consists of each customer's frequency, recency (the time
  of their last transactions) and total time observed - but the timing
  of each and every transaction (other than the last) is not needed by
  the model. If, however, you do have the granular data in the form of
  an event log (which contains at least columns for customer
  identification and the time of each transaction, and potentially more
  columns such as transaction amount), this package provides functions
  to convert it to a CBS. You can use \code{\link{dc.ReadLines}} to get
  your event log from a comma-delimited file to an event log usable by
  this package; it is possible to use read.table or read.csv, but
  formatting will be required afterwards. You can then convert the event
  log directly to a CBS (for both the calibration and holdout periods)
  using \code{\link{dc.ElogToCbsCbt}}. As the name suggests, this
  function also produces a customer-by-time matrix (CBT). This matrix
  consists of a row for every customer and a column for every date, and
  is populated by a statistic of your choice (reach, frequency, or
  spend). It is not necessary for any of the models presented in this
  package, but is used as a building block to produce the CBS.
  
  The BG/NBD model requires all the same inputs as the Pareto/NBD model.

  The BG/BB model requires the same information as the Pareto/NBD model,
  but as it models discrete transaction opportunities, this information
  can be condensed into a recency-frequency matrix. A recency-frequency
  matrix contains a row for every recency/frequency combination in the
  given time period, and each row contains the number of customers with
  that recency/frequency combination. Since frequency will always be
  less than or equal to recency, this matrix will contain (n)(n-1)/2 + 1
  rows at most, with n as the number of transaction opportunities (of
  course, the maximum number of rows for pooled data - for customers
  with varying numbers of transaction opportunities - will be the sum of
  the above equation for each unique number of transaction
  opportunities). You can convert a CBS to recency-frequency matrices
  using \code{\link{dc.MakeRFmatrixCal}} and
  \code{\link{dc.MakeRFmatrixHoldout}}.

  If you want to test the data contained in the package, or have data
  formatted as a customer-by-sufficient-statistic or recency-frequency
  matrix, a good starting place would be
  \code{\link{pnbd.EstimateParameters}}, 
  \code{\link{bgnbd.EstimateParameters}}, or
  \code{\link{bgbb.EstimateParameters}}.

  Following that, \code{\link{pnbd.PlotFrequencyInCalibration}}, 
  \code{\link{bgnbd.PlotFrequencyInCalibration}} and
  \code{\link{bgbb.PlotFrequencyInCalibration}} will give a check that
  the model fits the data in-sample. Further plotting functions,
  comparing actual and expected results, are labelled
  \dQuote{pnbd.Plot...}, \dQuote{bgnbd.Plot...} and \dQuote{bgbb.Plot...}. 
  The building blocks of these functions are also provided: 
  \code{\link{pnbd.LL}}, \code{\link{bgnbd.LL}}
  \code{\link{bgbb.LL}}, \code{\link{pnbd.pmf}},
  \code{\link{bgnbd.pmf}}, \code{\link{bgbb.pmf}},
  \code{\link{pnbd.Expectation}}, \code{\link{bgnbd.Expectation}}, 
  \code{\link{bgbb.Expectation}}, 
  \code{\link{pnbd.ConditionalExpectedTransactions}}, 
  \code{\link{bgnbd.ConditionalExpectedTransactions}}, and
  \code{\link{bgbb.ConditionalExpectedTransactions}} may be of
  particular interest.

  This package uses the following conventions:

  The time period used to estimate the model parameters is called the
  \emph{calibration period}. Users may be accustomed to this being
  called the estimation period, or simply being referred to as
  \dQuote{in-sample}. Function parameter names generally follow this
  convention: for example, \dQuote{n.cal} is used to refer to the number
  of transaction opportunities in the calibration period.

  The time period used to validate model performance is called the
  \emph{holdout period}. Users may be accustomed to this being called
  the validation period, or simply being referred to as
  \dQuote{out-of-sample}. Function parameters relating to this time
  period are generally appended with \dQuote{.star}. For example, n.star
  is used to refer to the number of transaction opportunities in the
  holdout period.

  As described in the papers referenced below, the BG/BB, BG/NBD and 
  Pareto/NBD   models are generally concerned with repeat transactions, 
  not total transactions. This means that a customer's first transaction 
  in the   calibration period is usually not part of the data being 
  modelled - this is due to the fact that a new customer generally does 
  not show up \dQuote{on the comapany's radar} until after their first 
  purchase has taken place. This means that the modal number of repeat 
  purchases tends to be zero. If your data does not have a relatively large 
  number of customers with zero transactions, but does have a relatively large
  number of customers with one transaction, and the estimation functions
  are struggling, the problem is most likely that you are including
  customers' very first transactions. Some of the data-conversion
  functions have examples illustrating how to work with data that
  includes this very first transaction. Note that this does not apply to
  the holdout period; in the holdout period, we already know about the
  customer and take all of their previous transactions into account.

}
\author{
  Lukasz Dziurzynski, Daniel McCarthy, Edward Wadsworth

  Contributors: Peter Fader, Elea McDonnell Feit, Bruce Hardie, Arun Gopalakrishnan, Eric Schwartz, Yao Zhang, Elea McDonnell Feit

  Maintainer: Daniel McCarthy <danielmc@wharton.upenn.edu>
}
\references{

See \url{www.brucehardie.com} for papers, notes, and datasets relating
to applied probability models in marketing.

Fader, Peter S., and Bruce G.S. Hardie. \dQuote{A Note on Deriving the
Pareto/NBD Model and Related Expressions.}
November. 2005. Web. \url{http://www.brucehardie.com/notes/008/}

Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. \dQuote{RFM and CLV:
  Using Iso-Value Curves for Customer Base Analysis.} \emph{Journal of
  Marketing Research} Vol.42, pp.415-430. November. 2005. \url{http://www.brucehardie.com/papers.html}

Fader, Peter S., and Bruce G.S. Hardie. \dQuote{Deriving an Expression
  for P (X(t) = x) Under the Pareto/NBD Model.}
September. 2006. Web. \url{http://www.brucehardie.com/notes/012/}

Fader, Peter S., and Bruce G.S. Hardie. \dQuote{Creating an RFM summary
  using Excel.} December. 2008. Web. \url{http://www.brucehardie.com/notes/022/}

Fader, Peter S., Bruce G.S. Hardie, and Jen Shang. \dQuote{Customer-Base
  Analysis in a Discrete-Time Noncontractual Setting.} \emph{Marketing
  Science} 29(6),
  pp. 1086-1108. 2010. INFORMS. \url{http://www.brucehardie.com/papers/020/}

  Jerath, Kinshuk, Peter S. Fader, and Bruce
G.S. Hardie. \dQuote{Customer-Base Analysis on a 'Data Diet': Model
Inference Using Repeated Cross-Sectional Summary (RCSS) Data.}
June. 2011. Available at SSRN: \url{http://ssrn.com/abstract=1708562} or
\url{http://dx.doi.org/10.2139/ssrn.1708562}

Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. \dQuote{``Counting Your Customers'' the Easy Way:
An Alternative to the Pareto/NBD Model.} \emph{Marketing Science} Vol.24, pp.275-284. Spring. 2005. 
\url{http://www.brucehardie.com/papers.html}

Fader, Peter S., Hardie, Bruce G.S., and Lee, Ka Lok. \dQuote{Computing P(alive) Using
the BG/NBD Model.} December. 2008. Web. \url{http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf}

}

\keyword{ package }

