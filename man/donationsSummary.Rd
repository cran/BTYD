\name{donationsSummary}
\alias{donationsSummary}
\docType{data}
\title{Discrete donation data summary}
\description{ 

  This dataset contains a recency-frequency matrix capturing the
  discrete transaction behavior of 11,104 customers over 6 transaction
  opportunities, summarized as a recency-frequency matrix and a vector
  of annual transactions.
  
}
\usage{data(donationsSummary)}
\format{
  \describe{
    \item{$rf.matrix}{A matrix with 22 rows (for each possible
      recency-frequency combination in 6 calibration period transaction
      opportunities) and 4 columns: number of transactions during the calibration
      period (\dQuote{x}), recency in the calibration period
      (\dQuote{t.x}), number of transaction opportunities in the
      calibration period (\dQuote{n.cal}), and number of customers with
      this recency-frequency combination in the calibration period
      (\dQuote{custs}).}
    \item{$rf.matrix.holdout}{A matrix with 15 rows (for each possible
      recency-frequency combination in 5 holdout period transaction
      opportunities) and 4 columns: number of transactions during the
      holdout period (\dQuote{x.star}), recency in the holdout period
      (\dQuote{t.x.star}), number of transaction opportunities in the
      holdout period (\dQuote{n.star}), and number of customers with the
      recency-frequency combination in the holdout period
      (\dQuote{custs}).}
    \item{$x.star}{A vector with 22 elements, containing the number of
      transactions made by each calibration period recency-frequency bin
      in the holdout period. It is in the same order as $rf.matrix.}
    \item{$annual.sales}{A vector with 11 elements, containing the
      number of transactions made by all customers in each time period
      in both the calibration and holdout periods.}
  }
}

\source{
  Data can be found online at \url{http://www.brucehardie.com/notes/010/}
  (Associated Excel spreadsheet)
}

\details{
  
  Data from \dQuote{a major nonprofit organization located in the midwestern
  United States that is funded in large part by donations from
  individuals. In 1995 the organization \dQuote{acquired} 11,104 first-time
  supporters; in each of the following six years, these individuals
  either did or did not support the organization.}
  
  This dataset contains, for each possible in-sample recency/frequency combination
  in the 1995 cohort, the number of customers and the number of
  transactions they made during the validation period.
  
}


\references{
  
  Fader, Peter S., Bruce G.S. Hardie, and Jen
  Shang. \dQuote{Customer-Base Analysis in a Discrete-Time
  Noncontractual Setting.} \emph{Marketing Science} 29(6),
  pp. 1086-1108. 2010. INFORMS. \url{http://www.brucehardie.com/papers/020/}

}
\keyword{datasets}
