################################################## BG/NBD estimation, visualization functions
library(hypergeo)

# Two things discovered in this script so far:
# -- bgnbd.cbs.LL should be called with the un-compressed version of cal.cbs, the 3-column one
# -- bgnbd.LL spec, as written, won't avoid the large x problem. Patched that, not tested yet.

#' Define general parameters
#'
#' This is to ensure consistency across all functions that require common bits
#' and bobs.
#' 
#' @inheritParams bgnbd.LL
#' @inheritParams bgnbd.ConditionalExpectedTransactions
#' @param func function calling dc.InputCheck
#' @param hardie if TRUE, use \code{\link{h2f1}} instead of
#'   \code{\link[hypergeo]{hypergeo}} when you call this function from within
#'   \code{\link{bgnbd.ConditionalExpectedTransactions}}.
#' @return a list with things you need for \code{\link{bgnbd.LL}},
#'   \code{\link{bgnbd.PAlive}} and
#'   \code{\link{bgnbd.ConditionalExpectedTransactions}}
#' @seealso \code{\link{bgnbd.LL}}
#' @seealso \code{\link{bgnbd.PAlive}}
#' @seealso \code{\link{bgnbd.ConditionalExpectedTransactions}}
bgnbd.generalParams <- function(params, 
                                func,
                                x, 
                                t.x, 
                                T.cal, 
                                T.star = NULL, 
                                hardie = NULL) {
  inputs <- try(dc.InputCheck(params = params, 
                              func = func, 
                              printnames = c("r", "alpha", "a", "b"),
                              x = x, 
                              t.x = t.x, 
                              T.cal = T.cal))
  if('try-error' == class(inputs)) return(inputs)
  
  x <- inputs$x
  t.x <- inputs$t.x
  T.cal <- inputs$T.cal
  
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  
  # last two components for the alt specification
  # to handle large values of x (Solution #2 in
  # http://brucehardie.com/notes/027/bgnbd_num_error.pdf, 
  # LL specification (4) on page 4):
  C3 = ((alpha + t.x)/(alpha + T.cal))^(r + x)
  C4 = a / (b + x - 1)
  
  # stuff you'll need in sundry places
  out <- list()
  out$PAlive <- 1/(1 + as.numeric(x > 0) * C4 / C3)
  
  # do these computations only if needed: that is,
  # if you call this function from bgnbd.LL
  if(func == 'bgnbd.LL') {
    
    # a helper for specifying the log form of the ratio of betas
    # in http://brucehardie.com/notes/027/bgnbd_num_error.pdf
    lb.ratio = function(a, b, x, y) {
      (lgamma(a) + lgamma(b) - lgamma(a + b)) - 
        (lgamma(x) + lgamma(y) - lgamma(x + y))
    }
    
    # First two components -- D1 and D2 -- for the alt spec
    # that can handle large values of x (Solution #2 in
    # http://brucehardie.com/notes/027/bgnbd_num_error.pdf)
    # Here is the D1 term of LL function (4) on page 4:
    D1 = lgamma(r + x) - 
      lgamma(r) + 
      lgamma(a + b) + 
      lgamma(b + x) - 
      lgamma(b) - 
      lgamma(a + b + x)
    D2 = r * log(alpha) - (r + x) * log(alpha + t.x)
    
    # original implementation of the log likelihood
    # A = D2 + lgamma(r + x) - lgamma(r)
    # B = exp(lb.ratio(a, b + x, a, b)) * 
    #   C3 + 
    #   as.numeric((x > 0)) * 
    #   exp(lb.ratio(a + 1, b + x - 1, a, b))
    # out$LL = sum(A + log(B))
    
    # with the corection for avoiding the NUM! problem: 
    out$LL = D1 + D2 + log(C3 + as.numeric((x > 0)) * C4)
  }
  
  # if T.star is not null, then this can produce 
  # conditional expected transactions too. this is 
  # another way of saying that you are calling this
  # function from bgnbd.ConditionalExpectedTransactions, 
  # in which case you also need to set hardie to TRUE or FALSE
  if(!is.null(T.star)) {
    stopifnot(hardie %in% c(TRUE, FALSE))
    term1 <- (a + b + x - 1) / (a - 1)
    if(hardie == TRUE) {
      hyper <- h2f1(r + x, 
                    b + x, 
                    a + b + x - 1, 
                    T.star/(alpha + T.cal + T.star))
    } else {
      hyper <- Re(hypergeo(r + x, 
                           b + x, 
                           a + b + x - 1, 
                           T.star/(alpha + T.cal + T.star)))
    }
    term2 <- 1 - 
      ((alpha + T.cal)/(alpha + T.cal + T.star))^(r + x) * 
      hyper
    out$CET <- term1 * term2 * out$PAlive
  }
  out
}

#' BG/NBD Log-Likelihood
#'
#' Calculates the log-likelihood of the BG/NBD model.
#'
#' \code{x}, \code{t.x} and \code{T.cal} may be vectors. The standard rules for
#' vector operations apply - if they are not of the same length, shorter vectors
#' will be recycled (start over at the first element) until they are as long as
#' the longest vector. It is advisable to keep vectors to the same length and to
#' use single values for parameters that are to be the same for all
#' calculations. If one of these parameters has a length greater than one, the
#' output will be also be a vector.
#'
#' @param params BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param x number of repeat transactions in the calibration period T.cal, or a
#'   vector of transaction frequencies.
#' @param t.x time of most recent repeat transaction, or a vector of recencies.
#' @param T.cal length of calibration period, or a vector of calibration period
#'   lengths.
#'
#' @seealso \code{\link{bgnbd.EstimateParameters}}
#' @seealso \code{\link{bgnbd.cbs.LL}}
#'
#' @return A vector of log-likelihoods as long as the longest input vector (x,
#'   t.x, or T.cal).
#'
#' @examples
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' # random assignment of parameters
#' params <- c(0.5, 6, 1.2, 3.3)
#' # returns the log-likelihood of the given parameters
#' bgnbd.cbs.LL (params, cal.cbs)
#'
#' # compare the speed and results to the following:
#' cal.cbs.compressed <- dc.compress.cbs(cal.cbs)
#' bgnbd.cbs.LL(params, cal.cbs.compressed)
#'
#' # Returns the log likelihood of the parameters for a customer who
#' # made 3 transactions in a calibration period that ended at t=6,
#' # with the last transaction occurring at t=4.
#' bgnbd.LL(params, x=3, t.x=4, T.cal=6)
#'
#' # We can also give vectors as function parameters:
#' set.seed(7)
#' x <- sample(1:4, 10, replace = TRUE)
#' t.x <- sample(1:4, 10, replace = TRUE)
#' T.cal <- rep(4, 10)
#' bgnbd.LL(params, x, t.x, T.cal)
bgnbd.LL <- function(params, 
                     x, 
                     t.x, 
                     T.cal) {
  bgnbd.generalParams(params, 
                      'bgnbd.LL', 
                      x, 
                      t.x, 
                      T.cal)$LL
}

#' BG/NBD Log-Likelihood Wrapper
#'
#' Calculates the log-likelihood sum of the BG/NBD model.
#'
#' Note: do not use a compressed \code{cal.cbs} matrix. It makes quicker work
#' for Pareto/NBD estimation as implemented in this package, but the opposite is
#' true for BG/NBD. For proof, compare the definition of the
#' \code{\link{bgnbd.cbs.LL}} to that of \code{\link{pnbd.cbs.LL}}.
#'
#' @param params BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period. If your data is compressed (see \code{\link{dc.compress.cbs}}),
#'   a fourth column labeled "custs" (number of customers with a specific
#'   combination of recency, frequency and length of calibration period) is
#'   available.
#'
#' @seealso \code{\link{bgnbd.EstimateParameters}}
#' @seealso \code{\link{bgnbd.LL}}
#'
#' @return The total log-likelihood of the provided data.
#'
#' @examples
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' # random assignment of parameters
#' params <- c(0.5, 6, 1.2, 3.3)
#' # returns the log-likelihood of the given parameters
#' bgnbd.cbs.LL(params, cal.cbs)
#'
#' # compare the speed and results to the following:
#' cal.cbs.compressed <- dc.compress.cbs(cal.cbs)
#' bgnbd.cbs.LL (params, cal.cbs.compressed)
#'
#' # Returns the log likelihood of the parameters for a customer who
#' # made 3 transactions in a calibration period that ended at t=6,
#' # with the last transaction occurring at t=4.
#' bgnbd.LL(params, x=3, t.x=4, T.cal=6)
#'
#' # We can also give vectors as function parameters:
#' set.seed(7)
#' x <- sample(1:4, 10, replace = TRUE)
#' t.x <- sample(1:4, 10, replace = TRUE)
#' T.cal <- rep(4, 10)
#' bgnbd.LL(params, x, t.x, T.cal)
bgnbd.cbs.LL <- function(params, 
                         cal.cbs) {
  dc.check.model.params(printnames = c("r", "alpha", "a", "b"), 
                        params = params, 
                        func = "bgnbd.cbs.LL")
  # Check that you have the right columns.
  # They should be 'x', 't.x', 'T.cal' and optionally 'custs.'
  # They stand for, respectively:
  # -- x: frequency
  # -- t.x: recency
  # -- T.cal: observed calendar time
  # -- custs: number of customers with this (x, t.x, T.cal) combo
  foo <- colnames(cal.cbs)
  stopifnot(all(c('x', 't.x', 'T.cal') %in% foo))
  x <- cal.cbs[,'x']
  t.x <- cal.cbs[,'t.x']
  T.cal <- cal.cbs[,'T.cal']
  
  # Avoid this unfurling exercise by calling bgnbd.cbs.LL 
  # with the uncompressed version of cal.cbs, which doesn't 
  # have a "custs" column.
  if ("custs" %in% colnames(cal.cbs)) {
    many_rows = function(vec, nreps) {
      return(rep(1, nreps) %*% t.default(vec))
    }
    custs <- cal.cbs[, "custs"]
    logvec = (1:length(custs)) * (custs > 1)
    logvec = logvec[logvec > 0]
    M = sum(logvec > 0)
    for (i in 1:M) {
      cal.cbs = rbind(cal.cbs, 
                      many_rows(cal.cbs[logvec[i], ], 
                                custs[logvec[i]] - 1))
    }
    x = cal.cbs[, "x"]
    t.x = cal.cbs[, "t.x"]
    T.cal = cal.cbs[, "T.cal"]
  }
  return(sum(bgnbd.LL(params, x, t.x, T.cal)))
}

#' BG/NBD Parameter Estimation
#'
#' Estimates parameters for the BG/NBD model.
#'
#' The best-fitting parameters are determined using the
#' \code{\link{bgnbd.cbs.LL}} function. The sum of the log-likelihood for each
#' customer (for a set of parameters) is maximized in order to estimate
#' parameters.
#'
#' A set of starting parameters must be provided for this method. If no
#' parameters are provided, (1,3,1,3) is used as a default. These values are
#' used because they provide good convergence across data sets. It may be useful
#' to use starting values for r and alpha that represent your best guess of the
#' heterogeneity in the buy and die rate of customers. It may be necessary to
#' run the estimation from multiple starting points to ensure that it converges.
#' To compare the log-likelihoods of different parameters, use
#' \code{\link{bgnbd.cbs.LL}}.
#'
#' The lower bound on the parameters to be estimated is always zero, since
#' BG/NBD parameters cannot be negative. The upper bound can be set with the
#' max.param.value parameter.
#'
#' This function may take some time to run.
#'
#' @param cal.cbs	calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period.
#' @param par.start	initial BG/NBD parameters - a vector with r, alpha, a, and
#'   b, in that order. r and alpha are unobserved parameters for the NBD
#'   transaction process. a and b are unobserved parameters for the Beta
#'   geometric dropout process.
#' @param max.param.value	the upper bound on parameters.
#' @param method the optimization method(s) passed along to
#'   \code{\link[optimx]{optimx}}.
#' @param hessian set it to TRUE if you want the Hessian matrix, and then you
#'   might as well have the complete  \code{\link[optimx]{optimx}} object
#'   returned.
#' @return Vector of estimated parameters.
#' @seealso \code{\link{bgnbd.cbs.LL}}
#' @references Fader, Peter S.; Hardie, and Bruce G.S.. "Overcoming the BG/NBD
#'   Model's #NUM! Error Problem." December. 2013. Web.
#'   \url{http://brucehardie.com/notes/027/bgnbd_num_error.pdf}
#'
#' @examples
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' # starting-point parameters
#' startingparams <- c(1.0, 3, 1.0, 3)
#'
#' # estimated parameters
#' est.params <- bgnbd.EstimateParameters(cal.cbs = cal.cbs,
#'                                        par.start = startingparams)
#'
#' # complete object returned by \code{\link[optimx]{optimx}}
#' optimx.set <- bgnbd.EstimateParameters(cal.cbs = cal.cbs,
#'                                        par.start = startingparams,
#'                                        hessian = TRUE)
#'
#' # log-likelihood of estimated parameters
#' bgnbd.cbs.LL(est.params, cal.cbs)
bgnbd.EstimateParameters <- function(cal.cbs, 
                                     par.start = c(1, 3, 1, 3), 
                                     max.param.value = 10000, 
                                     method = 'L-BFGS-B',
                                     hessian = FALSE) {
  dc.check.model.params(printnames = c("r", "alpha", "a", "b"), 
                        params = par.start, 
                        func = "bgnbd.EstimateParameters")
  bgnbd.eLL <- function(params, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] = max.param.value
    return(-1 * bgnbd.cbs.LL(params, cal.cbs))
  }
  logparams = log(par.start)
  
  results <- optimx(par = logparams, 
                    fn = bgnbd.eLL, 
                    cal.cbs = cal.cbs,
                    max.param.value = max.param.value, 
                    method = method, 
                    hessian = hessian)
  if(hessian == TRUE) {
    message('Your parameter estimates are now on a log scale. Exponentiate them before use.')
    return(results)
  }  
  unlist(exp(results[method, c('p1', 'p2', 'p3', 'p4')]))
}

#' BG/NBD Conditional Expected Transactions
#'
#' E\[X(T.cal, T.cal + T.star) | x, t.x, r, alpha, a, b\]
#'
#' \code{T.star}, \code{x}, \code{t.x} and \code{T.cal} may be vectors. The
#' standard rules for vector operations apply - if they are not of the same
#' length, shorter vectors will be recycled (start over at the first element)
#' until they are as long as the longest vector. It is advisable to keep vectors
#' to the same length and to use single values for parameters that are to be the
#' same for all calculations. If one of these parameters has a length greater
#' than one, the output will be a vector of probabilities.
#'
#' @inheritParams bgnbd.LL
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param hardie if TRUE, use \code{\link{h2f1}} instead of
#'   \code{\link[hypergeo]{hypergeo}}.
#' @return Number of transactions a customer is expected to make in a time
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected
#'   number of transactions.
#' @seealso \code{\link{bgnbd.Expectation}}
#' @references Fader, Peter S.; Hardie, Bruce G.S.and Lee, Ka Lok. “Computing
#'   P(alive) Using the BG/NBD Model.” December. 2008. Web.
#'   \url{http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf}
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' # Number of transactions a customer is expected to make in 2 time
#' # intervals, given that they made 10 repeat transactions in a time period
#' # of 39 intervals, with the 10th repeat transaction occurring in the 35th
#' # interval.
#' bgnbd.ConditionalExpectedTransactions(params, T.star=2, x=10, t.x=35, T.cal=39)
#'
#' # We can also compare expected transactions across different
#' # calibration period behaviors:
#' bgnbd.ConditionalExpectedTransactions(params, T.star=2, x=5:20, t.x=25, T.cal=39)
bgnbd.ConditionalExpectedTransactions <- function(params, 
                                                  T.star, 
                                                  x, 
                                                  t.x, 
                                                  T.cal, 
                                                  hardie = TRUE) {
  bgnbd.generalParams(params, 
                      'bgnbd.ConditionalExpectedTransactions', 
                      x, 
                      t.x, 
                      T.cal, 
                      T.star, 
                      hardie)$CET
}

#' BG/NBD Expectation
#'
#' Returns the number of repeat transactions that a randomly chosen customer
#' (for whom we have no prior information) is expected to make in a given time
#' period.
#'
#' E(X(t) | r, alpha, a, b)
#'
#' @param params BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param t length of time for which we are calculating the expected number of
#'   repeat transactions.
#' @param hardie if TRUE, use \code{\link{h2f1}} instead of
#'   \code{\link[hypergeo]{hypergeo}}.
#' @return Number of repeat transactions a customer is expected to make in a
#'   time period of length t.
#' @seealso \code{\link{bgnbd.ConditionalExpectedTransactions}}
#' @references Fader, Peter S.; Hardie, Bruce G.S.and Lee, Ka Lok. “Computing
#'   P(alive) Using the BG/NBD Model.” December. 2008. Web.
#'   \url{http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf}
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#'
#' # Number of repeat transactions a customer is expected to make in 2 time intervals.
#' bgnbd.Expectation(params, t=2, hardie = FALSE)
#'
#' # We can also compare expected transactions over time:
#' bgnbd.Expectation(params, t=1:10)
bgnbd.Expectation <- function(params, 
                              t, 
                              hardie = TRUE) {
  dc.check.model.params(printnames = c("r", "alpha", "a", "b"), 
                        params = params, 
                        func = "bgnbd.Expectation")
  if (any(t < 0) || !is.numeric(t)) 
    stop("t must be numeric and may not contain negative numbers.")
  r = params[1]
  alpha = params[2]
  a = params[3]
  b = params[4]
  
  term1 = (a + b - 1)/(a - 1)
  term2 = (alpha/(alpha + t))^r
  if(hardie == TRUE) {
    term3 = h2f1(r, b, a + b - 1, t/(alpha + t))
  } else {
    term3 = Re(hypergeo(r, b, a + b - 1, t/(alpha + t)))
  }
  output = term1 * (1 - term2 * term3)
  
  return(output)
}

#' BG/NBD Probability Mass Function
#'
#' Probability mass function for the BG/NBD.
#'
#' P(X(t)=x | r, alpha, a, b). Returns the probability that a customer makes x
#' repeat transactions in the time interval (0, t].
#'
#' Parameters t and x may be vectors. The standard rules for vector operations
#' apply - if they are not of the same length, the shorter vector will be
#' recycled (start over at the first element) until it is as long as the longest
#' vector. It is advisable to keep vectors to the same length and to use single
#' values for parameters that are to be the same for all calculations. If one of
#' these parameters has a length greater than one, the output will be a vector
#' of probabilities.
#'
#' @param params  BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param t 	length end of time period for which probability is being computed.
#'   May also be a vector.
#' @param x   number of repeat transactions by a random customer in the period
#'   defined by t. May also be a vector.
#' @return Probability of X(t)=x conditional on model parameters. If t and/or x
#'   has a length greater than one, a vector of probabilities will be returned.
#' @references Fader, Peter S.; Hardie, Bruce G.S.and Lee, Ka Lok. “Computing
#'   P(alive) Using the BG/NBD Model.” December. 2008. 
#'   [Web.](http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf)
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' # probability that a customer will make 10 repeat transactions in the
#' # time interval (0,2]
#' bgnbd.pmf(params, t=2, x=10)
#' # probability that a customer will make no repeat transactions in the
#' # time interval (0,39]
#' bgnbd.pmf(params, t=39, x=0)
#'
#' # Vectors may also be used as arguments:
#' bgnbd.pmf(params, t=30, x=11:20)
#' @md
bgnbd.pmf <- function(params, 
                      t, 
                      x) {
  inputs <- try(dc.InputCheck(params = params, 
                              func = 'bgnbd.pmf', 
                              printnames = c("r", "alpha", "a", "b"),
                              x = x, 
                              t = t))
  if('try-error' == class(inputs)) return(inputs)
  return(bgnbd.pmf.General(params, 
                           t.start = 0, 
                           t.end = inputs$t, 
                           x = inputs$x))
}

#' Generalized BG/NBD Probability Mass Function
#'
#' Generalized probability mass function for the BG/NBD.
#'
#' P(X(t.start, t.end)=x | r, alpha, a, b). Returns the probability that a
#' customer makes x repeat transactions in the time interval (t.start, t.end\].
#'
#' It is impossible for a customer to make a negative number of repeat
#' transactions. This function will return an error if it is given negative
#' times or a negative number of repeat transactions. This function will also
#' return an error if t.end is less than t.start.
#'
#' t.start, t.end, and x may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be a
#' vector of probabilities.
#'
#' @param params 	BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param t.start 	start of time period for which probability is being
#'   calculated. It can also be a vector of values.
#' @param t.end   end of time period for which probability is being calculated.
#'   It can also be a vector of values.
#' @param x 	number of repeat transactions by a random customer in the period
#'   defined by (t.start, t.end]. It can also be a vector of values.
#' @return Probability of x transaction occuring between t.start and t.end
#'   conditional on model parameters. If t.start, t.end, and/or x has a length
#'   greater than one, a vector of probabilities will be returned.
#' @references Fader, Peter S.; Hardie, Bruce G.S.and Lee, Ka Lok. “Computing
#'   P(alive) Using the BG/NBD Model.” December. 2008. 
#'   [Web.](http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf)
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' # probability that a customer will make 10 repeat transactions in the
#' # time interval (1,2]
#' bgnbd.pmf.General(params, t.start=1, t.end=2, x=10)
#' # probability that a customer will make no repeat transactions in the
#' # time interval (39,78]
#' bgnbd.pmf.General(params, t.start=39, t.end=78, x=0)
#' @md
bgnbd.pmf.General <- function(params, 
                              t.start, 
                              t.end, 
                              x) {
  inputs <- try(dc.InputCheck(params = params, 
                              func = 'bgnbd.pmf.General', 
                              printnames = c("r", "alpha", "a", "b"),
                              t.start = t.start, 
                              t.end = t.end, 
                              x = x))
  if('try-error' == class(inputs)) return(inputs)
  t.start = inputs$t.start
  t.end = inputs$t.end
  x = inputs$x
  max.length <- nrow(inputs)
  
  if (any(t.start > t.end)) {
    stop("Error in bgnbd.pmf.General: t.start > t.end.")
  }
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  equation.part.0 <- rep(0, max.length)
  t = t.end - t.start
  term3 = rep(0, max.length)
  term1 = beta(a, b + x)/beta(a, b) * 
    gamma(r + x)/gamma(r)/factorial(x) * 
    ((alpha/(alpha + t))^r) * ((t/(alpha + t))^x)
  
  for (i in 1:max.length) {
    if (x[i] > 0) {
      ii = c(0:(x[i] - 1))
      summation.term = sum(gamma(r + ii)/gamma(r)/factorial(ii) * 
                             ((t[i]/(alpha + t[i]))^ii))
      term3[i] = 1 - (((alpha/(alpha + t[i]))^r) * summation.term)
    }
  }
  term2 = as.numeric(x > 0) * beta(a + 1, b + x - 1)/beta(a, b) * term3
  return(term1 + term2)
}

#' BG/NBD P(Alive)
#'
#' Uses BG/NBD model parameters and a customer's past transaction behavior to
#' return the probability that they are still alive at the end of the
#' calibration period.
#'
#' P(Alive | X=x, t.x, T.cal, r, alpha, a, b)
#'
#' x, t.x, and T.cal may be vectors. The standard rules for vector operations
#' apply - if they are not of the same length, shorter vectors will be recycled
#' (start over at the first element) until they are as long as the longest
#' vector. It is advisable to keep vectors to the same length and to use single
#' values for parameters that are to be the same for all calculations. If one of
#' these parameters has a length greater than one, the output will be a vector
#' of probabilities.
#'
#' @inheritParams bgnbd.LL
#' @return Probability that the customer is still alive at the end of the
#'   calibration period. If x, t.x, and/or T.cal has a length greater than one,
#'   then this will be a vector of probabilities (containing one element
#'   matching each element of the longest input vector).
#' @references Fader, Peter S.; Hardie, Bruce G.S.and Lee, Ka Lok. “Computing
#'   P(alive) Using the BG/NBD Model.” December. 2008. 
#'   [Web.](http://www.brucehardie.com/notes/021/palive_for_BGNBD.pdf)
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' 
#' bgnbd.PAlive(params, x=23, t.x=39, T.cal=39)
#' # P(Alive) of a customer who has the same recency and total
#' # time observed.
#' 
#' bgnbd.PAlive(params, x=5:20, t.x=30, T.cal=39)
#' # Note the "increasing frequency paradox".
#' 
#' # To visualize the distribution of P(Alive) across customers:
#' 
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' params <- bgnbd.EstimateParameters(cbs, par.start = c(0.243, 4.414, 0.793, 2.426))
#' p.alives <- bgnbd.PAlive(params, cbs[,"x"], cbs[,"t.x"], cbs[,"T.cal"])
#' plot(density(p.alives))
#' @md
bgnbd.PAlive <- function(params, 
                         x, 
                         t.x, 
                         T.cal) {
  bgnbd.generalParams(params, 'bgnbd.PAlive', x, t.x, T.cal)$PAlive
}

#' BG/NBD Expected Cumulative Transactions
#'
#' Calculates the expected cumulative total repeat transactions by all customers
#' for the calibration and holdout periods.
#'
#' The function automatically divides the total period up into n.periods.final
#' time intervals. n.periods.final does not have to be in the same unit of time
#' as the T.cal data. For example: - if your T.cal data is in weeks, and you
#' want cumulative transactions per week, n.periods.final would equal T.star. -
#' if your T.cal data is in weeks, and you want cumulative transactions per day,
#' n.periods.final would equal T.star * 7.
#'
#' The holdout period should immediately follow the calibration period. This
#' function assume that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using max(T.cal) - T.cal rather than assuming that it is 0).
#'
#' @param params 	BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param T.cal   a vector to represent customers' calibration period lengths
#'   (in other words, the "T.cal" column from a customer-by-sufficient-statistic
#'   matrix).
#' @param T.tot   end of holdout period. Must be a single value, not a vector.
#' @param n.periods.final   number of time periods in the calibration and
#'   holdout periods. See details.
#' @param hardie  if TRUE, use h2f1 instead of hypergeo.
#' @return Vector of expected cumulative total repeat transactions by all
#'   customers.
#' @seealso [`bgnbd.Expectation`]
#' @examples
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' params <- c(0.243, 4.414, 0.793, 2.426)
#'
#' # Returns a vector containing cumulative repeat transactions for 273 days.
#' # All parameters are in weeks; the calibration period lasted 39 weeks.
#' bgnbd.ExpectedCumulativeTransactions(params,
#'                                      T.cal = cal.cbs[,"T.cal"],
#'                                      T.tot = 39,
#'                                      n.periods.final = 273,
#'                                      hardie = TRUE)
#' @md
bgnbd.ExpectedCumulativeTransactions <- function(params, 
                                                 T.cal, 
                                                 T.tot, 
                                                 n.periods.final, 
                                                 hardie) {
  
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = "bgnbd.ExpectedCumulativeTransactions")
  
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
    stop("T.cal must be a single numeric value and may not be negative.")
  if (length(n.periods.final) > 1 || n.periods.final < 0 || !is.numeric(n.periods.final)) 
    stop("n.periods.final must be a single numeric value and may not be negative.")
  
  intervals <- seq(T.tot/n.periods.final, 
                   T.tot, 
                   length.out = n.periods.final)
  
  cust.birth.periods <- max(T.cal) - T.cal
  
  expected.transactions <- sapply(intervals, 
                                  function(interval) {
                                    if (interval <= min(cust.birth.periods)) return(0)
                                    t <- interval - cust.birth.periods[cust.birth.periods <= interval]
                                    sum(bgnbd.Expectation(params = params, 
                                                          t = t, 
                                                          hardie = hardie))
                                  })
  
  return(expected.transactions)
}

#' BG/NBD Plot Frequency in Calibration Period
#'
#' Plots a histogram and returns a matrix comparing the actual and expected
#' number of customers who made a certain number of repeat transactions in the
#' calibration period, binned according to calibration period frequencies.
#'
#' This function requires a censor number, which cannot be higher than the
#' highest frequency in the calibration period CBS. The output matrix will have
#' (censor + 1) bins, starting at frequencies of 0 transactions and ending at a
#' bin representing calibration period frequencies at or greater than the censor
#' number. The plot may or may not include a bin for zero frequencies, depending
#' on the plotZero parameter.
#'
#' @param params BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x") and total time observed ("T.cal").
#' @param censor integer used to censor the data. See details.
#' @param plotZero If FALSE, the histogram will exclude the zero bin.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param title title placed on the top-center of the plot.
#' @return Calibration period repeat transaction frequency comparison matrix
#'   (actual vs. expected).
#' @examples 
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs 
#' # cal.cbs already has column names required by method
#'
#' # parameters estimated using bgnbd.EstimateParameters 
#' est.params <- c(0.243, 4.414, 0.793, 2.426) 
#' # the maximum censor number that can be used
#' max(cal.cbs[,"x"])
#'
#' bgnbd.PlotFrequencyInCalibration(est.params, cal.cbs, censor=7)
bgnbd.PlotFrequencyInCalibration <- function(params, 
                                             cal.cbs, 
                                             censor, 
                                             plotZero = TRUE, 
                                             xlab = "Calibration period transactions", 
                                             ylab = "Customers", 
                                             title = "Frequency of Repeat Transactions") {
  
  tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.PlotFrequencyInCalibration: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.PlotFrequencyInCalibration: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotFrequencyInCalibration")
  if (censor > max(x)) 
    stop("censor too big (> max freq) in PlotFrequencyInCalibration.")
  
  x = cal.cbs[, "x"]
  T.cal = cal.cbs[, "T.cal"]
  n.x <- rep(0, max(x) + 1)
  ncusts = nrow(cal.cbs)
  for (ii in unique(x)) {
    # Get number of customers to buy n.x times, over the grid of all possible n.x
    # values (no censoring)
    n.x[ii + 1] <- sum(ii == x)
  }
  n.x.censor <- sum(n.x[(censor + 1):length(n.x)])
  n.x.actual <- c(n.x[1:censor], n.x.censor)  # This upper truncates at censor (ie. if censor=7, 8 categories: {0, 1, ..., 6, 7+}).
  T.value.counts <- table(T.cal)  # This is the table of counts of all time durations from customer birth to end of calibration period.
  T.values <- as.numeric(names(T.value.counts))  # These are all the unique time durations from customer birth to end of calibration period.
  n.T.values <- length(T.values)  # These are the number of time durations we need to consider.
  n.x.expected <- rep(0, length(n.x.actual))  # We'll store the probabilities in here.
  n.x.expected.all <- rep(0, max(x) + 1)  # We'll store the probabilities in here.
  
  for (ii in 0:max(x)) {
    # We want to run over the probability of each transaction amount.
    this.x.expected = 0
    for (T.idx in 1:n.T.values) {
      # We run over all people who had all time durations.
      T = T.values[T.idx]
      if (T == 0) 
        next
      n.T = T.value.counts[T.idx]  # This is the number of customers who had this time duration.
      prob.of.this.x.for.this.T = bgnbd.pmf(params, T, ii)
      expected.given.x.and.T = n.T * prob.of.this.x.for.this.T
      this.x.expected = this.x.expected + expected.given.x.and.T
    }
    n.x.expected.all[ii + 1] = this.x.expected
  }
  n.x.expected[1:censor] = n.x.expected.all[1:censor]
  n.x.expected[censor + 1] = sum(n.x.expected.all[(censor + 1):(max(x) + 1)])
  
  col.names <- paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
  col.names[censor + 1] <- paste(col.names[censor + 1], "+", sep = "")
  censored.freq.comparison <- rbind(n.x.actual, n.x.expected)
  colnames(censored.freq.comparison) <- col.names
  cfc.plot <- censored.freq.comparison
  if (plotZero == FALSE) 
    cfc.plot <- cfc.plot[, -1]
  n.ticks <- ncol(cfc.plot)
  if (plotZero == TRUE) {
    x.labels <- 0:(n.ticks - 1)
    x.labels[n.ticks] <- paste(n.ticks - 1, "+", sep = "")
  }
  ylim <- c(0, ceiling(max(cfc.plot) * 1.1))
  barplot(cfc.plot, names.arg = x.labels, beside = TRUE, ylim = ylim, main = title, 
          xlab = xlab, ylab = ylab, col = 1:2)
  legend("topright", legend = c("Actual", "Model"), col = 1:2, lwd = 2)
  return(censored.freq.comparison)
}

#' BG/NBD Plot Frequency vs. Conditional Expected Frequency
#'
#' Plots the actual and conditional expected number transactions made by
#' customers in the holdout period, binned according to calibration period
#' frequencies. Also returns a matrix with this comparison and the number of
#' customers in each bin.
#'
#' This function requires a censor number, which cannot be higher than the
#' highest frequency in the calibration period CBS. The output matrix will have
#' (censor + 1) bins, starting at frequencies of 0 transactions and ending at a
#' bin representing calibration period frequencies at or greater than the censor
#' number.
#' 
#' @param params 	BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param T.star  length of then holdout period.
#' @param cal.cbs   calibration period CBS (customer by sufficient statistic).
#'   It must contain columns for frequency ("x"), recency ("t.x"), and total
#'   time observed ("T.cal"). Note that recency must be the time between the
#'   start of the calibration period and the customer's last transaction, not
#'   the time between the customer's last transaction and the end of the
#'   calibration period.
#' @param x.star  vector of transactions made by each customer in the holdout
#'   period.
#' @param censor  integer used to censor the data. See details.
#' @param xlab  descriptive label for the x axis.
#' @param ylab  descriptive label for the y axis.
#' @param xticklab  vector containing a label for each tick mark on the x axis.
#' @param title   title placed on the top-center of the plot.
#' @return  Holdout period transaction frequency comparison matrix (actual vs.
#'   expected).
#' @examples
#' data(cdnowSummary)
#' 
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#' 
#' # number of transactions by each customer in the 39 weeks
#' # following the calibration period
#' x.star <- cal.cbs[,"x.star"]
#' 
#' # parameters estimated using bgnbd.EstimateParameters
#' est.params <- c(0.243, 4.414, 0.793, 2.426)
#' # the maximum censor number that can be used
#' max(cal.cbs[,"x"])
#' 
#' # plot conditional expected holdout period frequencies,
#' # binned according to calibration period frequencies
#' bgnbd.PlotFreqVsConditionalExpectedFrequency(est.params, 
#'                                              T.star = 39, 
#'                                              cal.cbs, 
#'                                              x.star, 
#'                                              censor = 7)
bgnbd.PlotFreqVsConditionalExpectedFrequency <- function(params, 
                                                         T.star, 
                                                         cal.cbs, 
                                                         x.star, 
                                                         censor, 
                                                         xlab = "Calibration period transactions", 
                                                         ylab = "Holdout period transactions", 
                                                         xticklab = NULL, 
                                                         title = "Conditional Expectation") {
  tryCatch(x <- cal.cbs[, "x"], 
           error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs[, "t.x"], 
           error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], 
           error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotFreqVsConditionalExpectedFrequency")
  if (censor > max(x)) 
    stop("censor too big (> max freq) in PlotFreqVsConditionalExpectedFrequency.")
  
  if (any(T.star < 0) || !is.numeric(T.star)) 
    stop("T.star must be numeric and may not contain negative numbers.")
  if (any(x.star < 0) || !is.numeric(x.star)) 
    stop("x.star must be numeric and may not contain negative numbers.")
  
  n.bins = censor + 1
  transaction.actual = rep(0, n.bins)
  transaction.expected = rep(0, n.bins)
  bin.size = rep(0, n.bins)
  for (cc in 0:censor) {
    if (cc != censor) {
      this.bin = which(cc == x)
    } else if (cc == censor) {
      this.bin = which(x >= cc)
    }
    n.this.bin = length(this.bin)
    bin.size[cc + 1] = n.this.bin
    transaction.actual[cc + 1] = sum(x.star[this.bin])/n.this.bin
    transaction.expected[cc + 1] = sum(bgnbd.ConditionalExpectedTransactions(params, 
                                                                             T.star, x[this.bin], t.x[this.bin], T.cal[this.bin]))/n.this.bin
  }
  col.names = paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
  col.names[censor + 1] = paste(col.names[censor + 1], "+", sep = "")
  comparison = rbind(transaction.actual, transaction.expected, bin.size)
  colnames(comparison) = col.names
  if (is.null(xticklab) == FALSE) {
    x.labels = xticklab
  }
  if (is.null(xticklab) != FALSE) {
    if (censor < ncol(comparison)) {
      x.labels = 0:(censor)
      x.labels[censor + 1] = paste(censor, "+", sep = "")
    }
    if (censor >= ncol(comparison)) {
      x.labels = 0:(ncol(comparison))
    }
  }
  actual = comparison[1, ]
  expected = comparison[2, ]
  ylim = c(0, ceiling(max(c(actual, expected)) * 1.1))
  plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, 
       main = title)
  lines(expected, lty = 2, col = 2)
  axis(1, at = 1:ncol(comparison), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  return(comparison)
}

#' BG/NBD Plot Actual vs. Conditional Expected Frequency by Recency
#'
#' Plots the actual and conditional expected number of transactions made by
#' customers in the holdout period, binned according to calibration period
#' recencies. Also returns a matrix with this comparison and the number of
#' customers in each bin.
#'
#' This function does bin customers exactly according to recency; it bins
#' customers according to integer units of the time period of cal.cbs.
#' Therefore, if you are using weeks in your data, customers will be binned as
#' follows: customers with recencies between the start of the calibration period
#' (inclusive) and the end of week one (exclusive); customers with recencies
#' between the end of week one (inclusive) and the end of week two (exclusive);
#' etc.
#'
#' The matrix and plot will contain the actual number of transactions made by
#' each bin in the holdout period, as well as the expected number of
#' transactions made by that bin in the holdout period, conditional on that
#' bin's behavior during the calibration period.
#'
#' @inheritParams bgnbd.PlotFreqVsConditionalExpectedFrequency
#' @return Matrix comparing actual and conditional expected transactions in the
#'   holdout period.
#' @examples 
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs 
#' # cal.cbs already has column names required by method
#'
#' # number of transactions by each customer in the 39 weeks following 
#' # the calibration period 
#' x.star <- cal.cbs[,"x.star"]
#'
#' # parameters estimated using bgnbd.EstimateParameters 
#' est.params <- c(0.243, 4.414, 0.793, 2.426)
#'
#' # plot conditional expected holdout period transactions, 
#' # binned according to calibration period recencies
#' bgnbd.PlotRecVsConditionalExpectedFrequency(est.params, 
#'                                             cal.cbs, 
#'                                             T.star = 39,
#'                                             x.star)
bgnbd.PlotRecVsConditionalExpectedFrequency <- function(params, 
                                                        cal.cbs, 
                                                        T.star, 
                                                        x.star, 
                                                        xlab = "Calibration period recency", 
                                                        ylab = "Holdout period transactions", 
                                                        xticklab = NULL, 
                                                        title = "Actual vs. Conditional Expected Transactions by Recency") {
  
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotRecVsConditionalExpectedFrequency")
  
  if (any(T.star < 0) || !is.numeric(T.star)) 
    stop("T.star must be numeric and may not contain negative numbers.")
  if (any(x.star < 0) || !is.numeric(x.star)) 
    stop("x.star must be numeric and may not contain negative numbers.")
  
  tryCatch(x <- cal.cbs[, "x"], 
           error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs[, "t.x"], 
           error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], 
           error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  
  t.values <- sort(unique(t.x))
  n.recs <- length(t.values)
  transaction.actual <- rep(0, n.recs)
  transaction.expected <- rep(0, n.recs)
  rec.size <- rep(0, n.recs)
  
  for (tt in 1:n.recs) {
    this.t.x <- t.values[tt]
    this.rec <- which(t.x == this.t.x)
    n.this.rec <- length(this.rec)
    rec.size[tt] <- n.this.rec
    transaction.actual[tt] <- sum(x.star[this.rec])/n.this.rec
    transaction.expected[tt] <- sum(bgnbd.ConditionalExpectedTransactions(params, 
                                                                          T.star, x[this.rec], t.x[this.rec], T.cal[this.rec]))/n.this.rec
  }
  
  comparison <- rbind(transaction.actual, transaction.expected, rec.size)
  colnames(comparison) <- round(t.values, 3)
  
  bins <- seq(1, ceiling(max(t.x)))
  n.bins <- length(bins)
  actual <- rep(0, n.bins)
  expected <- rep(0, n.bins)
  bin.size <- rep(0, n.bins)
  
  x.labels <- NULL
  if (is.null(xticklab) == FALSE) {
    x.labels <- xticklab
  } else {
    x.labels <- 1:(n.bins)
  }
  point.labels <- rep("", n.bins)
  point.y.val <- rep(0, n.bins)
  for (ii in 1:n.bins) {
    if (ii < n.bins) {
      this.bin <- which(as.numeric(colnames(comparison)) >= (ii - 1) & as.numeric(colnames(comparison)) < 
                          ii)
    } else if (ii == n.bins) {
      this.bin <- which(as.numeric(colnames(comparison)) >= ii - 1)
    }
    actual[ii] <- sum(comparison[1, this.bin])/length(comparison[1, this.bin])
    expected[ii] <- sum(comparison[2, this.bin])/length(comparison[2, this.bin])
    bin.size[ii] <- sum(comparison[3, this.bin])
  }
  
  ylim <- c(0, ceiling(max(c(actual, expected)) * 1.1))
  plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, 
       main = title)
  lines(expected, lty = 2, col = 2)
  
  axis(1, at = 1:n.bins, labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  
  return(rbind(actual, expected, bin.size))
}

#' BG/NBD Plot Transaction Rate Heterogeneity
#'
#' Plots and returns the estimated gamma distribution of lambda (customers'
#' propensities to purchase).
#'
#' This returns the distribution of each customer's Poisson parameter, which
#' determines the rate at which each customer buys.
#'
#' @param params BG/NBD parameters - a vector with r, alpha, a, and b, in that
#'   order. r and alpha are unobserved parameters for the NBD transaction
#'   process. a and b are unobserved parameters for the Beta geometric dropout
#'   process.
#' @param lim upper-bound of the x-axis. A number is chosen by the function if
#'   none is provided.
#' @return Distribution of customers' propensities to purchase.
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' bgnbd.PlotTransactionRateHeterogeneity(params)
#' params <- c(0.53, 4.414, 0.793, 2.426)
#' bgnbd.PlotTransactionRateHeterogeneity(params)
bgnbd.PlotTransactionRateHeterogeneity <- function(params, 
                                                   lim = NULL) {
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotTransactionRateHeterogeneity")
  shape <- params[1]
  rate <- params[2]
  rate.mean <- round(shape/rate, 4)
  rate.var <- round(shape/rate^2, 4)
  if (is.null(lim)) {
    lim = qgamma(0.99, shape = shape, rate = rate)
  }
  x.axis.ticks <- seq(0, lim, length.out = 100)
  heterogeneity <- dgamma(x.axis.ticks, shape = shape, rate = rate)
  plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Transaction Rate", ylab = "Density", 
       main = "Heterogeneity in Transaction Rate")
  mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
  mtext(mean.var.label, side = 3)
  return(rbind(x.axis.ticks, heterogeneity))
}

#' BG/NBD Plot Dropout Probability Heterogeneity
#'
#' Plots and returns the estimated gamma distribution of p (customers'
#' probability of dropping out immediately after a transaction).
#'
#' @inheritParams bgnbd.PlotTransactionRateHeterogeneity
#' @return Distribution of customers' probabilities of dropping out.
#' @examples
#' params <- c(0.243, 4.414, 0.793, 2.426)
#' bgnbd.PlotDropoutRateHeterogeneity(params)
#' params <- c(0.243, 4.414, 1.33, 2.426)
#' bgnbd.PlotDropoutRateHeterogeneity(params)
bgnbd.PlotDropoutRateHeterogeneity <- function(params, 
                                               lim = NULL) {
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotDropoutRateHeterogeneity")
  alpha_param = params[3]
  beta_param = params[4]
  beta_param.mean = round(alpha_param/(alpha_param + beta_param), 4)
  beta_param.var = round(alpha_param * beta_param/((alpha_param + beta_param)^2)/(alpha_param + 
                                                                                    beta_param + 1), 4)
  if (is.null(lim)) {
    # get right end point of grid
    lim = qbeta(0.99, shape1 = alpha_param, shape2 = beta_param)
  }
  x.axis.ticks = seq(0, lim, length.out = 100)
  heterogeneity = dbeta(x.axis.ticks, shape1 = alpha_param, shape2 = beta_param)
  plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Dropout Probability p", 
       ylab = "Density", main = "Heterogeneity in Dropout Probability")
  mean.var.label = paste("Mean:", beta_param.mean, "    Var:", beta_param.var)
  mtext(mean.var.label, side = 3)
  return(rbind(x.axis.ticks, heterogeneity))
}

#' BG/NBD Tracking Cumulative Transactions Plot
#'
#' Plots the actual and expected cumulative total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' actual.cu.tracking.data does not have to be in the same unit of time as the
#' T.cal data. T.tot will automatically be divided into periods to match the
#' length of actual.cu.tracking.data. See
#' [bgnbd.ExpectedCumulativeTransactions].
#'
#' The holdout period should immediately follow the calibration period. This
#' function assume that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using max(T.cal) - T.cal rather than assuming that it is 0).
#'
#' @inheritParams bgnbd.ExpectedCumulativeTransactions
#' @inheritParams bgnbd.PlotFreqVsConditionalExpectedFrequency
#' @param actual.cu.tracking.data vector containing the cumulative number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods). See details.
#' @return Matrix containing actual and expected cumulative repeat transactions.
#' @examples
#' data(cdnowSummary)
#'
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' # Cumulative repeat transactions made by all customers across calibration
#' # and holdout periods
#' cu.tracking <- cdnowSummary$cu.tracking
#'
#' # parameters estimated using bgnbd.EstimateParameters
#' est.params <- c(0.243, 4.414, 0.793, 2.426)
#'
#' # All parameters are in weeks; the calibration period lasted 39
#' # weeks and the holdout period another 39.
#' bgnbd.PlotTrackingCum(est.params, 
#'                       T.cal = cal.cbs[,"T.cal"], 
#'                       T.tot = 78, 
#'                       actual.cu.tracking.data = cu.tracking, 
#'                       hardie = TRUE)
#' @md
bgnbd.PlotTrackingCum <- function(params, 
                                  T.cal, 
                                  T.tot, 
                                  actual.cu.tracking.data, 
                                  n.periods.final = NA,
                                  hardie,
                                  xlab = "Week", 
                                  ylab = "Cumulative Transactions", 
                                  xticklab = NULL, 
                                  title = "Tracking Cumulative Transactions") {
  
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.Plot.PlotTrackingCum")
  
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  if (any(actual.cu.tracking.data < 0) || !is.numeric(actual.cu.tracking.data)) 
    stop("actual.cu.tracking.data must be numeric and may not contain negative numbers.")
  
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
    stop("T.cal must be a single numeric value and may not be negative.")
  
  actual <- actual.cu.tracking.data
  if(is.na(n.periods.final)) n.periods.final <- length(actual)
  expected <- bgnbd.ExpectedCumulativeTransactions(params, 
                                                   T.cal, 
                                                   T.tot, 
                                                   n.periods.final, 
                                                   hardie)
  
  cu.tracking.comparison <- rbind(actual, expected)
  
  ylim <- c(0, max(c(actual, expected)) * 1.05)
  plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
       main = title)
  lines(expected, lty = 2, col = 2)
  if (is.null(xticklab) == FALSE) {
    if (ncol(cu.tracking.comparison) != length(xticklab)) {
      stop("Plot error, xticklab does not have the correct size")
    }
    axis(1, at = 1:ncol(cu.tracking.comparison), labels = xticklab)
  } else {
    axis(1, at = 1:length(actual), labels = 1:length(actual))
  }
  abline(v = max(T.cal), lty = 2)
  
  legend("bottomright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  
  return(cu.tracking.comparison)
}

#' BG/NBD Tracking Incremental Transactions Comparison
#'
#' Plots the actual and expected incremental total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' actual.inc.tracking.data does not have to be in the same unit of time as the
#' T.cal data. T.tot will automatically be divided into periods to match the
#' length of actual.inc.tracking.data. See
#' [bgnbd.ExpectedCumulativeTransactions].
#'
#' The holdout period should immediately follow the calibration period. This
#' function assume that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using max(T.cal) - T.cal rather than assuming that it is 0).
#'
#' @inheritParams bgnbd.PlotTrackingCum
#' @param actual.inc.tracking.data  vector containing the incremental number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods). See details.
#' @return Matrix containing actual and expected incremental repeat
#'   transactions.
#' @examples
#' data(cdnowSummary)
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#'
#' # Cumulative repeat transactions made by all customers across calibration
#' # and holdout periods
#' cu.tracking <- cdnowSummary$cu.tracking
#' # make the tracking data incremental
#' inc.tracking <- dc.CumulativeToIncremental(cu.tracking)
#'
#' # parameters estimated using bgnbd.EstimateParameters
#' est.params <- c(0.243, 4.414, 0.793, 2.426)
#'
#' # All parameters are in weeks; the calibration period lasted 39
#' # weeks and the holdout period another 39.
#' bgnbd.PlotTrackingInc(est.params, ,
#'                       T.cal = cal.cbs[,"T.cal"],
#'                       T.tot = 78,
#'                       actual.inc.tracking.data = inc.tracking, 
#'                       hardie = TRUE)
#' @md
bgnbd.PlotTrackingInc <- function(params, 
                                  T.cal, 
                                  T.tot, 
                                  actual.inc.tracking.data, 
                                  n.periods.final = NA,
                                  hardie, 
                                  xlab = "Week", 
                                  ylab = "Transactions", 
                                  xticklab = NULL, 
                                  title = "Tracking Weekly Transactions") {
  
  dc.check.model.params(printnames = c("r", "alpha", "a", "b"), 
                        params = params, 
                        func = "bgnbd.Plot.PlotTrackingCum")
  
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  if (any(actual.inc.tracking.data < 0) || !is.numeric(actual.inc.tracking.data)) 
    stop("actual.inc.tracking.data must be numeric and may not contain negative numbers.")
  
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
    stop("T.cal must be a single numeric value and may not be negative.")
  
  actual <- actual.inc.tracking.data
  if(is.na(n.periods.final)) n.periods.final <- length(actual)
  expected <- dc.CumulativeToIncremental(bgnbd.ExpectedCumulativeTransactions(params, 
                                                                              T.cal, 
                                                                              T.tot, 
                                                                              n.periods.final, 
                                                                              hardie))
  
  ylim <- c(0, max(c(actual, expected)) * 1.05)
  plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
       main = title)
  lines(expected, lty = 2, col = 2)
  if (is.null(xticklab) == FALSE) {
    if (length(actual) != length(xticklab)) {
      stop("Plot error, xticklab does not have the correct size")
    }
    axis(1, at = 1:length(actual), labels = xticklab)
  } else {
    axis(1, at = 1:length(actual), labels = 1:length(actual))
  }
  abline(v = max(T.cal), lty = 2)
  
  legend("topright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  
  return(rbind(actual, expected))
}



