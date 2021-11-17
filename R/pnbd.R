################################################## Pareto/NBD estimation, visualization functions

library(hypergeo)
library(optimx)

#' Define general parameters
#'
#' This is to ensure consistency across all functions that require the
#' likelihood function, or the log of it, and to make sure that the same
#' implementation of the hypergeometric function is used everywhere for building
#' \code{A0}.
#'
#' This function is only ever called by either \code{\link{pnbd.LL}} or
#' \code{\link{pnbd.PAlive}} so it returns directly the output that is expected
#' from those calling functions: either the log likelihood for a set of
#' customers, or the probability that a set of customers with characteristics
#' given by \code{x}, \code{t.x} and \code{T.cal}, having estimated a set of
#' \code{params}, is still alive. Either set of customers can be of size 1.
#' @inheritParams pnbd.LL
#' @param func name of the function calling dc.InputCheck; either \code{pnbd.LL}
#'   or \code{pnbd.PAlive}.
#' @return A vector of log likelihood values if \code{func} is \code{pnbd.LL},
#'   or a vector of probabilities that a customer is still alive if \code{func}
#'   is \code{pnbd.PAlive}.
#' @seealso \code{\link{pnbd.LL}}
#' @seealso \code{\link{pnbd.PAlive}}
#' @seealso \code{\link{pnbd.DERT}}
pnbd.generalParams <- function(params, 
                               x, 
                               t.x, 
                               T.cal, 
                               func,
                               hardie = TRUE) {
  # Since pnbd.LL and pnbd.pAlive are the only options
  # for func, we don't need a printnames argument
  # in the pnbd.generalParams wrapper.
  stopifnot(func %in% c('pnbd.LL', 'pnbd.PAlive'))
  inputs <- try(dc.InputCheck(params = params, 
                                    func = func, 
                                    printnames = c("r", "alpha", "s", "beta"),
                                    x = x, 
                                    t.x = t.x, 
                                    T.cal = T.cal))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  
  x <- inputs$x
  t.x <- inputs$t.x
  T.cal <- inputs$T.cal
  
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  
  maxab <- max(alpha, beta)
  absab <- abs(alpha - beta)
  
  param2 <- s + 1
  if (alpha < beta) {
    param2 <- r + x
  }
  
  a <- alpha + T.cal
  b <- maxab + t.x
  c <- beta + T.cal
  d <- maxab + T.cal
  w <- r + s + x
  
  if(hardie == TRUE) {
    F1 <- h2f1(a = w, 
               b = param2, 
               c = w + 1, 
               z = absab / b) 
    F2 <- h2f1(a = w, 
               b = param2, 
               c = w + 1, 
               z = absab / d) 
  } else {
    F1 <- Re(hypergeo(A = w, 
                      B = param2, 
                      C = w + 1, 
                      z = absab / b))
    F2 <- Re(hypergeo(A = w, 
                      B = param2, 
                      C = w + 1, 
                      z = absab / d))
  }
  A0 <- F1/(b^w) - F2/(d^w)
  
  # You only ever call this function from two other 
  # places: pnbd.LL or pnbd.PAlive.
  if(func == 'pnbd.LL') {
    # this returns the log likelihood for one random customer
    part1 <- r * log(alpha) + 
      s * log(beta) + 
      lgamma(r + x) - 
      lgamma(r)
    part2 <- 1 / (a^(w - s) * c^s)
    return(part1 + log(part2) + log(1 + (s/w) * A0 / part2))
  }  
  else if(func == 'pnbd.PAlive') {
    # This returns the probability that a random customer is still alive
    return(1 / (1 + s/w * a^(w - s) * c^s * A0))
  } else {
    return(NULL)
  }
}

#' Use Bruce Hardie's Gaussian hypergeometric implementation
#'
#' In benchmarking \code{\link{pnbd.LL}} runs more quickly and 
#' it returns the same results if it uses this helper instead of 
#' \code{\link[hypergeo]{hypergeo}}, which is the default. But \code{h2f1} 
#' is such a barebones function that in some edge cases it will keep
#' going until you get a segfault, where \code{\link[hypergeo]{hypergeo}} 
#' would have failed with a proper error message.
#'
#' @param a, counterpart to A in \code{\link[hypergeo]{hypergeo}}
#' @param b, counterpart to B in \code{\link[hypergeo]{hypergeo}}
#' @param c, counterpart to C in \code{\link[hypergeo]{hypergeo}}
#' @param z, counterpart to z in \code{\link[hypergeo]{hypergeo}}
#' @seealso \code{\link[hypergeo]{hypergeo}}
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the Pareto/NBD Model and 
#' Related Expressions." November. 2005. Web. \url{http://www.brucehardie.com/notes/008/}
h2f1 <- function(a, b, c, z) {
  lenz <- length(z)
  j = 0
  uj <- 1:lenz
  uj <- uj/uj
  y <- uj
  lteps <- 0
  
  while (lteps < lenz) {
    lasty <- y
    j <- j + 1
    uj <- uj * (a + j - 1) * (b + j - 1)/(c + j - 1) * z/j
    y <- y + uj
    lteps <- sum(y == lasty)
  }
  return(y)
}

#' Pareto/NBD Log-Likelihood
#' 
#' Calculates the log-likelihood of the Pareto/NBD model.
#'
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param x number of repeat transactions in the calibration period T.cal, or a
#'   vector of transaction frequencies.
#' @param t.x time of most recent repeat transaction, or a vector of recencies.
#' @param T.cal length of calibration period, or a vector of calibration period
#'   lengths.
#' @param hardie if TRUE, use \code{\link{h2f1}} instead of
#'   \code{\link[hypergeo]{hypergeo}}.
#'
#' @seealso \code{\link{pnbd.EstimateParameters}}
#'
#' @return A vector of log-likelihoods as long as the longest input vector (x,
#'   t.x, or T.cal).
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the
#'   Pareto/NBD Model and Related Expressions." November. 2005. Web.
#'   \url{http://www.brucehardie.com/notes/008/}
#' 
#' @examples
#' # Returns the log likelihood of the parameters for a customer who
#' # made 3 transactions in a calibration period that ended at t=6,
#' # with the last transaction occurring at t=4.
#' pnbd.LL(params, x=3, t.x=4, T.cal=6, hardie = TRUE)
#' 
#' # We can also give vectors as function parameters:
#' set.seed(7)
#' x <- sample(1:4, 10, replace = TRUE)
#' t.x <- sample(1:4, 10, replace = TRUE)
#' T.cal <- rep(4, 10)
#' pnbd.LL(params, x, t.x, T.cal, hardie = TRUE)
pnbd.LL <- function(params, 
                    x, 
                    t.x, 
                    T.cal, 
                    hardie = TRUE) {
  pnbd.generalParams(params = params, 
                     x = x, 
                     t.x = t.x, 
                     T.cal = T.cal, 
                     func = 'pnbd.LL', 
                     hardie = hardie)
}

#' Pareto/NBD P(Alive)
#'
#' Uses Pareto/NBD model parameters and a customer's past transaction behavior
#' to return the probability that they are still alive at the end of the
#' calibration period.
#'
#' P(Alive | X=x, t.x, T.cal, r, alpha, s, beta)
#'
#' x, t.x, and T.cal may be vectors. The standard rules for vector operations
#' apply - if they are not of the same length, shorter vectors will be recycled
#' (start over at the first element) until they are as long as the longest
#' vector. It is advisable to keep vectors to the same length and to use single
#' values for parameters that are to be the same for all calculations. If one of
#' these parameters has a length greater than one, the output will be a vector
#' of probabilities.
#'
#' @inheritParams pnbd.LL
#'
#' @return Probability that the customer is still alive at the end of the
#'   calibration period. If x, t.x, and/or T.cal has a length greater than one,
#'   then this will be a vector of probabilities (containing one element
#'   matching each element of the longest input vector).
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the
#'   Pareto/NBD Model and Related Expressions." November. 2005. Web.
#'   \url{http://www.brucehardie.com/notes/008/}
#' 
#' @examples 
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' params <- pnbd.EstimateParameters(cbs, hardie = TRUE)
#' 
#' pnbd.PAlive(params, x=0, t.x=0, T.cal=39, TRUE)
#' # 0.2941633; P(Alive) of a customer who made no repeat transactions.
#' 
#' pnbd.PAlive(params, x=23, t.x=39, T.cal=39, TRUE)
#' # 1; P(Alive) of a customer who has the same recency and total
#' # time observed.
#' 
#' pnbd.PAlive(params, x=5:20, t.x=30, T.cal=39, TRUE)
#' # Note the "increasing frequency paradox".
#' 
#' # To visualize the distribution of P(Alive) across customers:
#' p.alives <- pnbd.PAlive(params, cbs[,"x"], cbs[,"t.x"], cbs[,"T.cal"], TRUE)
#' plot(density(p.alives))
pnbd.PAlive <- function(params, 
                        x, 
                        t.x, 
                        T.cal, 
                        hardie = TRUE) {
  pnbd.generalParams(params = params, 
                     x = x, 
                     t.x = t.x, 
                     T.cal = T.cal, 
                     func = 'pnbd.PAlive', 
                     hardie = hardie)
}

#' Pareto/NBD Discounted Expected Residual Transactions
#'
#' Calculates the discounted expected residual transactions of a customer, given
#' their behavior during the calibration period.
#'
#' DERT(d | r, alpha, s, beta, X = x, t.x, T.cal)
#'
#' x, t.x, T.cal may be vectors. The standard rules for vector operations apply
#' - if they are not of the same length, shorter vectors will be recycled (start
#' over at the first element) until they are as long as the longest vector. It
#' is advisable to keep vectors to the same length and to use single values for
#' parameters that are to be the same for all calculations. If one of these
#' parameters has a length greater than one, the output will be also be a
#' vector.
#'
#' @inheritParams pnbd.LL
#' @param d the discount rate to be used. Make sure that it matches up with your
#'   chosen time period (do not use an annual rate for monthly data, for
#'   example).
#'
#' @return The number of discounted expected residual transactions for a
#'   customer with a particular purchase pattern during the calibration period.
#' @references Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. "RFM and CLV:
#'   Using Iso-Value Curves for Customer Base Analysis." Journal of Marketing
#'   Research Vol.42, pp.415-430. November. 2005.
#'   \url{http://www.brucehardie.com/papers.html}
#' @references See equation 2.
#' @references Note that this paper refers to what this package is calling
#'   discounted expected residual transactions (DERT) simply as discounted
#'   expected transactions (DET).
#' 
#' @examples
#' # elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD2"),2,3)
#' # elog[, 'date'] <- as.Date(elog[, 'date'], format = '%Y%m%d')
#' # cal.cbs <- dc.ElogToCbsCbt(elog)$cal$cbs
#' # params <- pnbd.EstimateParameters(cal.cbs, hardie = TRUE)
#' params <- c(0.5629966, 12.5590370, 0.4081095, 10.5148048)
#' 
#' # 15% compounded annually has been converted to 0.0027 compounded continuously,
#' # as we are dealing with weekly data and not annual data.
#' d <- 0.0027
#' 
#' # calculate the discounted expected residual transactions of a customer
#' # who made 7 transactions in a calibration period that was 77.86
#' # weeks long, with the last transaction occurring at the end of
#' # the 35th week.
#' pnbd.DERT(params, 
#'           x = 7, 
#'           t.x = 35, 
#'           T.cal = 77.86, 
#'           d, 
#'           hardie = TRUE)
#' 
#' # We can also use vectors to compute DERT for several customers:
#' pnbd.DERT(params, 
#'           x = 1:10, 
#'           t.x = 30, 
#'           T.cal = 77.86, 
#'           d, 
#'           hardie = TRUE)
pnbd.DERT <- function(params, 
                      x, 
                      t.x, 
                      T.cal, 
                      d, 
                      hardie = TRUE) {
  loglike <- try(pnbd.LL(params = params, 
                         x = x, 
                         t.x = t.x, 
                         T.cal = T.cal, 
                         hardie = hardie))
  if('try-error' %in% class(loglike)) return(loglike)
  
  # This is the remainder of the original pnbd.DERT function def. 
  # No need to get too clever here. Revert to explicit assignment 
  # of params to r, alpha, s, beta the old-school way.
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  z <- d * (beta + T.cal)
  
  tricomi.part.1 = ((z)^(1 - s))/(s - 1) * 
    genhypergeo(U = c(1), 
                L = c(2 - s), 
                z = z, 
                check_mod = FALSE)
  
  tricomi.part.2 = gamma(1 - s) * 
    genhypergeo(U = c(s), 
                L = c(s), 
                z = z, 
                check_mod = FALSE)
  
  tricomi = tricomi.part.1 + tricomi.part.2
  
  result <- exp(r * log(alpha) + 
                  s * log(beta) + 
                  (s - 1) * log(d) + 
                  lgamma(r + x + 1) + 
                  log(tricomi) - 
                  lgamma(r) - 
                  (r + x + 1) * log(alpha + T.cal) - 
                  loglike)
  return(result)
}

#' Pareto/NBD Log-Likelihood
#'
#' Calculates the log-likelihood of the Pareto/NBD model.
#'
#' @inheritParams pnbd.LL
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period. If your data is compressed (see \code{\link{dc.compress.cbs}}), a
#'   fourth column labelled "custs" (number of customers with a specific
#'   combination of recency, frequency and length of calibration period) will
#'   make this function faster.
#' @seealso \code{\link{pnbd.EstimateParameters}}
#' @seealso \code{\link{pnbd.LL}}
#' @return The log-likelihood of the provided data.
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the
#'   Pareto/NBD Model and Related Expressions." November. 2005. Web.
#'   \url{http://www.brucehardie.com/notes/008/}
#'   
#' @examples
#' data(cdnowSummary)
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#' 
#' # random assignment of parameters
#' params <- c(0.5, 8, 0.7, 10)
#' # returns the log-likelihood of the given parameters
#' pnbd.cbs.LL (params, cal.cbs, TRUE)
#' 
#' # compare the speed and results to the following:
#' cal.cbs.compressed <- dc.compress.cbs(cal.cbs)
#' pnbd.cbs.LL (params, cal.cbs.compressed, TRUE)
pnbd.cbs.LL <- function(params, 
                        cal.cbs, 
                        hardie = TRUE) {
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = "pnbd.cbs.LL")
  # Check that you have the right columns.
  # They should be 'x', 't.x', 'T.cal' and optionally 'custs'
  # in this order. They stand for, respectively
  # -- x: frequency
  # -- t.x: recency
  # -- T.cal: observed calendar time
  # -- custs: number of customers with this (x, t.x, T.cal) combo
  foo <- colnames(cal.cbs)
  stopifnot(foo[1] == 'x' & 
              foo[2] == 't.x' & 
              foo[3] == 'T.cal')
  x <- cal.cbs[,'x']
  t.x <- cal.cbs[,'t.x']
  T.cal <- cal.cbs[,'T.cal']

  if ("custs" %in% foo) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * pnbd.LL(params, x, t.x, T.cal, hardie)))
}

#' Pareto/NBD Parameter Estimation
#'
#' The best-fitting parameters are determined using the
#' \code{\link{pnbd.cbs.LL}} function. The sum of the log-likelihood for each
#' customer (for a set of parameters) is maximized in order to estimate
#' parameters.
#'
#' A set of starting parameters must be provided for this method. If no
#' parameters are provided, (1,1,1,1) is used as a default. It may be useful to
#' use starting values for r and s that represent your best guess of the
#' heterogeneity in the buy and die rate of customers. It may be necessary to
#' run the estimation from multiple starting points to ensure that it converges.
#' To compare the log-likelihoods of different parameters, use
#' \code{\link{pnbd.cbs.LL}}.
#'
#' The lower bound on the parameters to be estimated is always zero, since
#' Pareto/NBD parameters cannot be negative. The upper bound can be set with the
#' max.param.value parameter.
#'
#' This function may take some time to run. It uses \code{\link[optimx]{optimx}} 
#' for maximum likelihood estimation, not \code{\link[stats]{optim}}.
#'
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period. If your data is compressed (see \code{\link{dc.compress.cbs}}), a
#'   fourth column labelled "custs" (number of customers with a specific
#'   combination of recency, frequency and length of calibration period) will
#'   make this function faster.
#' @param par.start initial Pareto/NBD parameters - a vector with r, alpha, s,
#'   and beta, in that order. r and alpha are unobserved parameters for the NBD
#'   transaction process. s and beta are unobserved parameters for the Pareto
#'   (exponential gamma) dropout process.
#' @param max.param.value the upper bound on parameters.
#' @param method the optimization method(s).
#' @param hardie if TRUE, have \code{\link{pnbd.LL}} use \code{\link{h2f1}}
#'   instead of \code{\link[hypergeo]{hypergeo}}.
#' @param hessian set it to TRUE if you want the Hessian matrix, and then you
#'   might as well have the complete  \code{\link[optimx]{optimx}} object
#'   returned.
#'
#' @return Unnamed vector of estimated parameters by default, \code{optimx}
#'   object with everything if \code{hessian} is TRUE.
#' @seealso \code{\link{pnbd.cbs.LL}}
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
#' startingparams <- c(0.5, 6, 0.9, 8)
#' 
#' # estimated parameters
#' est.params <- pnbd.EstimateParameters(cal.cbs = cal.cbs, 
#'                                       par.start = startingparams, 
#'                                       method = 'L-BFGS-B',
#'                                       hardie = TRUE)
#'                                       
#' # complete object returned by \code{\link[optimx]{optimx}}
#' optimx.set <- pnbd.EstimateParameters(cal.cbs = cal.cbs, 
#'                                       par.start = startingparams, 
#'                                       hardie = TRUE, 
#'                                       hessian = TRUE)
#' 
#' # log-likelihood of estimated parameters
#' pnbd.cbs.LL(est.params, cal.cbs, TRUE)
pnbd.EstimateParameters <- function(cal.cbs, 
                                    par.start = c(1, 1, 1, 1), 
                                    max.param.value = 10000, 
                                    method = 'L-BFGS-B', 
                                    hardie = TRUE, 
                                    hessian = FALSE) {
    dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                          params = par.start, 
                          func = "pnbd.EstimateParameters")
    
    ## helper function to be optimized
    pnbd.eLL <- function(params, cal.cbs, hardie) {
        params <- exp(params)
        params[params > max.param.value] <- max.param.value
        return(-1 * pnbd.cbs.LL(params = params, 
                                cal.cbs = cal.cbs, 
                                hardie = hardie))
    }
    logparams <- log(par.start)
    if(hessian == TRUE) {
      return(optimx(par = logparams, 
                    fn = pnbd.eLL, 
                    method = method,
                    cal.cbs = cal.cbs, 
                    hardie = hardie, 
                    hessian = hessian))
    }
    results <- optimx(par = logparams, 
                      fn = pnbd.eLL, 
                      method = method,
                      cal.cbs = cal.cbs, 
                      hardie = hardie)
    params <- exp(unname(unlist(results[method, c('p1', 'p2', 'p3', 'p4')])))
    return(params)
}

#' Generalized Pareto/NBD Probability Mass Function
#'
#' Generalized probability mass function for the Pareto/NBD.
#'
#' P(X(t.start, t.end)=x | r, alpha, s, beta). Returns the probability that a
#' customer makes x repeat transactions in the time interval (t.start, t.end].
#'
#' It is impossible for a customer to make a negative number of repeat
#' transactions. This function will return an error if it is given negative
#' times or a negative number of repeat transactions. This function will also
#' return an error if t.end is less than t.start.
#'
#' \code{t.start}, \code{t.end}, and \code{x} may be vectors. The standard rules
#' for vector operations apply - if they are not of the same length, shorter
#' vectors will be recycled (start over at the first element) until they are as
#' long as the longest vector. It is advisable to keep vectors to the same
#' length and to use single values for parameters that are to be the same for
#' all calculations. If one of these parameters has a length greater than one,
#' the output will be a vector of probabilities.
#'
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param t.start start of time period for which probability is being
#'   calculated. It can also be a vector of values.
#' @param t.end end of time period for which probability is being calculated. It
#'   can also be a vector of values.
#' @param x number of repeat transactions by a random customer in the period
#'   defined by (t.start, t.end]. It can also be a vector of values.
#' @param hardie if TRUE, use \code{\link{h2f1}} instead of
#'   \code{\link[hypergeo]{hypergeo}}.
#'
#' @return Probability of x transaction occuring between t.start and t.end
#'   conditional on model parameters. If t.start, t.end, and/or x has a length
#'   greater than one, a vector of probabilities will be returned.
#' @references Fader, Peter S., and Bruce G.S. Hardie. "Deriving an Expression
#'   for P (X(t) = x) Under the Pareto/NBD Model." Sept. 2006. Web.
#'   \url{http://www.brucehardie.com/notes/012/}
#' @references Fader, Peter S., Bruce G.S. Hardie, and Kinshuk Jerath. "Deriving
#'   an Expression for P (X(t, t + tau) = x) Under the Pareto/NBD Model." Sept.
#'   2006. Web. \url{http://www.brucehardie.com/notes/013/}
#' 
#' @examples
#' # probability that a customer will make 10 repeat transactions in the
#' # time interval (1,2]
#' data("cdnowSummary")
#' cal.cbs <- cdnowSummary$cbs
#' params <- pnbd.EstimateParameters(cal.cbs = cal.cbs, 
#'                                   method = "L-BFGS-B",
#'                                   hardie = TRUE)
#' pnbd.pmf.General(params, t.start=1, t.end=2, x=10, hardie = TRUE)
#' # probability that a customer will make no repeat transactions in the
#' # time interval (39,78]
#' pnbd.pmf.General(params, 
#'                  t.start = 39, 
#'                  t.end = 78, 
#'                  x = 0, 
#'                  hardie = TRUE)
pnbd.pmf.General <- function(params, 
                             t.start, 
                             t.end, 
                             x, 
                             hardie = TRUE) {
  if (any(t.start > t.end)) {
    stop("Error in pnbd.pmf.General: t.start > t.end.")
  }
  inputs <- try(dc.InputCheck(params = params, 
                              func = "pnbd.pmf.General", 
                              printnames = c("r", "alpha", "s", "beta"), 
                              t.start = t.start, 
                              t.end = t.end, 
                              x = x))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  
  t.start <- inputs$t.start
  t.end <- inputs$t.end
  x <- inputs$x
  max.length <- nrow(inputs)
  
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  
  equation.part.0 <- rep(0, max.length)
  equation.part.0[x == 0] <- 1 - exp(s * log(beta) - s * log(beta + t.start))
  
  ## (t.end - t.start)^x is left outside the exp() in equation.part.1 
  ## because when t.end = t.start and x = 0, exp() gets us into trouble: 
  ##  -- exp(0 * log(0)) = NaN; 
  ##  -- doing directly 0^0 instead gets us 0^0=1
  ##  -- 1 is much better than NaN.
  equation.part.1 <- exp(lgamma(r + x) - 
                           lgamma(r) - 
                           lfactorial(x) + 
                           r * log(alpha) - 
                           r * log(alpha + t.end - t.start) - 
                           x * log(alpha + t.end - t.start) + 
                           s * log(beta) - 
                           s * log(beta + t.end)) * 
    (t.end - t.start)^x
  
  equation.part.2 <- r * log(alpha) + s * log(beta) + lbeta(r + x, s + 1) - lbeta(r, s)
  
  # Marshal the parameters of the hypergeometric and the
  # denominator in the expressions of B1 and B2 shown in
  # http://www.brucehardie.com/notes/013/Pareto_NBD_interval_pmf_rev.pdf
  B1B2 <- function(hardie, 
                   r, 
                   alpha, 
                   s, 
                   beta, 
                   x, 
                   t.start, 
                   t.end = 0, 
                   ii = 0) {
    myalpha <- alpha
    mybeta <- beta + t.start
    maxab <- max(myalpha, mybeta) + t.end
    absab <- abs(myalpha - mybeta)
    
    param2 <- s + 1
    if (myalpha < mybeta) {
      # < same as <= in the case of
      # the hypergeometric, because = 
      # case is trivial.
      param2 <- r + x
    }
    w <- r + s + ii
    
    a <- w
    b <- param2
    c <- r + s + x + 1
    z <- absab/maxab
    den <- maxab^w
    if(hardie == TRUE) return(h2f1(a, b, c, z)/den)
    return(Re(hypergeo(a, b, c, z))/den) 
  }
  
  B.1 <- mapply(x = inputs$x,
                t.start = inputs$t.start, 
                B1B2, 
                hardie = hardie, 
                r = r, 
                alpha = alpha, 
                s = s, 
                beta = beta)
  #B.1 <- B1B2(hardie, r, alpha, s, beta, x, t.start)

  equation.part.2.summation <- rep(NA, max.length)
  ## In the paper, for i=0 we have t^i / i * B(r+s, i). 
  ## the denominator reduces to:
  ## i * Gamma (r+s) * Gamma(i) / Gamma (r+s+i) : 
  ## Gamma (r+s) * Gamma(i+1) / Gamma(r+s+i) : 
  ## Gamma (r+s) * Gamma(1) / Gamma(r+s) : 
  ## 1 
  ## The 1 represents this reduced piece of the equation.
  
  for (i in 1:max.length) {
    ii <- c(1:x[i])
    equation.part.2.summation[i] <- B1B2(hardie, r, alpha, s, beta, x[i], t.start[i], t.end[i], 0)
    if (x[i] > 0) {
      equation.part.2.summation[i] <- equation.part.2.summation[i] + 
        sum((t.end[i] - t.start[i])^ii/(ii * beta(r + s, ii)) * 
              B1B2(hardie, r, alpha, s, beta, x[i], t.start[i], t.end[i], ii))
    }
  }
  return(equation.part.0 + 
           equation.part.1 + 
           exp(equation.part.2 + 
                 log(B.1 - equation.part.2.summation)))
}

#' Pareto/NBD Conditional Expected Transactions
#'
#' Uses Pareto/NBD model parameters and a customer's past transaction behavior
#' to return the number of transactions they are expected to make in a given
#' time period.
#'
#' E\[X(T.cal, T.cal + T.star) | x, t.x, r, alpha, s, beta\]
#'
#' \code{T.star}, \code{x}, \code{t.x}, and \code{T.cal} may be vectors. The
#' standard rules for vector operations apply - if they are not of the same
#' length, shorter vectors will be recycled (start over at the first element)
#' until they are as long as the longest vector. It is advisable to keep vectors
#' to the same length and to use single values for parameters that are to be the
#' same for all calculations. If one of these parameters has a length greater
#' than one, the output will be a vector of probabilities.
#'
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param x number of repeat transactions in the calibration period T.cal, or a
#'   vector of calibration period frequencies.
#' @param t.x time of most recent repeat transaction, or a vector of recencies.
#' @param T.cal length of calibration period, or a vector of calibration period
#'   lengths.
#' @param hardie if TRUE, have \code{\link{pnbd.PAlive}} use \code{\link{h2f1}}
#'   instead of \code{\link[hypergeo]{hypergeo}}.
#'
#' @return Number of transactions a customer is expected to make in a time
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected
#'   number of transactions.
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the
#'   Pareto/NBD Model and Related Expressions." November. 2005. Web.
#'   \url{http://www.brucehardie.com/notes/008/}
#' @seealso \code{\link{pnbd.Expectation}}
#' 
#' @examples
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' # Number of transactions a customer is expected to make in 2 time
#' # intervals, given that they made 10 repeat transactions in a time period
#' # of 39 intervals, with the 10th repeat transaction occurring in the 35th
#' # interval.
#' pnbd.ConditionalExpectedTransactions(params, 
#'                                      T.star = 2, 
#'                                      x = 10, 
#'                                      t.x = 35, 
#'                                      T.cal = 39, 
#'                                      hardie = TRUE)
#' 
#' # We can also compare expected transactions across different
#' # calibration period behaviors:
#' pnbd.ConditionalExpectedTransactions(params, 
#'                                      T.star = 2, 
#'                                      x = 5:20, 
#'                                      t.x = 25, 
#'                                      T.cal = 39, 
#'                                      hardie = TRUE)
pnbd.ConditionalExpectedTransactions <- function(params, 
                                                 T.star, 
                                                 x, 
                                                 t.x, 
                                                 T.cal, 
                                                 hardie = TRUE) {
  inputs <- try(dc.InputCheck(params = params, 
                              func = "pnbd.ConditionalExpectedTransactions", 
                              printnames = c("r", "alpha", "s", "beta"), 
                              T.star = T.star, 
                              x = x, 
                              t.x = t.x, 
                              T.cal = T.cal))
  if('try-error' == class(inputs)) return(str(inputs)$message)

  T.star <- inputs$T.star
  x <- inputs$x
  t.x <- inputs$t.x
  T.cal <- inputs$T.cal
    
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
    
  P1 <- (r + x) * (beta + T.cal)/((alpha + T.cal) * (s - 1))
  P2 <- (1 - ((beta + T.cal)/(beta + T.cal + T.star))^(s - 1))
  P3 <- pnbd.PAlive(params = params, 
                    x = x, 
                    t.x = t.x, 
                    T.cal = T.cal, 
                    hardie = hardie)
  return(P1 * P2 * P3)
}

#' Pareto/NBD Probability Mass Function
#'
#' Probability mass function for the Pareto/NBD.
#'
#' P(X(t)=x | r, alpha, s, beta). Returns the probability that a customer makes
#' x repeat transactions in the time interval (0, t].
#'
#' Parameters \code{t} and \code{x} may be vectors. The standard rules for
#' vector operations apply - if they are not of the same length, the shorter
#' vector will be recycled (start over at the first element) until it is as long
#' as the longest vector. It is advisable to keep vectors to the same length and
#' to use single values for parameters that are to be the same for all
#' calculations. If one of these parameters has a length greater than one, the
#' output will be a vector of probabilities.
#'
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param t length end of time period for which probability is being computed.
#'   May also be a vector.
#' @param x number of repeat transactions by a random customer in the period
#'   defined by t. May also be a vector.
#' @param hardie if TRUE, have \code{\link{pnbd.pmf.General}} use
#'   \code{\link{h2f1}} instead of \code{\link[hypergeo]{hypergeo}}.
#'
#' @return Probability of X(t)=x conditional on model parameters. If t and/or x
#'   has a length greater than one, a vector of probabilities will be returned.
#' @references Fader, Peter S., and Bruce G.S. Hardie. “Deriving an Expression
#'   for P (X(t) = x) Under the Pareto/NBD Model.” Sept. 2006. Web.
#'   \url{http://www.brucehardie.com/notes/012/}
#' 
#' @examples
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' # probability that a customer will make 10 repeat transactions in the
#' # time interval (0,2]
#' pnbd.pmf(params, t=2, x=10, hardie = TRUE)
#' # probability that a customer will make no repeat transactions in the
#' # time interval (0,39]
#' pnbd.pmf(params, t=39, x=0, hardie = TRUE)
#' 
#' # Vectors may also be used as arguments:
#' pnbd.pmf(params = params, 
#'          t = 30, 
#'          x = 11:20, 
#'          hardie = TRUE)
pnbd.pmf <- function(params, 
                     t, 
                     x, 
                     hardie = TRUE) {
  inputs <- try(dc.InputCheck(params = params, 
                              func = "pnbd.pmf", 
                              printnames = c("r", "alpha", "s", "beta"), 
                              t = t, 
                              x = x))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  pnbd.pmf.General(params = params, 
                   t.start = 0, 
                   t.end = inputs$t, 
                   x = inputs$x, 
                   hardie = hardie)
}

#' Pareto/NBD Expectation
#'
#' Returns the number of repeat transactions that a randomly chosen customer
#' (for whom we have no prior information) is expected to make in a given time
#' period.
#'
#' E(X(t) | r, alpha, s, beta)
#'
#' @inheritParams pnbd.LL
#' @param t The length of time for which we are calculating the expected number
#'   of repeat transactions.
#' @return Number of repeat transactions a customer is expected to make in a
#'   time period of length t.
#' @references Fader, Peter S., and Bruce G.S. Hardie. "A Note on Deriving the
#'   Pareto/NBD Model and Related Expressions." November. 2005. Web.
#'   \url{http://www.brucehardie.com/notes/008/}
#' @seealso [`pnbd.ConditionalExpectedTransactions`]
#' @examples 
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' 
#' # Number of repeat transactions a customer is expected to make in 2 time intervals.
#' pnbd.Expectation(params = params, 
#'                  t = 2)
#' 
#' # We can also compare expected transactions over time:
#' pnbd.Expectation(params = params, 
#'                  t = 1:10)
#' @md
pnbd.Expectation <- function(params, t) {
  
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = "pnbd.Expectation")
  
  if (any(t < 0) || !is.numeric(t)) 
    stop("t must be numeric and may not contain negative numbers.")
  
  r = params[1]
  alpha = params[2]
  s = params[3]
  beta = params[4]
  
  return((r * beta)/(alpha * (s - 1)) * (1 - (beta/(beta + t))^(s - 1)))
}

#' Pareto/NBD Expected Cumulative Transactions
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
#' @inheritParams pnbd.LL
#' @param  T.tot End of holdout period. Must be a single value, not a vector.
#' @param n.periods.final Number of time periods in the calibration and holdout
#'   periods. See details.
#' @return Vector of expected cumulative total repeat transactions by all
#'   customers.
#' @seealso [`pnbd.Expectation`]
#' @examples 
#' data(cdnowSummary)
#' 
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#' 
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' 
#' # Returns a vector containing cumulative repeat transactions for 546 days.
#' # All parameters are in weeks; the calibration period lasted 39 weeks
#' # and the holdout period another 39.
#' pnbd.ExpectedCumulativeTransactions(params = params, 
#'                                     T.cal = cal.cbs[,"T.cal"], 
#'                                     T.tot = 78, 
#'                                     n.periods.final = 546)
#' @md
pnbd.ExpectedCumulativeTransactions <- function(params, 
                                                T.cal, 
                                                T.tot, 
                                                n.periods.final) {
  inputs <- try(dc.InputCheck(params = params, 
                                func = "pnbd.ExpectedCumulativeTransactions", 
                                printnames = c("r", "alpha", "s", "beta"), 
                                T.cal = T.cal, 
                                T.tot = T.tot, 
                                n.periods.final = n.periods.final))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  stopit <- "must be a single numeric value and may not be negative."
  if (length(T.tot) > 1) stop(paste("T.tot", stopit, sep = " "))
  if (length(n.periods.final) > 1) stop(paste("n.periods.final", stopit, sep = " "))

  ## Divide up time into equal intervals
  intervals <- seq(T.tot/n.periods.final, 
                   T.tot, 
                   length.out = n.periods.final)

  cust.birth.periods <- max(T.cal) - T.cal

  expected.transactions <- sapply(intervals, 
                                  function(interval) {
                                    if (interval <= min(cust.birth.periods)) 
                                      return(0)
                                    sum(pnbd.Expectation(params, 
                                                         interval - cust.birth.periods[cust.birth.periods <= interval]))
                                  })

  return(expected.transactions)
}

#' Pareto/NBD Plot Frequency in Calibration Period
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
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x") and total time observed ("T.cal").
#' @param censor integer used to censor the data. See details.
#' @param hardie if TRUE, have \code{\link{pnbd.pmf}} use \code{\link{h2f1}}
#'   instead of \code{\link[hypergeo]{hypergeo}}.
#' @param plotZero if FALSE, the histogram will exclude the zero bin.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param title title placed on the top-center of the plot.
#' 
#' @return Calibration period repeat transaction frequency comparison matrix
#'   (actual vs. expected).
#'   
#' @examples 
#' data(cdnowSummary)
#' cal.cbs <- cdnowSummary$cbs
#' # cal.cbs already has column names required by method
#' 
#' # parameters estimated using pnbd.EstimateParameters
#' est.params <- cdnowSummary$est.params
#' # the maximum censor number that can be used
#' max(cal.cbs[,"x"])
#' 
#' pnbd.PlotFrequencyInCalibration(params = est.params, 
#'                                 cal.cbs = cal.cbs, 
#'                                 censor = 7, 
#'                                 hardie = TRUE)
pnbd.PlotFrequencyInCalibration <- function(params, 
                                            cal.cbs, 
                                            censor, 
                                            hardie = TRUE,
                                            plotZero = TRUE, 
                                            xlab = "Calibration period transactions", 
                                            ylab = "Customers", 
                                            title = "Frequency of Repeat Transactions") {
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = "pnbd.PlotFrequencyInCalibration")
    stopifnot("x" %in% colnames(cal.cbs) | "T.cal" %in% colnames(cal.cbs))
    x <- cal.cbs[, "x"]
    T.cal <- cal.cbs[, "T.cal"]
  
    if (censor > max(x)) 
        stop("censor too big (> max freq) in PlotFrequencyInCalibration.")
    
    n.x <- rep(0, max(x) + 1)
    custs = nrow(cal.cbs)
    
    for (ii in unique(x)) {
        n.x[ii + 1] <- sum(ii == x)
    }
    
    n.x.censor <- sum(n.x[(censor + 1):length(n.x)])
    n.x.actual <- c(n.x[1:censor], n.x.censor)
    
    T.value.counts <- table(T.cal)
    T.values <- as.numeric(names(T.value.counts))
    n.T.values <- length(T.values)
    
    total.probability <- 0
    
    n.x.expected <- rep(0, length(n.x.actual))
    
    for (ii in 1:(censor)) {
        this.x.expected <- 0
        for (T.idx in 1:n.T.values) {
            T <- T.values[T.idx]
            if (T == 0) 
                next
            n.T <- T.value.counts[T.idx]
            expected.given.x.and.T <- n.T * pnbd.pmf(params, T, ii - 1, hardie)
            this.x.expected <- this.x.expected + expected.given.x.and.T
            total.probability <- total.probability + expected.given.x.and.T/custs
        }
        n.x.expected[ii] <- this.x.expected
    }
    n.x.expected[censor + 1] <- custs * (1 - total.probability)
    
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
    } else {
        x.labels <- 1:(n.ticks)
        x.labels[n.ticks] <- paste(n.ticks, "+", sep = "")
    }
    
    ylim <- c(0, ceiling(max(cfc.plot) * 1.1))
    barplot(cfc.plot, names.arg = x.labels, beside = TRUE, ylim = ylim, main = title, 
        xlab = xlab, ylab = ylab, col = 1:2)
    
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lwd = 2)
    
    return(censored.freq.comparison)
}

#' Pareto/NBD Plot Frequency vs. Conditional Expected Frequency
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
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param T.star length of the holdout period. It must be a scalar for this
#'   plot's purposes: you have one holdout period of a given length.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period.
#' @param x.star vector of transactions made by each customer in the holdout
#'   period.
#' @param censor integer used to censor the data. See details.
#' @param hardie if TRUE, have
#'   \code{\link{pnbd.ConditionalExpectedTransactions}} use \code{\link{h2f1}}
#'   instead of \code{\link[hypergeo]{hypergeo}}.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @param title title placed on the top-center of the plot.
#' 
#' @return Holdout period transaction frequency comparison matrix (actual vs. expected).
#' 
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
#' # parameters estimated using pnbd.EstimateParameters
#' est.params <- cdnowSummary$est.params
#' # the maximum censor number that can be used
#' max(cal.cbs[,"x"])
#' 
#' # plot conditional expected holdout period frequencies,
#' # binned according to calibration period frequencies
#' pnbd.PlotFreqVsConditionalExpectedFrequency(params = est.params, 
#'                                             T.star = 39, 
#'                                             cal.cbs = cal.cbs, 
#'                                             x.star = x.star, 
#'                                             censor = 7, 
#'                                             hardie = TRUE)
pnbd.PlotFreqVsConditionalExpectedFrequency <- function(params, 
                                                        T.star, 
                                                        cal.cbs, 
                                                        x.star, 
                                                        censor, 
                                                        hardie = TRUE,
                                                        xlab = "Calibration period transactions", 
                                                        ylab = "Holdout period transactions", 
                                                        xticklab = NULL, 
                                                        title = "Conditional Expectation") {
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = "pnbd.PlotFreqVsConditionalExpectedFrequency")
  # Check that you have the right columns.
  # They should be 'x', 't.x', 'T.cal' and optionally 'custs'
  # in this order. They stand for, respectively
  # -- x: frequency
  # -- t.x: recency
  # -- T.cal: observed calendar time
  # -- custs: number of customers with this (x, t.x, T.cal) combo
  foo <- colnames(cal.cbs)
  stopifnot(foo[1] == 'x' & 
              foo[2] == 't.x' & 
              foo[3] == 'T.cal')
  x <- cal.cbs[,'x']
  t.x <- cal.cbs[,'t.x']
  T.cal <- cal.cbs[,'T.cal']  
  
  if (censor > max(x)) 
    stop("censor too big (> max freq) in PlotFreqVsConditionalExpectedFrequency.")
    
  if (length(T.star) > 1 || T.star < 0 || !is.numeric(T.star)) 
    stop("T.star must be a positive scalar.")
  if (any(x.star < 0) || !is.numeric(x.star)) 
    stop("x.star must be numeric and may not contain negative numbers.")
  
  n.bins <- censor + 1
  bin.size <- rep(0, n.bins)
  # First, find the right sequence of 
  # transaction counts. It may have gaps,
  # for which n.this.bin = 0, which is 
  # kinda hard to divide by. So:
  for (cc in 0:censor) {
    if (cc < censor) {
      this.bin <- which(x == cc)
    } else {
      this.bin <- which(x >= cc)
    }
    n.this.bin <- length(this.bin)
    bin.size[cc + 1] <- n.this.bin
  }
  
  # Now you got the right list of net bins:
  # those with at least 1 customer in them.
  bin.size <- bin.size[bin.size > 0]
  n.bins <- length(bin.size)
  names(bin.size) <- names(table(x))[1:n.bins]
  xvals <- as.integer(names(bin.size))
  transaction.actual <- rep(0, n.bins)
  transaction.expected <- rep(0, n.bins)
  
  for(cc in 1:n.bins) {
    n.this.bin <- bin.size[cc]
    if (cc < n.bins) {
      this.bin <- which(x == xvals[cc])
    } else {
      this.bin <- which(x >= xvals[cc])
    }
    transaction.actual[cc] <- sum(x.star[this.bin])/n.this.bin
    transaction.expected[cc] <- sum(pnbd.ConditionalExpectedTransactions(params, 
                                                                         T.star, 
                                                                         x[this.bin], 
                                                                         t.x[this.bin], 
                                                                         T.cal[this.bin], 
                                                                         hardie))/n.this.bin
  }
  col.names <- paste("freq", names(bin.size), sep = ".")
  col.names[n.bins] <- paste(col.names[n.bins], '+', sep = '')
  comparison <- rbind(transaction.actual, transaction.expected, bin.size)
  colnames(comparison) <- col.names
  
  x.labels <- sapply(colnames(comparison), function(x) gsub('freq.', '', x))
  actual <- comparison[1, ]
  expected <- comparison[2, ]
  
  ylim <- c(0, ceiling(max(c(actual, expected)) * 1.1))
  plot(actual, type = "l", xaxt = "n", 
       col = 1, ylim = ylim, xlab = xlab, ylab = ylab, 
       main = title)
  lines(expected, lty = 2, col = 2)
  
  axis(1, at = 1:ncol(comparison), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  return(comparison)
}

#' Pareto/NBD Plot Actual vs. Conditional Expected Frequency by Recency
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
#' between the end of week one (inclusive) and the end of week two (exlusive);
#' etc.
#'
#' The matrix and plot will contain the actual number of transactions made by
#' each bin in the holdout period, as well as the expected number of
#' transactions made by that bin in the holdout period, conditional on that
#' bin's behavior during the calibration period.
#'
#' @param params Pareto/NBD parameters - a vector with r, alpha, s, and beta, in
#'   that order. r and alpha are unobserved parameters for the NBD transaction
#'   process. s and beta are unobserved parameters for the Pareto (exponential
#'   gamma) dropout process.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period.
#' @param T.star length of the holdout period. It must be a scalar for this
#'   plot's purposes: you have one holdout period of a given length.
#' @param x.star vector of transactions made by each customer in the holdout
#'   period.
#' @param hardie if TRUE, have
#'   \code{\link{pnbd.ConditionalExpectedTransactions}} use \code{\link{h2f1}}
#'   instead of \code{\link[hypergeo]{hypergeo}}.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @param title title placed on the top-center of the plot.
#' 
#' @return Matrix comparing actual and conditional expected transactions in the holdout period.
#' 
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
#' # parameters estimated using pnbd.EstimateParameters
#' est.params <- cdnowSummary$est.params
#' 
#' # plot conditional expected holdout period transactions, binned according to
#' # calibration period recencies
#' pnbd.PlotRecVsConditionalExpectedFrequency(params = est.params, 
#'                                            cal.cbs = cal.cbs, 
#'                                            T.star = 39, 
#'                                            x.star = x.star, 
#'                                            hardie = TRUE)
pnbd.PlotRecVsConditionalExpectedFrequency <- function(params, 
                                                       cal.cbs, 
                                                       T.star, 
                                                       x.star, 
                                                       hardie = TRUE,
                                                       xlab = "Calibration period recency", 
                                                       ylab = "Holdout period transactions", 
                                                       xticklab = NULL, 
                                                       title = "Actual vs. Conditional Expected Transactions by Recency") {
  # No use for inputs, other than as error check.
  inputs <- try(dc.InputCheck(params = params, 
                                func = "pnbd.PlotRecVsConditionalExpectedFrequency", 
                                printnames = c("r", "alpha", "s", "beta"), 
                                T.star = T.star, 
                                x.star = x.star))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  if (length(T.star) > 1) 
    stop("T.star must be a scalar.")
  
  # Check that you have the right columns.
  # They should be 'x', 't.x', 'T.cal' and optionally 'custs'
  # in this order. They stand for, respectively
  # -- x: frequency
  # -- t.x: recency
  # -- T.cal: observed calendar time
  # -- custs: number of customers with this (x, t.x, T.cal) combo
  foo <- colnames(cal.cbs)
  stopifnot(foo[1] == 'x' & 
              foo[2] == 't.x' & 
              foo[3] == 'T.cal')
  x <- cal.cbs[,'x']
  t.x <- cal.cbs[,'t.x']
  T.cal <- cal.cbs[,'T.cal']
    
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
    transaction.expected[tt] <- sum(pnbd.ConditionalExpectedTransactions(params, 
                                                                         T.star, 
                                                                         x[this.rec], 
                                                                         t.x[this.rec], 
                                                                         T.cal[this.rec], 
                                                                         hardie))/n.this.rec
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
      this.bin <- which(as.numeric(colnames(comparison)) >= (ii - 1) & 
                          as.numeric(colnames(comparison)) < ii)
    } else if (ii == n.bins) {
      this.bin <- which(as.numeric(colnames(comparison)) >= ii - 1)
    }
    actual[ii] <- sum(comparison[1, this.bin])/length(comparison[1, this.bin])
    expected[ii] <- sum(comparison[2, this.bin])/length(comparison[2, this.bin])
    bin.size[ii] <- sum(comparison[3, this.bin])
  }
  
  ylim <- c(0, ceiling(max(c(actual, expected)) * 1.1))
  plot(actual, 
       type = "l", 
       xaxt = "n", 
       col = 1, 
       ylim = ylim, 
       xlab = xlab, 
       ylab = ylab, 
       main = title)
  lines(expected, 
        lty = 2, 
        col = 2)
  
  axis(1, at = 1:n.bins, labels = x.labels)
  legend("topleft", 
         legend = c("Actual", "Model"), 
         col = 1:2, 
         lty = 1:2, 
         lwd = 1)
  
  return(rbind(actual, expected, bin.size))
}

#' Pareto/NBD Plot Discounted Expected Residual Transactions
#'
#' Plots discounted expected residual transactions for different combinations of
#' calibration period frequency and recency.
#'
#' The length of the calibration period \code{T.cal} must be a single value, not
#' a vector.
#'
#' @inheritParams pnbd.DERT
#' @param type must be either "persp" (perspective - 3 dimensional) or
#'   "contour". Determines the type of plot produced by this function.
#'
#' @return A matrix with discounted expected residual transaction values for
#'   every combination of calibration period frequency \code{x} and calibration
#'   period recency \code{t.x}.
#'
#' @references Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. "RFM and CLV:
#'   Using Iso-Value Curves for Customer Base Analysis." Journal of Marketing
#'   Research Vol.42, pp.415-430. November. 2005.
#'   \url{http://www.brucehardie.com/papers.html}
#' @references Note that this paper refers to what this package is calling
#'   discounted expected residual transactions (DERT) simply as discounted
#'   expected transactions (DET).
#' 
#' @examples
#' # The RFM and CLV paper uses all 78 weeks of the cdnow data to
#' # estimate parameters. These parameters can be estimated as follows:

#' # elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD2"),2,3)
#' # elog[, 'date'] <- as.Date(elog[, 'date'], format = '%Y%m%d')
#' # cal.cbs <- dc.ElogToCbsCbt(elog)$cal$cbs
#' # pnbd.EstimateParameters(cal.cbs, hardie = TRUE)
#' 
#' # (The final function was run several times with its own output as
#' # input for starting parameters, to ensure that the result converged).
#' 
#' params <- c(0.5629966, 12.5590370, 0.4081095, 10.5148048)
#' 
#' # 15% compounded annually has been converted to 0.0027 compounded continously,
#' # as we are dealing with weekly data and not annual data.
#' d <- 0.0027
#' 
#' pnbd.Plot.DERT(params = params, 
#'                x = 0:14, 
#'                t.x = 0:77, 
#'                T.cal = 77.86, 
#'                d = d, 
#'                hardie = TRUE, 
#'                type = "persp")
#' pnbd.Plot.DERT(params = params, 
#'                x = 0:14, 
#'                t.x = 0:77, 
#'                T.cal = 77.86, 
#'                d = d, 
#'                hardie = TRUE, 
#'                type="contour")
pnbd.Plot.DERT <- function(params, 
                           x, 
                           t.x, 
                           T.cal, 
                           d, 
                           hardie = TRUE, 
                           type = "persp") {
  # No use for inputs, other than as error check.
  inputs <- try(dc.InputCheck(params = params, 
                                func = "pnbd.Plot.DERT", 
                                printnames = c("r", "alpha", "s", "beta"), 
                                x = x, 
                                t.x = t.x, 
                                T.cal = T.cal))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  if (length(T.cal) > 1) stop("T.cal must be a single numeric value and may not be negative.")
  if (!(type %in% c("persp", "contour"))) {
    stop("The plot type in pnbd.Plot.DERT must be either 'persp' or 'contour'.")
  }
    
  DERT <- matrix(NA, length(t.x), length(x))
  rownames(DERT) <- t.x
  colnames(DERT) <- x
  for (i in 1:length(t.x)) {
    for (j in 1:length(x)) {
      DERT[i, j] <- pnbd.DERT(params, 
                              x[j], 
                              t.x[i], 
                              T.cal, 
                              d, 
                              hardie)
    }
  }
    
  if (type == "contour") {
    if (max(DERT, na.rm = TRUE) <= 10) {
      levels <- 1:max(DERT, na.rm = TRUE)
    } else if (max(DERT, na.rm = TRUE) <= 20) {
      levels <- c(1, seq(2, max(DERT, na.rm = TRUE), 2))
    } else {
      levels <- c(1, 2, seq(5, max(DERT, na.rm = TRUE), 5))
    }
    contour(x = t.x, 
            y = x, 
            z = DERT, 
            levels = levels, 
            xlab = "Recency", 
            ylab = "Frequency", 
            main = "Iso-Value Representation of DERT")
  }
  
  if (type == "persp") {
    persp(x = t.x, 
          y = x, 
          z = DERT, 
          theta = -30, 
          phi = 20, 
          axes = TRUE, 
          ticktype = "detailed", 
          nticks = 5, 
          main = "DERT as a Function of Frequency and Recency", 
          shade = 0.5, 
          xlab = "Recency", 
          ylab = "Frequency", 
          zlab = "Discounted expected residual transactions")
  }
  return(DERT)
}

#' Pareto/NBD Tracking Cumulative Transactions Plot
#'
#' Plots the actual and expected cumulative total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' actual.cu.tracking.data does not have to be in the same unit of time as the
#' T.cal data. T.tot will automatically be divided into periods to match the
#' length of actual.cu.tracking.data. See
#' [`pnbd.ExpectedCumulativeTransactions`].
#'
#' The holdout period should immediately follow the calibration period. This
#' function assume that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using max(T.cal) - T.cal rather than assuming that it is 0).
#'
#' @inheritParams pnbd.ExpectedCumulativeTransactions
#' @param actual.cu.tracking.data A vector containing the cumulative number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods). See details.
#' @param xlab Descriptive label for the x axis.
#' @param ylab Descriptive label for the y axis.
#' @param xticklab Vector containing a label for each tick mark on the x axis.
#' @param title Title placed on the top-center of the plot.
#' @return Matrix containing actual and expected cumulative repeat transactions.
#' 
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
#' # parameters estimated using pnbd.EstimateParameters
#' est.params <- cdnowSummary$est.params
#' 
#' # All parameters are in weeks; the calibration period lasted 39
#' # weeks and the holdout period another 39.
#' pnbd.PlotTrackingCum(params = est.params, 
#'                      T.cal = cal.cbs[,"T.cal"], 
#'                      T.tot = 78, 
#'                      actual.cu.tracking.data = cu.tracking)
#' @md
pnbd.PlotTrackingCum <- function(params, 
                                 T.cal, 
                                 T.tot, 
                                 actual.cu.tracking.data, 
                                 n.periods.final = NA,
                                 xlab = "Week", 
                                 ylab = "Cumulative Transactions", 
                                 xticklab = NULL, 
                                 title = "Tracking Cumulative Transactions") {
  # No use for inputs, other than as error check, so suppress
  # any warnings about incompatible vector lengths here:
  inputs <- suppressWarnings(try(dc.InputCheck(params = params, 
                                               func = "pnbd.PlotTrackingCum", 
                                               printnames = c("r", "alpha", "s", "beta"), 
                                               T.cal = T.cal, 
                                               T.tot = T.tot,
                                               actual.cu.tracking.data = actual.cu.tracking.data, 
                                               n.periods.final = n.periods.final)))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  inputs <- NULL
  if (length(T.tot) > 1) stop("T.tot must be a single numeric value and may not be negative.")
  
  actual <- actual.cu.tracking.data
  if(is.na(n.periods.final)) n.periods.final <- length(actual)
  expected <- pnbd.ExpectedCumulativeTransactions(params = params, 
                                                  T.cal = T.cal, 
                                                  T.tot = T.tot, 
                                                  n.periods.final = n.periods.final)
  
  cu.tracking.comparison <- rbind(actual, expected)
  
  ylim <- c(0, max(c(actual, expected)) * 1.05)
  plot(actual, 
       type = "l", 
       xaxt = "n", 
       xlab = xlab, 
       ylab = ylab, 
       col = 1, 
       ylim = ylim, 
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
  abline(v = max(T.cal), 
         lty = 2)
  
  legend("bottomright", 
         legend = c("Actual", "Model"), 
         col = 1:2, 
         lty = 1:2, 
         lwd = 1)
  
  return(cu.tracking.comparison)
}

#' Pareto/NBD Tracking Incremental Transactions Comparison
#'
#' Plots the actual and expected incremental total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' actual.inc.tracking.data does not have to be in the same unit of time as the
#' T.cal data. T.tot will automatically be divided into periods to match the
#' length of actual.inc.tracking.data. See
#' [`pnbd.ExpectedCumulativeTransactions`].
#'
#' The holdout period should immediately follow the calibration period. This
#' function assume that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using max(T.cal) - T.cal rather than assuming that it is 0).
#'
#' @inheritParams pnbd.PlotTrackingCum
#' @param actual.inc.tracking.data A vector containing the incremental number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods). See details.
#' @return Matrix containing actual and expected incremental repeat transactions.
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
#' # parameters estimated using pnbd.EstimateParameters
#' est.params <- cdnowSummary$est.params
#' 
#' # All parameters are in weeks; the calibration period lasted 39
#' # weeks and the holdout period another 39.
#' pnbd.PlotTrackingInc(params = est.params, 
#'                      T.cal = cal.cbs[,"T.cal"], 
#'                      T.tot = 78, 
#'                      actual.inc.tracking.data = inc.tracking)
#' @md
pnbd.PlotTrackingInc <- function(params, 
                                 T.cal, 
                                 T.tot, 
                                 actual.inc.tracking.data, 
                                 n.periods.final = NA,
                                 xlab = "Week", 
                                 ylab = "Transactions", 
                                 xticklab = NULL, 
                                 title = "Tracking Weekly Transactions") {
  # No use for inputs, other than as error check, so suppress
  # any warnings about incompatible vector lengths here:
  inputs <- suppressWarnings(try(dc.InputCheck(params = params, 
                                               func = "pnbd.PlotTrackingInc", 
                                               printnames = c("r", "alpha", "s", "beta"), 
                                               T.cal = T.cal, 
                                               T.tot = T.tot,
                                               actual.inc.tracking.data = actual.inc.tracking.data, 
                                               n.periods.final = n.periods.final)))
  if('try-error' == class(inputs)) return(str(inputs)$message)
  inputs <- NULL
  if (length(T.tot) > 1) stop("T.tot must be a single numeric value and may not be negative.")
  
  actual <- actual.inc.tracking.data
  if(is.na(n.periods.final)) n.periods.final <- length(actual)
  expected <- pnbd.ExpectedCumulativeTransactions(params = params, 
                                                  T.cal = T.cal, 
                                                  T.tot = T.tot, 
                                                  n.periods.final = n.periods.final)
  expected <- dc.CumulativeToIncremental(expected)

  ylim <- c(0, max(c(actual, expected)) * 1.05)
  plot(actual, 
       type = "l", 
       xaxt = "n", 
       xlab = xlab, 
       ylab = ylab, 
       col = 1, 
       ylim = ylim, 
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
  abline(v = max(T.cal), 
         lty = 2)
  
  legend("topright", 
         legend = c("Actual", "Model"), 
         col = 1:2, 
         lty = 1:2, 
         lwd = 1)
  return(rbind(actual, expected))
}

#' Plot Pareto/NBD Rate Heterogeneity
#'
#' A helper for plotting either the estimated gamma distribution of mu
#' (customers' propensities to drop out), or the estimated gamma distribution of
#' lambda (customers' propensities to purchase).
#'
#' @inheritParams pnbd.LL
#' @param func A string that is either "pnbd.PlotDropoutRateHeterogeneity" or
#'   "pnbd.PlotTransactionRateHeterogeneity".
#' @param lim The upper-bound of the x-axis. A number is chosen by the function
#'   if none is provided.
#' @return Depending on the value of `func`, either the distribution of
#'   customers' propensities to purchase or the distribution of customers'
#'   propensities to drop out.
#' @seealso [`pnbd.PlotDropoutRateHeterogeneity`]
#' @seealso [`pnbd.PlotTransactionRateHeterogeneity`]
#' @md
pnbd.PlotRateHeterogeneity <- function(params, 
                                       func, 
                                       lim = NULL) {
  stopifnot(func %in% c("pnbd.PlotDropoutRateHeterogeneity", 
                        "pnbd.PlotTransactionRateHeterogeneity"))
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, 
                        func = func)   
  shape_rate <- list(pnbd.PlotTransactionRateHeterogeneity = c(shape = params[1], 
                                                               rate = params[2]), 
                     pnbd.PlotDropoutRateHeterogeneity = c(shape = params[3], 
                                                           rate = params[4]))
  xlab_main <- list(pnbd.PlotTransactionRateHeterogeneity = c(xlab = "Transaction Rate", 
                                                              main = "Heterogeneity in Transaction Rate"), 
                    pnbd.PlotDropoutRateHeterogeneity = c(xlab = "Dropout Rate", 
                                                          main = "Heterogeneity in Dropout Rate"))
  shape <- shape_rate[[func]]['shape']
  rate <- shape_rate[[func]]['rate']
  rate.mean <- round(shape/rate, 4)
  rate.var <- round(shape/rate^2, 4)
  if (is.null(lim)) {
    lim = qgamma(0.99, shape = shape, rate = rate)
  }
  x.axis.ticks <- seq(0, lim, length.out = 100)
  heterogeneity <- dgamma(x.axis.ticks, 
                          shape = shape, 
                          rate = rate)
  plot(x.axis.ticks, 
       heterogeneity, 
       type = "l", 
       xlab = xlab_main[[func]]['xlab'], 
       ylab = "Density", 
       main = xlab_main[[func]]['main'])
  mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
  mtext(mean.var.label, side = 3)
  return(rbind(x.axis.ticks, heterogeneity))  
}

#' Pareto/NBD Plot Transaction Rate Heterogeneity
#'
#' Plots and returns the estimated gamma distribution of lambda (customers'
#' propensities to purchase).
#'
#' This returns the distribution of each customer's Poisson parameter, which
#' determines the level of their purchasing (using the Pareto/NBD assumption
#' that purchasing on the individual level can be modeled with a Poisson
#' distribution).
#'
#' @inheritParams pnbd.PlotRateHeterogeneity
#' @return Distribution of customers' propensities to purchase.
#' @examples
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' pnbd.PlotTransactionRateHeterogeneity(params)
#' params <- c(3, 10.56, 0.61, 11.64)
#' pnbd.PlotTransactionRateHeterogeneity(params)
pnbd.PlotTransactionRateHeterogeneity <- function(params, 
                                                  lim = NULL) {
  pnbd.PlotRateHeterogeneity(params = params, 
                             func = "pnbd.PlotTransactionRateHeterogeneity", 
                             lim = lim)
}

#' Pareto/NBD Plot Dropout Rate Heterogeneity
#'
#' Plots and returns the estimated gamma distribution of mu (customers'
#' propensities to drop out).
#'
#' This returns the distribution of each customer's exponential parameter that
#' determines their lifetime (using the Pareto/NBD assumption that a customer's
#' lifetime can be modeled with an exponential distribution).
#'
#' @inheritParams pnbd.PlotRateHeterogeneity
#' @return Distribution of customers' propensities to drop out.
#' @examples
#' params <- c(0.55, 10.56, 0.61, 11.64)
#' pnbd.PlotDropoutRateHeterogeneity(params)
#' params <- c(0.55, 10.56, 3, 11.64)
#' pnbd.PlotDropoutRateHeterogeneity(params)
pnbd.PlotDropoutRateHeterogeneity <- function(params, 
                                              lim = NULL) {
  pnbd.PlotRateHeterogeneity(params = params, 
                             func = "pnbd.PlotDropoutRateHeterogeneity", 
                             lim = lim)
}


 
