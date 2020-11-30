## Methods to model and forecast the amount that members are spending during
## transactions.

library(hypergeo)
library(lattice)

# Now trying Markdown + Roxygen (https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html)
#' Define general parameters
#'
#' This is to ensure consistency across all spend functions.
#'
#' This function is only ever called by functions defined in the original BTYD
#' package, such as [`spend.LL`], [`spend.marginal.likelihood`] or
#' [`spend.expected.value`] so it returns directly the output that is expected
#' from those calling functions.
#'
#' @inheritParams spend.LL
#' @param func name of the function calling [`dc.InputCheck`].
#' @return That depends on `func`: 1. If `func` is `spend.marginal.likelihood`,
#'   the marginal distribution of a customer's average transaction value (if m.x
#'   or x has a length greater than 1, a vector of marginal likelihoods will be
#'   returned). 2. If `func` is `spend.LL`, the log-likelihood of the
#'   gamma-gamma model; if m.x or x has a length greater than 1, this is a
#'   vector of log-likelihoods. 3. If `func` is `spend.expected.value`, the
#'   expected transaction value for a customer conditional on their transaction
#'   behavior during the calibration period. If m.x or x has a length greater
#'   than one, then a vector of expected transaction values will be returned.
#' @seealso [`spend.LL`]
#' @seealso [`spend.marginal.likelihood`]
#' @md
spend.generalParams  <- function(params, 
                                 func,
                                 m.x, 
                                 x) {
    stopifnot(func %in% c('spend.marginal.likelihood', 
                          'spend.LL', 
                          'spend.expected.value'))    
    inputs <- try(dc.InputCheck(params = params, 
                                printnames = c("p", "q", "gamma"),
                                func = func, 
                                m.x = m.x, 
                                x = x))
    if('try-error' == class(inputs)) return(inputs)
    x <- inputs$x
    m.x <- inputs$m.x
    max.length <- nrow(inputs)
    if (any(x == 0) || any(m.x == 0)) {
        warning("Customers with 0 transactions or 0 average spend in spend.marginal.likelihood")
    }
    
    p <- params[1]
    q <- params[2]
    gamma <- params[3]
    if(func == 'spend.expected.value') {
        M <- (gamma + m.x * x) * p/(p * x + q - 1)
        return(M)
    }
    
    result <- rep(0, max.length)
    
    ## non.zero: a vector indicating which elements have neither x == 0 or m.x == 0
    non.zero <- which(x > 0 & m.x > 0)
    common_piece <- q * log(gamma) + 
        (p * x[non.zero] - 1) * log(m.x[non.zero]) + 
        (p * x[non.zero]) * log(x[non.zero]) - 
        (p * x[non.zero] + q) * log(gamma + m.x[non.zero] * x[non.zero])
    if(func == 'spend.marginal.likelihood') {
        # return the marginal likelihood of a customer's average transaction value.
        result[non.zero] <- exp(lgamma(p * x[non.zero] + q) - 
                                    lgamma(p * x[non.zero]) - 
                                    lgamma(q) + 
                                    common_piece)
    } else if(func == 'spend.LL') {  
        result[non.zero] <- (-lbeta(p * x[non.zero], q) + common_piece)
    } else {
        return(NULL)
    }
    result
}

#' Gamma-gamma marginal likelihood
#'
#' Calculates the marginal likelihood of a customer's average transaction value.
#'
#' m.x and x may be vectors. The standard rules for vector operations apply - if
#' they are not of the same length, the shorter vector will be recycled (start
#' over at the first element) until it is as long as the longest vector. It is
#' advisable to keep vectors to the same length and to use single values for
#' parameters that are to be the same for all calculations. If one of these
#' parameters has a length greater than one, the output will be a vector of
#' probabilities.
#'
#' This function will issue a warning if any of m.x or x is 0, and will return a
#' marginal likelihood of 0 for those values.
#'
#' f(m.x | p, q, gamma, x).
#'
#' @param params	a vector of gamma-gamma parameters: p, q, and gamma, in that
#'   order. p is the shape parameter for each transaction. The scale parameter
#'   for each transaction is distributed across customers according to a gamma
#'   distribution with parameters q (shape) and gamma (scale).
#' @param m.x	the customer's average observed transaction value in the
#'   calibration period. May also be a vector of average observed transaction
#'   values - see details.
#' @param x	the number of transactions the customer made in the calibration
#'   period. May also be a vector of frequencies - see details.
#' @return The marginal distribution of a customer's average transaction value.
#'   If m.x or x has a length greater than 1, a vector of marginal likelihoods
#'   will be returned.
#' @references Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. “RFM and CLV:
#'   Using Iso-Value Curves for Customer Base Analysis.” Journal of Marketing
#'   Research Vol.42, pp.415-430. November. 2005.
#'   [Web.](http://www.brucehardie.com/papers/rfm_clv_2005-02-16.pdf)
#' @references See equation 3.
#' @examples
#' params <- c(6, 4, 16)
#' 
#' # calculate the marginal distribution of the average transaction value
#' # of a customer who spent an average of $35 over 3 transactions.
#' spend.marginal.likelihood(params, m.x=35, x=3)
#' 
#' # Several values can also be computed at once:
#' spend.marginal.likelihood(params, m.x=30:40, x=3)
#' spend.marginal.likelihood(params, m.x=35, x=1:10)
#' spend.marginal.likelihood(params, m.x=30:40, x=1:11)
#' @md
spend.marginal.likelihood <- function(params, 
                                      m.x, 
                                      x) {
    spend.generalParams(params = params,
                        func = 'spend.marginal.likelihood', 
                        m.x = m.x, 
                        x = x)
}

#' Spend Log-Likelihood
#'
#' Calculates the log-likelihood of the gamma-gamma model for customer spending.
#'
#' m.x and x may be vectors. The standard rules for vector operations apply - if
#' they are not of the same length, the shorter vector will be recycled (start
#' over at the first element) until it is as long as the longest vector. It is
#' advisable to keep vectors to the same length and to use single values for
#' parameters that are to be the same for all calculations. If one of these
#' parameters has a length greater than one, the output will be a vector of
#' log-likelihoods.
#'
#' @inheritParams spend.marginal.likelihood
#' @return The log-likelihood of the gamma-gamma model. If m.x or x has a length
#'   greater than 1, this is a vector of log-likelihoods.
#' @references Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. “RFM and CLV:
#'   Using Iso-Value Curves for Customer Base Analysis.” Journal of Marketing
#'   Research Vol.42, pp.415-430. November. 2005.
#'   [Web.](http://www.brucehardie.com/papers/rfm_clv_2005-02-16.pdf)
#' @examples \dontrun{
#' data(cdnowSummary)
#' ave.spend <- cdnowSummary$m.x
#' tot.trans <- cdnowSummary$cbs[,"x"]
#' # params <- c(6.25, 3.74, 15.44) # in original documentation. check below:
#' params <- spend.EstimateParameters(m.x.vector = ave.spend, x.vector = tot.trans)
#' # get the total log-likelihood of the data and parameters
#' # above. There will be many warnings due to the zeroes that are
#' # included in the data. If you wish to avoid these warnings, use:
#' 
#' # ave.spend <- ave.spend[which(tot.trans > 0)]
#' # tot.trans <- tot.trans[which(tot.trans > 0)]
#' 
#' # Note that we used tot.trans to remove the zeroes from ave.spend.
#' # This is because we need the vectors to be the same length, and it
#' # is possible that your data include customers who made transactions
#' # worth zero dollars (in which case the vector lengths would differ
#' # if we used ave.spend to remove the zeroes from ave.spend).
#' 
#' sum(spend.LL(params, ave.spend, tot.trans))
#' 
#' # This log-likelihood may be different than mentioned in the
#' # referenced paper; in the paper, a slightly different function
#' # which relies on total spend (not average spend) is used.
#' }
#' @md
spend.LL <- function(params, 
                     m.x, 
                     x) {
    spend.generalParams(params = params,
                        func = 'spend.LL', 
                        m.x = m.x, 
                        x = x)
}

#' Conditional expected transaction value
#'
#' Calculates the expected transaction value for a customer, conditional on the
#' number of transaction and average transaction value during the calibration
#' period.
#'
#' E(M | p, q, gamma, m.x, x).
#'
#' m.x and x may be vectors. The standard rules for vector operations apply - if
#' they are not of the same length, the shorter vector will be recycled (start
#' over at the first element) until it is as long as the longest vector. It is
#' advisable to keep vectors to the same length and to use single values for
#' parameters that are to be the same for all calculations. If one of these
#' parameters has a length greater than one, the output will be a vector of
#' probabilities.
#'
#' @inheritParams spend.marginal.likelihood
#' @return The expected transaction value for a customer conditional on their
#'   transaction behavior during the calibration period. If m.x or x has a
#'   length greater than one, then a vector of expected transaction values will
#'   be returned.
#' @references Fader, Peter S., Bruce G.S. Hardie, and Ka L. Lee. “RFM and CLV:
#'   Using Iso-Value Curves for Customer Base Analysis.” Journal of Marketing
#'   Research Vol.42, pp.415-430. November. 2005.
#'   [Web.](http://www.brucehardie.com/papers/rfm_clv_2005-02-16.pdf)
#' @examples \dontrun{
#' data(cdnowSummary)
#' ave.spend <- cdnowSummary$m.x
#' tot.trans <- cdnowSummary$cbs[,"x"]
#' # params <- c(6, 4, 16); # in original documentation. rounded values of:
#' params <- spend.EstimateParameters(m.x.vector = ave.spend, x.vector = tot.trans);
#' # calculate the expected transaction value of a customer
#' # who spent an average of $35 over 3 transactions.
#' spend.expected.value(params, m.x=35, x=3)
#' 
#' # m.x and x may be vectors:
#' spend.expected.value(params, m.x=30:40, x=3)
#' spend.expected.value(params, m.x=35, x=1:10)
#' spend.expected.value(params, m.x=30:40, x=1:11)
#' }
#' @md
spend.expected.value <- function(params, 
                                 m.x, 
                                 x) {
    spend.generalParams(params = params,
                        func = 'spend.expected.value', 
                        m.x = m.x, 
                        x = x)    
}

#' Spend Parameter Estimation
#'
#' Estimates parameters for the gamma-gamma spend model.
#'
#' The best-fitting parameters are determined using the spend.LL function. The
#' sum of the log-likelihood for each customer (for a set of parameters) is
#' maximized in order to estimate parameters.
#'
#' A set of starting parameters must be provided for this method. If no
#' parameters are provided, (1,1,1,1) is used as a default. It may be necessary
#' to run the estimation from multiple starting points to ensure that it
#' converges. To compare the log-likelihoods of different parameters, use
#' \link{spend.LL}.
#'
#' The lower bound on the parameters to be estimated is always zero, since
#' gamma-gamma parameters cannot be negative. The upper bound can be set with
#' the max.param.value parameter.
#'
#' @param m.x.vector  a vector with each customer's average observed transaction
#'   value in the calibration period.
#' @param x.vector  a vector with the number of transactions each customer made
#'   in the calibration period. Must correspond to m.x.vector in terms of
#'   ordering of customers and length of the vector.
#' @param par.start  initial vector of gamma-gamma parameters: p, q, and gamma,
#'   in that order. p is the shape parameter for each transaction. The scale
#'   parameter for each transaction is distributed across customers according to
#'   a gamma distribution with parameters q (shape) and gamma (scale).
#' @param max.param.value  the upper bound on parameters.
#' @return Vector of estimated parameters.
#' @examples \dontrun{
#' data(cdnowSummary)
#' ave.spend <- cdnowSummary$m.x
#' tot.trans <- cdnowSummary$cbs[,"x"]
#' 
#' # There will be many warnings due to the zeroes that are
#' # included in the data above. To avoid them, use the following:
#' # (see example for spend.LL)
#' 
#' ave.spend <- ave.spend[which(tot.trans > 0)]
#' tot.trans <- tot.trans[which(tot.trans > 0)]
#' 
#' # We will let the spend function use default starting parameters
#' spend.EstimateParameters(ave.spend, tot.trans)
#' }
spend.EstimateParameters <- function(m.x.vector, 
                                     x.vector, 
                                     par.start = c(1, 1, 1), 
                                     max.param.value = 10000) {
    
    if (any(m.x.vector < 0) || !is.numeric(m.x.vector)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x.vector < 0) || !is.numeric(x.vector)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (length(m.x.vector) != length(x.vector))
        stop("m.x.vector and x.vector must be the same length.")
    if (any(x.vector == 0) || any(m.x.vector == 0)) 
        warning("Customers with 0 transactions or 0 average spend in spend.LL")
    
    spend.eLL <- function(params, 
                          m.x.vector, 
                          x.vector, 
                          max.param.value) {
        params <- exp(params)
        params[params > max.param.value] <- max.param.value
        return(-1 * sum(spend.LL(params = params, 
                                 m.x = m.x.vector, 
                                 x = x.vector)))
    }
    logparams <- log(par.start)
    results <- optim(logparams, 
                     spend.eLL, 
                     m.x.vector = m.x.vector, 
                     x.vector = x.vector, 
                     max.param.value = max.param.value, 
                     method = "L-BFGS-B")
    estimated.params <- exp(results$par)
    estimated.params[estimated.params > max.param.value] <- max.param.value
    return(estimated.params)
}

#' Plot Actual vs. Expected Average Transaction Value
#'
#' Plots the actual and expected densities of average transaction values, and
#' returns a vector with each customer's average transaction value probability.
#'
#' @inheritParams spend.EstimateParameters
#' @param params	a vector of gamma-gamma parameters: p, q, and gamma, in that
#'   order. p is the shape parameter for each transaction. The scale parameter
#'   for each transaction is distributed across customers according to a gamma
#'   distribution with parameters q (shape) and gamma (scale).
#' @param xlab  descriptive label for the x axis.
#' @param ylab  descriptive label for the y axis.
#' @param title  title placed on the top-center of the plot.
#' @return a vector with the probability of each customer's average transaction
#'   value.
#' @seealso [`spend.marginal.likelihood`]
#' @examples \dontrun{
#' data(cdnowSummary)
#' ave.spend <- cdnowSummary$m.x
#' tot.trans <- cdnowSummary$cbs[,"x"]
#' # params <- c(6.25, 3.74, 15.44) # in original documentation. check below:
#' params <- spend.EstimateParameters(m.x.vector = ave.spend, x.vector = tot.trans)
#' 
#' # Plot the actual and expected average transaction value across customers.
#' f.m.x <- spend.plot.average.transaction.value(params, ave.spend, tot.trans)
#' }
#' @md
spend.plot.average.transaction.value <- function(params, 
                                                 m.x.vector, 
                                                 x.vector, 
                                                 xlab = "Average Transaction Value", 
                                                 ylab = "Marginal Distribution of Average Transaction Value", 
                                                 title = "Actual vs. Expected Average Transaction Value Across Customers") {
    
    if (any(m.x.vector < 0) || !is.numeric(m.x.vector)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x.vector < 0) || !is.numeric(x.vector)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (length(m.x.vector) != length(x.vector))
        stop("m.x.vector and x.vector must be the same length.")
    if (any(x.vector == 0) || any(m.x.vector == 0)) {
        warning(paste("There are customers with 0 transactions or 0 average spend.", 
                      "spend.plot.average.transaction.value removed them before plotting."))
    }
    
    # remove any customers with zero repeat transactions
    ave.spending <- m.x.vector[which(x.vector > 0)]
    tot.transactions <- x.vector[which(x.vector > 0)]
    
    f.m.x <- spend.marginal.likelihood(params, 
                                       ave.spending, 
                                       tot.transactions)
    plot(ave.spending, 
         y = f.m.x, 
         pch = 16, 
         type = "n", 
         xlab = xlab, 
         ylab = ylab, 
         main = title)
    lines(density(ave.spending, 
                  bw = "nrd", 
                  adjust = 0.6), 
          col = 1, 
          lty = 1)
    lines(smooth.spline(ave.spending, 
                        y = f.m.x, 
                        w = tot.transactions, 
                        df = 15), 
          col = 2, 
          lty = 2)
    legend("topright", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lty = 1:2, 
           lwd = 1)
    return(f.m.x)
}

