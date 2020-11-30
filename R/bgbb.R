################################################################################ Beta-Geometric Beta-Binomial Functions

library(hypergeo)

#' BG/BB Log-Likelihood using a recency-frequency matrix
#'
#' Calculates the log-likelihood of the BG/BB model.
#'
#' @param params BG/BB parameters - a vector with alpha, beta, gamma, and delta,
#'   in that order. Alpha and beta are unobserved parameters for the
#'   beta-Bernoulli transaction process. Gamma and delta are unobserved
#'   parameters for the beta-geometric dropout process.
#' @param rf.matrix recency-frequency matrix. It must contain columns for
#'   frequency ("x"), recency ("t.x"), number of transaction opportunities in
#'   the calibration period ("n.cal"), and the number of customers with this
#'   combination of recency, frequency and transaction opportunities in the
#'   calibration period ("custs"). Note that recency must be the time between
#'   the start of the calibration period and the customer's last transaction,
#'   not the time between the customer's last transaction and the end of the
#'   calibration period.
#' @seealso [`bgbb.LL`]
#' @return The total log-likelihood of the provided data in rf.matrix.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' bgbb.rf.matrix.LL(params, rf.matrix)
#' @md
bgbb.rf.matrix.LL <- function(params, 
                              rf.matrix) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.rf.matrix.LL")
    
    tryCatch(x.vec <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x.vec <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.periods.vec <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(n.custs.vec <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    LL.sum <- sum(n.custs.vec * 
                      bgbb.LL(params, 
                              x.vec, 
                              t.x.vec, 
                              n.periods.vec))
    
    return(LL.sum)
}

#' BG/BB Log-Likelihood
#'
#' Calculates the log-likelihood of the BG/BB model.
#'
#' x, t.x, and n.cal may be vectors. The standard rules for vector operations
#' apply - if they are not of the same length, shorter vectors will be recycled
#' (start over at the first element) until they are as long as the longest
#' vector. It is advisable to keep vectors to the same length and to use single
#' values for parameters that are to be the same for all calculations. If one of
#' these parameters has a length greater than one, the output will be also be a
#' vector.
#'
#' @param params    BG/BB parameters - a vector with alpha, beta, gamma, and
#'   delta, in that order. Alpha and beta are unobserved parameters for the
#'   beta-Bernoulli transaction process. Gamma and delta are unobserved
#'   parameters for the beta-geometric dropout process.
#' @param x     the number of repeat transactions made by the customer in the
#'   calibration period. Can also be vector of frequencies - see details.
#' @param t.x   recency - the transaction opportunity in which the customer made
#'   their last transaction. Can also be a vector of recencies - see details.
#' @param n.cal     number of transaction opportunities in the calibration
#'   period. Can also be a vector of calibration period transaction
#'   opportunities - see details.
#' @return A vector of log-likelihoods as long as the longest input vector (x,
#'   t.x, or n.cal).
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' 
#' # Returns the log likelihood of the parameters for a customer who
#' # made 3 transactions in a calibration period with 6 transaction opportunities,
#' # with the last transaction occurring during the 4th transaction opportunity.
#' bgbb.LL(params, x=3, t.x=4, n.cal=6)
#' 
#' # We can also give vectors as function parameters:
#' set.seed(7)
#' x <- sample(1:3, 10, replace = TRUE)
#' t.x <- sample(3:5, 10, replace = TRUE)
#' n.cal <- rep(5, 10)
#' bgbb.LL(params, x, t.x, n.cal)  
#' @md
bgbb.LL <- function(params, 
                    x, 
                    t.x, 
                    n.cal) {
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.LL', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs)
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    max.length <- nrow(inputs)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    denom.ab <- lbeta(alpha, beta)
    denom.gd <- lbeta(gamma, delta)
    
    indiv.LL.sum <- lbeta(alpha + x, 
                          beta + n.cal - x) - 
        denom.ab + 
        lbeta(gamma, delta + n.cal) - 
        denom.gd
    
    check <- n.cal - t.x - 1
    addition <- function(alpha, 
                         beta, 
                         gamma, 
                         delta, 
                         denom.ab, 
                         denom.gd, 
                         x, 
                         t.x, 
                         check) {
        ii <- 0:check
        # implement log-sum-exp trick as shown on Wikipedia:
        # https://en.wikipedia.org/wiki/LogSumExp
        xset <- lbeta(alpha + x, 
                      beta + t.x - x + ii) - 
            denom.ab + 
            lbeta(gamma + 1, 
                  delta + t.x + ii) - 
            denom.gd
        # was:
        # log(sum(exp(xset)))
        # now:
        xstar <- max(xset)
        xdiff <- xset - xstar
        xstar + log(sum(exp(xdiff)))
    }
    # for every element of vectors for which t.x<n.cal, add the result of 'addition'
    # in logspace.  addLogs defined in dc.R. addLogs(a,b) = log(exp(a) + exp(b))
    for (i in 1:max.length) {
        if (check[i] >= 0) 
            indiv.LL.sum[i] <- addLogs(indiv.LL.sum[i], 
                                       addition(alpha, 
                                                beta, 
                                                gamma, 
                                                delta, 
                                                denom.ab, 
                                                denom.gd, 
                                                x[i], 
                                                t.x[i], 
                                                check[i]))
    }
    return(indiv.LL.sum)
}

#' BG/BB Parameter estimation
#'
#' Estimates parameters for the BG/BB model.
#'
#' The best-fitting parameters are determined using the [`bgbb.rf.matrix.LL`]
#' function. The sum of the log-likelihood for each customer (for a set of
#' parameters) is maximized in order to estimate paramaters.
#'
#' A set of starting parameters must be provided for this method. If no
#' parameters are provided, (1,1,1,1) is used as a default. It may be useful to
#' use starting values for parameters that represent your best guess of the
#' heterogeneity in the transaction and dropout rates of customers. It may be
#' necessary to run the estimation from multiple starting points to ensure that
#' it converges. To compare the log-likelihoods of different parameters, use
#' [`bgbb.rf.matrix.LL`].
#'
#' The lower bound on the parameters to be estimated is always zero, since BG/BB
#' parameters cannot be negative. The upper bound can be set with the
#' max.param.value parameter.
#'
#' @inheritParams bgbb.rf.matrix.LL
#' @param par.start initial BG/BB parameters - a vector with alpha, beta, gamma,
#'   and delta, in that order. Alpha and beta are unobserved parameters for the
#'   beta-Bernoulli transaction process. Gamma and delta are unobserved
#'   parameters for the beta-geometric dropout process.
#' @param max.param.value the upper bound on parameters.
#' @return Vector of estimated paramaters.
#' @seealso [`bgbb.rf.matrix.LL`]
#' @examples
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' # log-likelihood of estimated parameters
#' bgbb.rf.matrix.LL(est.params, rf.matrix)
#' @md
bgbb.EstimateParameters <- function(rf.matrix, 
                                    par.start = c(1, 1, 1, 1), 
                                    max.param.value = 1000) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = par.start, 
                          func = "bgbb.EstimateParameters")
    
    bgbb.eLL <- function(params, 
                         rf.matrix, 
                         max.param.value) {
        params <- exp(params)
        params[params > max.param.value] <- max.param.value
        return(-1 * bgbb.rf.matrix.LL(params, rf.matrix))
    }
    logparams <- log(par.start)
    results <- optim(logparams, 
                     bgbb.eLL, 
                     rf.matrix = rf.matrix, 
                     max.param.value = max.param.value, 
                     method = "L-BFGS-B")
    estimated.params <- exp(results$par)
    estimated.params[estimated.params > max.param.value] <- max.param.value
    return(estimated.params)
}

#' BG/BB Probability Mass Function
#'
#' Probability mass function for the BG/BB.
#'
#' P(X(n)=x | alpha, beta, gamma, delta). Returns the probability that a
#' customer makes x transactions in the first n transaction opportunities.
#'
#' Parameters `n` and `x` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, the shorter vector
#' will be recycled (start over at the first element) until it is as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be a
#' vector of probabilities.
#' 
#' @inheritParams bgbb.rf.matrix.LL
#' @param n     number of transaction opportunities; may also be a vector.
#' @param x 	number of transactions; may also be a vector.
#' @return Probability of X(n)=x, conditional on model parameters.
#' @seealso [`bgbb.pmf.General`]
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples 
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' # The probability that a customer made 3 transactions in the first
#' # 6 transaction opportunities.
#' bgbb.pmf(params, n=6, x=3)
#' 
#' # Vectors may also be used as arguments:
#' bgbb.pmf(params, n=6, x=0:6)
#' @md
bgbb.pmf <- function(params, 
                     n, 
                     x) {
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.pmf', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                n = n, 
                                x = x))
    if('try-error' == class(inputs)) return(inputs)
    if (any(inputs$x > inputs$n)) {
        stop("bgbb.pmf was given x > n")
    }
    
    return(bgbb.pmf.General(params = params, 
                            n.cal = 0, 
                            n.star = inputs$n, 
                            x.star = inputs$x))
}

#' BG/BB General Probability Mass Function 
#' 
#' Calculates the probability that a customer will make `x.star` transactions in
#' the first `n.star` transaction opportunities following the calibration
#' period.
#' 
#' P(X(n, n + n*) = x* | alpha, beta, gamma, delta). This is a more generalized
#' version of the bgbb.pmf. Setting `n.cal` to 0 reduces this function to the
#' probability mass function in its usual format - the probability that a user
#' will make x.star transactions in the first n.star transaction opportunities.
#'
#' It is impossible for a customer to make a negative number of transactions, or
#' to make more transactions than there are transaction opportunities. This
#' function will throw an error if such inputs are provided.
#'
#' `n.cal`, `n.star`, and `x.star` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be a
#' vector of probabilities.
#' 
#' @inheritParams bgbb.LL
#' @param n.star number of transaction opportunities in the holdout period, or a
#'   vector of holdout period transaction opportunities.
#' @param x.star number of transactions in the holdout period, or a vector of
#'   transaction frequencies.
#' @return Probability of X(n, n + n*) = x*, given BG/BB model parameters.
#' @seealso [`bgbb.pmf`]
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples 
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' # Probability that a customer will make 3 transactions in the 10
#' # transaction opportunities following the 6 transaction opportunities
#' # in the calibration period, given BG/BB parameters.
#' bgbb.pmf.General(params, n.cal=6, n.star=10, x.star=3)
#' 
#' # Vectors may also be provided as input:
#' # Comparison between different frequencies:
#' bgbb.pmf.General(params, n.cal=6, n.star=10, x.star=1:10)
#' # Comparison between different holdout transaction opportunities:
#' bgbb.pmf.General(params, n.cal=6, n.star=5:15, x.star=3)
#' @md
bgbb.pmf.General <- function(params, 
                             n.cal, 
                             n.star, 
                             x.star) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.pmf.General', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                n.cal = n.cal, 
                                n.star = n.star,
                                x.star = x.star))
    if('try-error' == class(inputs)) return(inputs) 
    
    n.cal <- inputs$n.cal
    n.star <- inputs$n.star
    x.star <- inputs$x.star
    max.length <- nrow(inputs)
    
    if (any(x.star > n.star)) {
        stop("bgbb.pmf.General was given x.star > n.star")
    }
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- rep(0, max.length)
    piece.1[x.star == 0] <- 1 - exp(lbeta(gamma, 
                                          delta + n.cal[x.star == 0]) - 
                                        lbeta(gamma, 
                                              delta))
    
    piece.2 <- exp(lchoose(n.star, x.star) + 
                       lbeta(alpha + x.star, 
                             beta + n.star - x.star) - 
                       lbeta(alpha, beta) + 
                       lbeta(gamma, delta + n.cal + n.star) - 
                       lbeta(gamma, delta))
    
    piece.3 = rep(0, max.length)
    rows.to.sum <- which(x.star <= n.star - 1)
    piece.3[rows.to.sum] <- unlist(sapply(rows.to.sum, 
                                          function(index) {
                                              ii <- x.star[index]:(n.star[index] - 1)
                                              sum(exp(lchoose(ii, x.star[index]) + 
                                                          lbeta(alpha + x.star[index], 
                                                                beta + ii - x.star[index]) - 
                                                          lbeta(alpha, beta) + 
                                                          lbeta(gamma + 1, 
                                                                delta + n.cal[index] + ii) - 
                                                          lbeta(gamma, delta)))
                                          }))
    
    expectation <- piece.1 + piece.2 + piece.3
    return(expectation)
}

#' BG/BB Expectation
#' 
#' Returns the number of transactions that a randomly chosen customer (for whom
#' we have no prior information) is expected to make in the first n transaction
#' opportunities.
#' 
#' E(X(n) | alpha, beta, gamma, delta)
#' 
#' @inheritParams bgbb.pmf
#' @return Mean of the BG/BB probability mass function.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples 
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' # Expected number of transactions that a randomly chosen customer
#' # will make in the first 10 transaction opportunities.
#' bgbb.Expectation(params, n=10)
#' 
#' # We can also compare expected transactions over time:
#' bgbb.Expectation(params, n=1:10)
#' @md
bgbb.Expectation <- function(params, 
                             n) {
    
    # we don't need inputs here but the parameter checks
    # we need for this function can be done with dc.InputCheck()    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.Expectation', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                n = n))
    if('try-error' == class(inputs)) return(inputs)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.one <- (alpha/(alpha + beta)) * (delta/(gamma - 1))
    piece.two <- exp(lgamma(gamma + delta) - 
                         lgamma(gamma + delta + n) + 
                         lgamma(1 + delta + n) - 
                         lgamma(1 + delta))
    return(piece.one * (1 - piece.two))
}

#' BG/BB P(Alive)
#' 
#' Uses BG/BB model parameters and a customer's past transaction behavior to
#' return the probability that they will be alive in the transaction opportunity
#' following the calibration period.
#'
#' `x`, `t.x`, and `n.cal` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be a
#' vector of probabilities.
#' 
#' P(alive at n+1 | alpha, beta, gamma, delta, x, t.x, n)
#' 
#' @inheritParams bgbb.LL
#' @return Probability that the customer is alive at the (n+1)th transaction
#'   opportunity. If `x`, `t.x`, and/or `n.cal` are of length greater than one,
#'   then this will be a vector of probabilities (containing one element
#'   matching each element of the longest input vector).
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' 
#' # The probability that a customer who made 3 transactions in
#' # the calibration period (which consisted of 6 transaction
#' # opportunities), with the last transaction occurring at the
#' # 4th transaction opportunity, is alive at the 7th transaction
#' # opportunity
#' bgbb.PAlive(params, x=3, t.x=4, n.cal=6)
#' 
#' # The input parameters may also be vectors:
#' bgbb.PAlive(params, x=1, t.x=1:6, n.cal=6)
#' @md
bgbb.PAlive <- function(params, 
                        x, 
                        t.x, 
                        n.cal) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PAlive', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + x, beta + n.cal - x) - 
                       lbeta(alpha, beta) + 
                       lbeta(gamma, delta + n.cal + 1) - 
                       lbeta(gamma, delta))
    piece.2 <- 1/exp(bgbb.LL(params, 
                             x, 
                             t.x, 
                             n.cal))
    p.alive <- piece.1 * piece.2
    
    return(p.alive)
}

#' BG/BB Discounted Expected Residual Transactions
#'
#' Computes the number of discounted expected residual transactions by a
#' customer, conditional on their behavior in the calibration period.
#'
#' DERT(d | alpha, beta, gamma, delta, x, t.x, n). This is the present value of
#' the expected future transaction stream for a customer with x transactions and
#' a recency of t.x in n.cal transaction opportunities, discounted by a rate d.
#'
#' `x`, `t.x`, and `n.cal` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be
#' also be a vector.
#'
#' @inheritParams bgbb.PAlive
#' @param d discount rate.
#' @return The present value of the expected future transaction stream for a particular customer.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' See equation 14.
#' @examples 
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' # Compute DERT for a customer who made 3 transactions
#' # in the calibration period(consisting of 6 transaction
#' # opportunities), with the last transaction occurring
#' # during the 4th transaction opportunity, discounted at
#' # 10%.
#' bgbb.DERT(params, x=3, t.x=4, n.cal=6, d=0.1)
#' 
#' # We can also compare DERT for several customers:
#' bgbb.DERT(params, x=1:6, t.x=6, n.cal=6, d=0.1)
#' @md
bgbb.DERT <- function(params, 
                      x, 
                      t.x, 
                      n.cal, 
                      d) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.DERT', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + x + 1, 
                         beta + n.cal - x) - 
                       lbeta(alpha, 
                             beta))
    piece.2 <- exp(lbeta(gamma, 
                         delta + n.cal + 1) - 
                       lbeta(gamma, 
                             delta))/(1 + d)
    piece.3 <- Re(hypergeo(A = 1, 
                           B = delta + n.cal + 1, 
                           C = gamma + delta + n.cal + 1, 
                           z = 1/(1 + d)))
    piece.4 <- exp(bgbb.LL(params, x, t.x, n.cal))
    
    dert <- piece.1 * piece.2 * (piece.3/piece.4)
    
    return(dert)
}

#' BG/BB Discounted Expected Residual Transactions using a recency-frequency matrix
#' 
#' Computes the number of discounted expected residual transactions by a
#' customer, conditional on their behavior in the calibration period.
#' 
#' @inheritParams bgbb.DERT
#' @inheritParams bgbb.rf.matrix.LL
#' @return The present value of the expected future transaction stream for a particular customer.
#' @seealso [`bgbb.DERT`]
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' See equation 14.
#' @examples 
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # compute DERT for a customer from every row in rf.matrix,
#' # discounted at 10%.
#' bgbb.rf.matrix.DERT(est.params, rf.matrix, d = 0.1)
#' @md
bgbb.rf.matrix.DERT <- function(params, 
                                rf.matrix, 
                                d) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.rf.matrix.DERT")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.DERT(params, 
                     x, 
                     t.x, 
                     n.cal, 
                     d))
}

#' BG/BB Conditional Expected Transactions
#'
#' Calculates the number of expected transactions in the holdout period,
#' conditional on a customer's behavior in the calibration period.
#'
#' E(X(n, n+n*) | alpha, beta, gamma, delta, x, t.x, n). This function requires
#' the holdout period to immediately follow the calibration period.
#'
#' `n.cal`, `n.star`, `x`, and `t.x` may be vectors. The standard rules for
#' vector operations apply - if they are not of the same length, shorter vectors
#' will be recycled (start over at the first element) until they are as long as
#' the longest vector. It is advisable to keep vectors to the same length and to
#' use single values for parameters that are to be the same for all
#' calculations. If one of these parameters has a length greater than one, the
#' output will be a vector of probabilities.
#'
#' @inheritParams bgbb.LL
#' @inheritParams bgbb.pmf.General
#' @return The number of transactions a customer is expected to make in the
#'   `n.star` transaction opportunities following the calibration period,
#'   conditional on their behavior during the calibration period.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples  
#' params <- c(1.20, 0.75, 0.66, 2.78)
#' # the number of transactions a customer is expected
#' # to make in the 10 transaction opportunities
#' # following the calibration period, which consisted
#' # of 6 transaction opportunities (during which they
#' # made 3 transactions, the last of which occurred
#' # in the 4th opportunity)
#' bgbb.ConditionalExpectedTransactions(params, n.cal=6, n.star=10, x=3, t.x=4)
#' 
#' # We can also use vectors as input:
#' bgbb.ConditionalExpectedTransactions(params, n.cal=6, n.star=1:10, x=3, t.x=4)
#' bgbb.ConditionalExpectedTransactions(params, n.cal=6, n.star=10, x=1:4, t.x=4)
#' @md
bgbb.ConditionalExpectedTransactions <- function(params, 
                                                 n.cal, 
                                                 n.star, 
                                                 x, 
                                                 t.x) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.ConditionalExpectedTransactions', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                n.cal = n.cal, 
                                n.star = n.star,
                                x = x, 
                                t.x = t.x))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    n.star <- inputs$n.star
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- 1/exp(bgbb.LL(params, x, t.x, n.cal))
    piece.2 <- exp(lbeta(alpha + x + 1, 
                         beta + n.cal - x) - 
                       lbeta(alpha, 
                             beta))
    piece.3 <- delta/(gamma - 1)
    piece.4 <- exp(lgamma(gamma + delta) - 
                       lgamma(1 + delta))
    piece.5 <- exp(lgamma(1 + delta + n.cal) - 
                       lgamma(gamma + delta + n.cal))
    piece.6 <- exp(lgamma(1 + delta + n.cal + n.star) - 
                       lgamma(gamma + delta + n.cal + n.star))
    
    expected.frequency <- piece.1 * piece.2 * piece.3 * piece.4 * (piece.5 - piece.6)
    
    which.is.nan <- is.nan(expected.frequency)
    if (sum(which.is.nan) > 0) {
        error.msg.long <- paste("numerical error, parameters exploded in bgbb.ConditionalExpectedTransactions", 
                                "params:", alpha, beta, gamma, delta, "n.cal:", n.cal[which.is.nan], 
                                "n.star:", n.star[which.is.nan], "x:", x[which.is.nan], "t.x:", t.x[which.is.nan], 
                                "...", "piece.1:", piece.1[which.is.nan], "piece.2:", piece.2[which.is.nan], 
                                "piece.3:", piece.3[which.is.nan], "piece.4:", piece.4[which.is.nan], 
                                "piece.5:", piece.5[which.is.nan], "piece.6:", piece.6[which.is.nan])
        stop(error.msg.long)
    }
    return(expected.frequency)
}

#' BG/BB Plot Frequency in Calibration Period
#'
#' Plots the actual and expected number of customers who made a certain number
#' of repeat transactions in the calibration period. Also returns a matrix with
#' this comparison.
#'
#' @inheritParams bgbb.rf.matrix.LL
#' @param censor optional. Any calibration period frequency at this number, or
#'   above it, will be binned together. If the censor number is greater than the
#'   maximum recency in the recency-frequency matrix, the maximum recency will
#'   be used as the censor number.
#' @param plotZero If FALSE, the histogram will exclude the zero bin.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param title title placed on the top-center of the plot.
#' @return Calibration period repeat transaction frequency comparison matrix,
#'   actual vs. expected.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples 
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # plot actual vs. expected frequencies in the calibration period
#' bgbb.PlotFrequencyInCalibration(est.params, rf.matrix)
#' @md
bgbb.PlotFrequencyInCalibration <- function(params, 
                                            rf.matrix, 
                                            censor = NULL, 
                                            plotZero = TRUE, 
                                            xlab = "Calibration period transactions", 
                                            ylab = "Customers", 
                                            title = "Frequency of Repeat Transactions") {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.PlotFrequencyInCalibration")
    
    tryCatch(x <- rf.matrix[, "x"], 
             error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], 
             error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], 
             error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    max.n.cal <- max(n.cal)
    if (is.null(censor)) 
        censor <- max.n.cal
    total.custs <- sum(custs)
    actual.frequency <- rep(0, max.n.cal + 1)
    expected.frequency <- rep(0, max.n.cal + 1)
    
    for (ii in 0:max.n.cal) {
        actual.frequency[ii + 1] <- sum(custs[x == ii])
        expected.frequency[ii + 1] <- sum(unlist(sapply(unique(n.cal[n.cal >= ii]), 
                                                        function(this.n.cal) {
                                                            sum(custs[n.cal == this.n.cal]) * bgbb.pmf(params, this.n.cal, ii)
                                                        })))
    }
    
    freq.comparison <- rbind(actual.frequency, expected.frequency)
    colnames(freq.comparison) <- 0:max.n.cal
    
    if (ncol(freq.comparison) <= censor) {
        censored.freq.comparison <- freq.comparison
    } else {
        ## Rename for easier coding
        fc <- freq.comparison
        ## Build censored freq comparison (cfc)
        cfc <- fc
        cfc <- cfc[, 1:(censor + 1)]
        cfc[1, (censor + 1)] <- sum(fc[1, (censor + 1):ncol(fc)])
        cfc[2, (censor + 1)] <- sum(fc[2, (censor + 1):ncol(fc)])
        
        censored.freq.comparison <- cfc
    }
    
    if (plotZero == FALSE) 
        censored.freq.comparison <- censored.freq.comparison[, -1]
    
    n.ticks <- ncol(censored.freq.comparison)
    if (plotZero == TRUE) {
        x.labels <- 0:(n.ticks - 1)
        if (censor < ncol(freq.comparison) - 1) 
            x.labels[n.ticks] <- paste(n.ticks - 1, "+", sep = "")
    } else {
        x.labels <- 1:(n.ticks)
        if (censor < ncol(freq.comparison) - 1) 
            x.labels[n.ticks] <- paste(n.ticks, "+", sep = "")
    }
    colnames(censored.freq.comparison) <- x.labels
    
    ylim <- c(0, 
              ceiling(max(c(censored.freq.comparison[1, ], 
                            censored.freq.comparison[2, ])) * 1.1))
    
    barplot(censored.freq.comparison, 
            beside = TRUE, 
            ylim = ylim, 
            main = title, 
            xlab = xlab, 
            ylab = ylab, 
            col = 1:2)
    legend("top", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lwd = 2, 
           cex = 0.75)
    
    return(censored.freq.comparison)
}

#' BG/BB Plot Frequency in Holdout
#'
#' Plots the actual and expected number of customers who made a certain number
#' of transactions in the holdout period, binned according to holdout period
#' frequencies. Also returns a matrix with this comparison and the number of
#' customers in each bin.
#'
#' @inheritParams bgbb.PlotFrequencyInCalibration
#' @param n.cal number of transaction opportunities in the calibration period.
#' @param rf.matrix.holdout holdout period recency-frequency matrix. It must
#'   contain columns for frequency in the holdout period ("x.star"), the number
#'   of transaction opportunities in the holdout period ("n.star"), and the
#'   number of customers with each frequency ("custs").
#' @return Holdout period repeat transaction frequency comparison matrix (actual vs. expected).
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples 
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' rf.matrix.holdout <- donationsSummary$rf.matrix.holdout
#' # donationsSummary$rf.matrix and donationsSummary$rf.matrix.holdout already
#' # have appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # number of periods in the calibration period
#' n.cal = max(rf.matrix[,"n.cal"])
#' 
#' bgbb.PlotFrequencyInHoldout(est.params, n.cal, rf.matrix.holdout)
#' @md
bgbb.PlotFrequencyInHoldout <- function(params, 
                                        n.cal, 
                                        rf.matrix.holdout, 
                                        censor = NULL, 
                                        plotZero = TRUE, 
                                        title = "Frequency of Repeat Transactions", 
                                        xlab = "Holdout period transactions", 
                                        ylab = "Customers") {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.PlotFrequencyInHoldout")
    if (n.cal < 0 || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not be negative.")
    
    tryCatch(x.star <- rf.matrix.holdout[, "x.star"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a frequency column labelled \"x.star\""))
    tryCatch(n.star <- rf.matrix.holdout[, "n.star"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a column with the number of transaction opportunities for each group, labelled \"n.star\""))
    tryCatch(custs <- rf.matrix.holdout[, "custs"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    max.n.star <- max(n.star)
    if (is.null(censor)) 
        censor <- max.n.star
    total.custs <- sum(custs)
    actual.frequency <- rep(0, max.n.star + 1)
    expected.frequency <- rep(0, max.n.star + 1)
    
    for (ii in 0:max.n.star) {
        actual.frequency[ii + 1] <- sum(custs[x.star == ii])
        expected.frequency[ii + 1] <- sum(unlist(sapply(unique(n.star[n.star >= ii]), 
                                                        function(this.n.star) {
                                                            sum(custs[n.star == this.n.star]) * bgbb.pmf.General(params, n.cal, 
                                                                                                                 this.n.star, ii)
                                                        })))
    }
    
    freq.comparison <- rbind(actual.frequency, expected.frequency)
    colnames(freq.comparison) <- 0:max.n.star
    
    if (ncol(freq.comparison) <= censor) {
        censored.freq.comparison <- freq.comparison
    } else {
        ## Rename for easier coding
        fc <- freq.comparison
        ## Build censored freq comparison (cfc)
        cfc <- fc
        cfc <- cfc[, 1:(censor + 1)]
        cfc[1, (censor + 1)] <- sum(fc[1, (censor + 1):ncol(fc)])
        cfc[2, (censor + 1)] <- sum(fc[2, (censor + 1):ncol(fc)])
        
        if (plotZero == FALSE) {
            cfc <- cfc[, -1]
        }
        
        censored.freq.comparison <- cfc
    }
    
    if (plotZero == TRUE) {
        x.labels <- 0:(ncol(censored.freq.comparison) - 1)
    } else {
        x.labels <- 1:(ncol(censored.freq.comparison))
    }
    
    if (censor < ncol(freq.comparison) - 1) {
        x.labels[(censor + 1)] <- paste(censor, "+", sep = "")
    }
    colnames(censored.freq.comparison) <- x.labels
    
    barplot(censored.freq.comparison, 
            beside = TRUE, 
            main = title, 
            xlab = xlab, 
            ylab = ylab, 
            col = 1:2)
    legend("topright", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lwd = 2)
    
    return(censored.freq.comparison)
}

#' BG/BB Tracking Cumulative Transactions Plot
#'
#' Plots the actual and expected cumulative total repeat transactions by all
#' customers for the calibration and holdout periods. Also returns a matrix with
#' this comparison.
#'
#' The holdout period should immediately follow the calibration period. This
#' function assumes that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using `max(n.cal) - n.cal` rather than assuming that they are all 0).
#'
#' @inheritParams bgbb.rf.matrix.LL
#' @inheritParams bgbb.PlotFrequencyInHoldout
#' @param actual.cum.repeat.transactions vector containing the cumulative number
#'   of repeat transactions made by customers in all transaction opportunities
#'   (both calibration and holdout periods). Its unit of time should be the same
#'   as the units of the recency-frequency matrix used to estimate the model
#'   parameters.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @return Matrix containing actual and expected cumulative repeat transactions.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples 
#' data(donationsSummary)
#' # donationsSummary$rf.matrix already has appropriate column names
#' rf.matrix <- donationsSummary$rf.matrix
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # get the annual repeat transactions, and transform them into
#' # a cumulative form
#' actual.inc.repeat.transactions <- donationsSummary$annual.trans
#' actual.cum.repeat.transactions <- cumsum(actual.inc.repeat.transactions)
#' 
#' # set appropriate x-axis
#' x.tickmarks <- c( "'96","'97","'98","'99","'00","'01","'02","'03","'04","'05","'06" )
#' 
#' # plot actual vs. expected transactions. The calibration period was 6 periods long.
#' bgbb.PlotTrackingCum(est.params, rf.matrix, actual.cum.repeat.transactions, xticklab=x.tickmarks)
#' @md
bgbb.PlotTrackingCum <- function(params, 
                                 rf.matrix, 
                                 actual.cum.repeat.transactions, 
                                 xlab = "Time", 
                                 ylab = "Cumulative Transactions", 
                                 xticklab = NULL, 
                                 title = "Tracking Cumulative Transactions") {
    
    # we don't need inputs here but the parameter checks
    # we need for this plot can be done with dc.InputCheck()
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PlotTrackingCum', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                actual.cum.repeat.transactions = actual.cum.repeat.transactions))
    if('try-error' == class(inputs)) return(inputs)
    
    tryCatch(n.cal <- rf.matrix[, "n.cal"], 
             error = function(e) stop("Error in bgbb.PlotTrackingCum: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], 
             error = function(e) stop("Error in bgbb.PlotTrackingCum: rf.matrix must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    actual <- actual.cum.repeat.transactions
    n.periods <- length(actual)
    
    cust.birth.periods <- max(n.cal) - n.cal
    
    expected <- sapply(1:n.periods, function(interval) {
        if (interval <= min(cust.birth.periods)) 
            return(0)
        sum(bgbb.Expectation(params, 
                             interval - 
                                 cust.birth.periods[cust.birth.periods <= interval]) * 
                custs[cust.birth.periods <= interval])
    })
    
    pur.comparison <- rbind(actual, expected)
    
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
        if (ncol(pur.comparison) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size in bgbb.PlotTrackingCum")
        }
        axis(1, at = 1:ncol(pur.comparison), labels = xticklab)
    }
    if (is.null(n.cal) == FALSE) {
        abline(v = max(n.cal), lty = 2)
    }
    legend("bottomright", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lty = 1:2, 
           lwd = 1, 
           cex = 0.75)
    
    return(pur.comparison)
}

#' BG/BB Tracking Incremental Transactions Plot
#'
#' Plots the actual and expected incremental total repeat transactions by all
#' customers for the calibration and holdout periods. Also returns a matrix of
#' this comparison.
#'
#' The holdout period should immediately follow the calibration period. This
#' function assumes that all customers' calibration periods end on the same date,
#' rather than starting on the same date (thus customers' birth periods are
#' determined using `max(n.cal) - n.cal` rather than assuming that they are all 0).
#'
#' @inheritParams bgbb.PlotTrackingCum
#' @param actual.inc.repeat.transactions vector containing the incremental
#'   number of repeat transactions made by customers in all transaction
#'   opportunities (both calibration and holdout periods). Its unit of time
#'   should be the same as the units of the recency-frequency matrix used to
#'   estimate the model parameters.
#' @return Matrix containing actual and expected incremental repeat transactions.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples 
#' data(donationsSummary)
#' # donationsSummary$rf.matrix already has appropriate column names
#' rf.matrix <- donationsSummary$rf.matrix
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # get the annual repeat transactions
#' actual.inc.repeat.transactions <- donationsSummary$annual.trans
#' 
#' # Set appropriate x-axis
#' x.tickmarks <- c( "'96","'97","'98","'99","'00","'01","'02","'03","'04","'05","'06" )
#' 
#' # Plot actual vs. expected transactions. The calibration period was 6 periods long.
#' bgbb.PlotTrackingInc(est.params, rf.matrix, actual.inc.repeat.transactions, xticklab=x.tickmarks)
#' @md
bgbb.PlotTrackingInc <- function(params, 
                                 rf.matrix, 
                                 actual.inc.repeat.transactions, 
                                 xlab = "Time", 
                                 ylab = "Transactions", 
                                 xticklab = NULL, 
                                 title = "Tracking Incremental Transactions") {
    
    # we don't need inputs here but the parameter checks
    # we need for this plot can be done with dc.InputCheck()
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PlotTrackingInc', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                actual.inc.repeat.transactions = actual.inc.repeat.transactions))
    if('try-error' == class(inputs)) return(inputs)    
    
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotTrackingInc: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotTrackingInc: rf.matrix must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    actual <- actual.inc.repeat.transactions
    n.periods <- length(actual)
    
    cust.birth.periods <- max(n.cal) - n.cal
    
    expected.cumulative <- sapply(1:n.periods, function(interval) {
        if (interval <= min(cust.birth.periods)) 
            return(0)
        sum(bgbb.Expectation(params, interval - cust.birth.periods[cust.birth.periods <= 
                                                                       interval]) * custs[cust.birth.periods <= interval])
    })
    
    expected <- dc.CumulativeToIncremental(expected.cumulative)
    
    pur.comparison <- rbind(actual, expected)
    
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
        if (ncol(pur.comparison) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size in bgbb.PlotTrackingInc")
        }
        axis(1, at = 1:ncol(pur.comparison), labels = xticklab)
    }
    if (is.null(n.cal) == FALSE) {
        abline(v = max(n.cal), lty = 2)
    }
    legend("topright", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lty = 1:2, 
           lwd = 1, 
           cex = 0.75)
    
    return(pur.comparison)
}

#' BG/BB Plot Frequency vs Conditional Expected Frequency
#'
#' Plots the actual and conditional expected number of transactions made by
#' customers in the holdout period, binned according to calibration period
#' frequencies. Also returns a matrix with this comparison and the number of
#' customers in each bin.
#'
#' @inheritParams bgbb.PlotTrackingCum
#' @param n.star number of transaction opportunities in the holdout period.
#' @param x.star a vector containing the number of transactions made in the
#'   holdout period by the groups of customers with the same recency and
#'   frequency in the calibration period. It must be in the same order as the
#'   rf.matrix.
#' @param trunc optional integer used to truncate the plot. In the plot, all
#'   calibration period frequencies above the truncation number will be removed.
#'   If the truncation number is greater than the maximum frequency, R will warn
#'   you and change it to the maximum frequency.
#' @return Holdout period transaction frequency comparison matrix (actual vs.
#'   expected), binned by calibration period frequency.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples
#' data(donationsSummary)
#'
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#'
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#'
#' # get the holdout period transactions
#' x.star <- donationsSummary$x.star
#'
#' # number of transaction opportunities in the holdout period
#' n.star <- 5
#'
#' # Plot holdout period transactions
#' bgbb.PlotFreqVsConditionalExpectedFrequency(est.params, n.star, rf.matrix, x.star, trunc=6)
#' @md
bgbb.PlotFreqVsConditionalExpectedFrequency <- function(params, 
                                                        n.star, 
                                                        rf.matrix, 
                                                        x.star, 
                                                        trunc = NULL, 
                                                        xlab = "Calibration period transactions", 
                                                        ylab = "Holdout period transactions", 
                                                        xticklab = NULL, 
                                                        title = "Conditional Expectation") {
    
    if (length(x.star) != nrow(rf.matrix)) 
        stop("x.star must have the same number of entries as rows in rf.matrix")
    if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
        stop("n.star must be a single value or have as many entries as rows in rf.matrix")
    
    dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                          params = params, 
                          func = "bgbb.PlotFreqVsConditionalExpectedFrequency")
    
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    
    n.star <- rep(n.star, length.out = nrow(rf.matrix))
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    if (is.null(trunc)) 
        trunc <- max(n.cal)
    
    if (trunc > max(n.cal)) {
        warning("The truncation number provided in bgbb.PlotFreqVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
                max(n.cal))
        trunc = max(n.cal)
    }
    
    actual.freq <- rep(0, max(n.cal) + 1)
    expected.freq <- rep(0, max(n.cal) + 1)
    bin.size <- rep(0, max(n.cal) + 1)
    
    for (ii in 0:max(n.cal)) {
        bin.size[ii + 1] <- sum(custs[x == ii])
        actual.freq[ii + 1] <- sum(x.star[x == ii])
        expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
                                                                          n.cal[x == ii], n.star[x == ii], ii, t.x[x == ii]) * custs[x == ii])
    }
    
    comparison <- rbind(actual.freq/bin.size, expected.freq/bin.size, bin.size)
    colnames(comparison) <- paste("freq.", 0:max(n.cal), sep = "")
    
    if (is.null(xticklab) == FALSE) {
        x.labels <- xticklab
    } else {
        if ((trunc + 1) < ncol(comparison)) {
            x.labels <- 0:(trunc)
        } else {
            x.labels <- 0:(ncol(comparison) - 1)
        }
    }
    
    actual.freq <- comparison[1, 1:(trunc + 1)]
    expected.freq <- comparison[2, 1:(trunc + 1)]
    
    custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
    if (custs.in.plot < 0.9 * sum(custs)) {
        warning("Less than 90% of customers are represented in your plot (", custs.in.plot, 
                " of ", sum(custs), " are plotted).")
    }
    
    ylim <- c(0, ceiling(max(c(actual.freq, expected.freq)) * 1.1))
    plot(actual.freq, 
         type = "l", 
         xaxt = "n", 
         col = 1, 
         ylim = ylim, 
         xlab = xlab, 
         ylab = ylab, 
         main = title)
    lines(expected.freq, 
          lty = 2, 
          col = 2)
    
    axis(1, 
         at = 1:(trunc + 1), 
         labels = x.labels)
    legend("topleft", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lty = 1:2, 
           lwd = 1)
    
    return(comparison)
}

#' BG/BB Plot Recency vs Conditional Expected Frequency
#'
#' Plots the actual and conditional expected number of transactions made by
#' customers in the holdout period, binned according to calibration period
#' recencies. Also returns a matrix with this comparison and the number of
#' customers in each bin.
#'
#' @inheritParams bgbb.PlotFreqVsConditionalExpectedFrequency
#' @return Holdout period transaction frequency comparison matrix (actual vs.
#'   expected), binned by calibration period recency.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @examples
#' data(donationsSummary)
#'
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#'
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#'
#' # get the holdout period transactions
#' x.star <- donationsSummary$x.star
#'
#' # number of transaction opportunities in the holdout period
#' n.star <- 5
#'
#' # Compare holdout period transactions.
#' bgbb.PlotRecVsConditionalExpectedFrequency(est.params, n.star, rf.matrix, x.star, trunc=6)
#' @md
bgbb.PlotRecVsConditionalExpectedFrequency <- function(params, 
                                                       n.star, 
                                                       rf.matrix, 
                                                       x.star, 
                                                       trunc = NULL, 
                                                       xlab = "Calibration period recency", 
                                                       ylab = "Holdout period transactions", 
                                                       xticklab = NULL, 
                                                       title = "Conditional Expected Transactions by Recency") {
    
    if (length(x.star) != nrow(rf.matrix)) 
        stop("x.star must have the same number of entries as rows in rf.matrix")
    if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
        stop("n.star must be a single value or have as many entries as rows in rf.matrix")
    
    dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                          params = params, 
                          func = "bgbb.PlotRecVsConditionalExpectedFrequency")
    
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    
    n.star <- rep(n.star, length.out = nrow(rf.matrix))
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    if (is.null(trunc)) 
        trunc <- max(n.cal)
    
    if (trunc > max(n.cal)) {
        warning("The truncation number provided in bgbb.PlotRecVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
                max(n.cal))
        trunc = max(n.cal)
    }
    
    actual.freq <- rep(0, max(n.cal) + 1)
    expected.freq <- rep(0, max(n.cal) + 1)
    bin.size <- rep(0, max(n.cal) + 1)
    
    for (ii in 0:max(n.cal)) {
        bin.size[ii + 1] <- sum(custs[t.x == ii])
        actual.freq[ii + 1] <- sum(x.star[t.x == ii])
        expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
                                                                          n.cal[t.x == ii], n.star[t.x == ii], x[t.x == ii], ii) * custs[t.x == 
                                                                                                                                             ii])
    }
    
    comparison <- rbind(actual.freq/bin.size, 
                        expected.freq/bin.size, 
                        bin.size)
    colnames(comparison) <- paste("rec.", 0:max(n.cal), sep = "")
    
    custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
    if (custs.in.plot < 0.9 * sum(custs)) {
        warning("Less than 90% of customers are represented in your plot (", custs.in.plot, 
                " of ", sum(custs), " are plotted).")
    }
    
    if (is.null(xticklab) == FALSE) {
        x.labels <- xticklab
    } else {
        if ((trunc + 1) < ncol(comparison)) {
            x.labels <- 0:(trunc)
        } else {
            x.labels <- 0:(ncol(comparison) - 1)
        }
    }
    
    actual.freq <- comparison[1, 1:(trunc + 1)]
    expected.freq <- comparison[2, 1:(trunc + 1)]
    
    ylim <- c(0, ceiling(max(c(actual.freq, expected.freq)) * 1.1))
    plot(actual.freq, 
         type = "l", 
         xaxt = "n", 
         col = 1, 
         ylim = ylim, 
         xlab = xlab, 
         ylab = ylab, 
         main = title)
    lines(expected.freq, lty = 2, col = 2)
    
    axis(1, 
         at = 1:(trunc + 1), 
         labels = x.labels)
    legend("topleft", 
           legend = c("Actual", "Model"), 
           col = 1:2, 
           lty = 1:2, 
           lwd = 1)
    
    return(comparison)
}

#' BG/BB Posterior Mean (l,m)th Product Moment
#'
#' Computes the `(l,m)`th product moment of the joint posterior distribution of
#' P (the Bernoulli transaction process parameter) and Theta (the geometric
#' dropout process parameter).
#'
#' E((P)^l(Theta)^m | alpha, beta, gamma, delta, x, t.x, n)
#'
#' `x`, `t.x`, and `n.cal` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be
#' also be a vector.
#' 
#' @inheritParams bgbb.LL
#' @param l moment degree of P
#' @param m moment degree of Theta
#' @return The expected posterior `(l,m)`th product moment.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#'   
#' See equation 17.
#' @md
bgbb.PosteriorMeanLmProductMoment <- function(params, 
                                              l, 
                                              m, 
                                              x, 
                                              t.x, 
                                              n.cal) {
    if (l < 0 || length(l) != 1 || !is.numeric(l)) 
        stop("l must be a single numeric value and may not be less than 0.")
    if (m < 0 || length(m) != 1 || !is.numeric(m)) 
        stop("m must be a single numeric value and may not be less than 0.")
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PosteriorMeanLmProductMoment', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + l, beta) - 
                       lbeta(alpha, beta) + 
                       lbeta(gamma + m, delta) - 
                       lbeta(gamma, delta))
    piece.2 <- exp(bgbb.LL(c(alpha + l, beta, gamma + m, delta), 
                           x, 
                           t.x, 
                           n.cal))
    piece.3 <- exp(bgbb.LL(params, 
                           x, 
                           t.x, 
                           n.cal))
    
    mean <- piece.1 * (piece.2/piece.3)
    
    return(mean)
}

#' BG/BB Posterior Mean Transaction Rate
#'
#' Computes the mean value of the marginal posterior value of P, the Bernoulli
#' transaction process parameter.
#'
#' E(P | alpha, beta, gamma, delta, x, t.x, n). This is calculated by setting `l
#' = 1` and `m = 0` in [`bgbb.PosteriorMeanLmProductMoment`].
#'
#' `x`, `t.x`, and `n.cal` may be vectors. The standard rules for vector
#' operations apply - if they are not of the same length, shorter vectors will
#' be recycled (start over at the first element) until they are as long as the
#' longest vector. It is advisable to keep vectors to the same length and to use
#' single values for parameters that are to be the same for all calculations. If
#' one of these parameters has a length greater than one, the output will be
#' also be a vector.
#' 
#' @inheritParams bgbb.LL
#' @return The posterior mean transaction rate.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @seealso [`bgbb.rf.matrix.PosteriorMeanTransactionRate`]  
#' @examples  
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # return the posterior mean transaction rate vector
#' bgbb.rf.matrix.PosteriorMeanTransactionRate(est.params, rf.matrix)
#' @md  
bgbb.PosteriorMeanTransactionRate <- function(params, 
                                              x, 
                                              t.x, 
                                              n.cal) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PosteriorMeanTransactionRate', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    
    mean.transaction.rate <- bgbb.PosteriorMeanLmProductMoment(params = params, 
                                                               l = 1, 
                                                               m = 0, 
                                                               x = x, 
                                                               t.x = t.x, 
                                                               n.cal = n.cal)
    return(mean.transaction.rate)
}

#' BG/BB Posterior Mean Transaction Rate using a recency-frequency matrix
#'
#' Computes the mean value of the marginal posterior value of P, the Bernoulli
#' transaction process parameter.
#'
#' `rf.matrix` has columns x`, `t.x`, and `n.cal`. 
#' 
#' @inheritParams bgbb.rf.matrix.LL
#' @return The posterior mean transaction rate.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @seealso [`bgbb.PosteriorMeanTransactionRate`]  
#' @examples  
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # return the posterior mean transaction rate vector
#' bgbb.rf.matrix.PosteriorMeanTransactionRate(est.params, rf.matrix)
#' @md  
bgbb.rf.matrix.PosteriorMeanTransactionRate <- function(params, 
                                                        rf.matrix) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.rf.matrix.PosteriorMeanTransactionRate")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.PosteriorMeanTransactionRate(params, 
                                             x, 
                                             t.x, 
                                             n.cal))
}

#' BG/BB Posterior Mean Dropout Rate
#' 
#' Computes the mean value of the marginal posterior value of Theta, the geometric dropout process parameter.
#' 
#' E(Theta | alpha, beta, gamma, delta, x, t.x, n). This is calculated by setting `l = 0` and `m = 1` in [`bgbb.PosteriorMeanLmProductMoment`].
#' 
#' `x`, `t.x`, and `n.cal` may be vectors. The standard rules for vector operations apply - if they are not of the same length, shorter vectors will be recycled (start over at the first element) until they are as long as the longest vector. It is advisable to keep vectors to the same length and to use single values for parameters that are to be the same for all calculations. If one of these parameters has a length greater than one, the output will be also be a vector.
#' 
#' @inheritParams bgbb.LL
#' @return The posterior mean dropout rate.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @seealso [`bgbb.rf.matrix.PosteriorMeanDropoutRate`]  
#' @examples  
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # return the posterior mean dropout rate vector
#' bgbb.rf.matrix.PosteriorMeanDropoutRate(est.params, rf.matrix)
#' @md
bgbb.PosteriorMeanDropoutRate <- function(params, 
                                          x, 
                                          t.x, 
                                          n.cal) {
    
    inputs <- try(dc.InputCheck(params = params, 
                                func = 'bgbb.PosteriorMeanDropoutRate', 
                                printnames = c("alpha", "beta", "gamma", "delta"),
                                x = x, 
                                t.x = t.x, 
                                n.cal = n.cal))
    if('try-error' == class(inputs)) return(inputs) 
    
    x <- inputs$x
    t.x <- inputs$t.x
    n.cal <- inputs$n.cal
    
    mean.dropout.rate <- bgbb.PosteriorMeanLmProductMoment(params = params, 
                                                           l = 0, 
                                                           m = 1, 
                                                           x = x, 
                                                           t.x = t.x, 
                                                           n.cal = n.cal)
    return(mean.dropout.rate)
}

#' BG/BB Posterior Mean Dropout Rate using a recency-frequency matrix
#' 
#' Computes the mean value of the marginal posterior value of Theta, the geometric dropout process parameter.
#' 
#' E(Theta | alpha, beta, gamma, delta, x, t.x, n). This is calculated by setting `l = 0` and `m = 1` in [`bgbb.PosteriorMeanLmProductMoment`].
#' 
#' `rf.matrix` has columns x`, `t.x`, and `n.cal`. 
#' 
#' @inheritParams bgbb.rf.matrix.LL
#' @return The posterior mean dropout rate.
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/)
#' @seealso [`bgbb.PosteriorMeanDropoutRate`]  
#' @examples 
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # return the posterior mean dropout rate vector
#' bgbb.rf.matrix.PosteriorMeanDropoutRate(est.params, rf.matrix)
#' @md
bgbb.rf.matrix.PosteriorMeanDropoutRate <- function(params, 
                                                    rf.matrix) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.rf.matrix.PosteriorMeanDropoutRate")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.PosteriorMeanDropoutRate(params, 
                                         x, 
                                         t.x, 
                                         n.cal))
}

#' BG/BB Heatmap of Holdout Period Expected Transactions
#'
#' Plots a heatmap based on the conditional expected holdout period frequency
#' for each recency-frequency combination in the calibration period.
#'
#' E(X(n, n+n*) | alpha, beta, gamma, delta, x, t.x, n). This function requires
#' the holdout period to immediately follow the calibration period.
#'
#' @inheritParams bgbb.pmf
#' @param n.cal number of transaction opportunities in the calibration period.
#' @param n.star number of transaction opportunities in the holdout period.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @param title title placed on the top-center of the plot.
#' @return A matrix containing the conditional expected transactions in the
#'   holdout period for each recency-frequency combination in the calibration
#'   period. The rows represent calibration period frequencies, and the columns
#'   represent calibration period recencies.
#' @seealso [`bgbb.ConditionalExpectedTransactions`]
#' @references  Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
#'   "Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
#'   _Marketing Science_ 29(6), pp. 1086-1108. 2010. INFORMS.
#'   [Web.](http://www.brucehardie.com/papers/020/) 
#' @examples 
#' data(donationsSummary)
#' 
#' rf.matrix <- donationsSummary$rf.matrix
#' # donationsSummary$rf.matrix already has appropriate column names
#' 
#' # starting-point parameters
#' startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)
#' 
#' # Plot a heatmap of conditional expected transactions in
#' # a holdout period of 5 transaction opportunities, given
#' # that the calibration period consisted of 6 transaction
#' # opportunities.
#' bgbb.HeatmapHoldoutExpectedTrans(est.params, n.cal=6, n.star=5)
#' @md
bgbb.HeatmapHoldoutExpectedTrans <- function(params, 
                                             n.cal, 
                                             n.star, 
                                             xlab = "Recency", 
                                             ylab = "Frequency", 
                                             xticklab = NULL, 
                                             title = "Heatmap of Conditional Expected Transactions") {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.HeatmapHoldoutExpectedTrans")
    if (n.cal < 0 || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not be negative.")
    if (n.star < 0 || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not be negative.")
    
    heatmap.mx <- matrix(0, n.cal + 1, n.cal + 1)
    heatmap.mx[1, 1] <- bgbb.ConditionalExpectedTransactions(params, 
                                                             n.cal, 
                                                             n.star, 
                                                             0, 0)
    for (xx in 1:n.cal) {
        for (tt in 1:n.cal) {
            if (xx <= tt) {
                expected.trans <- bgbb.ConditionalExpectedTransactions(params, 
                                                                       n.cal, 
                                                                       n.star, 
                                                                       xx, 
                                                                       tt)
                heatmap.mx[xx + 1, tt + 1] <- expected.trans
            }
        }
    }
    if (is.null(xticklab) == TRUE) {
        xticklab <- 0:n.cal
    }
    colnames(heatmap.mx) <- xticklab
    rownames(heatmap.mx) <- 0:n.cal
    heatmap(heatmap.mx, 
            Rowv = NA, 
            Colv = NA, 
            col = gray(8:2/9), 
            scale = "none", 
            ylab = ylab, 
            xlab = xlab, 
            main = title, 
            verbose = TRUE)
    return(heatmap.mx)
}

#' BG/BB Plot Transaction Rate Heterogeneity
#' 
#' Plots and returns the estimated beta distribution of P (customers' propensities to purchase).
#' 
#' This returns the distribution of each customer's Bernoulli parameter, which determines the level of their purchasing (using the BG/BB assumption that purchasing on the individual level can be modeled with a Bernoulli distribution).
#' 
#' @inheritParams bgbb.pmf
#' @return Distribution of customers' propensities to purchase.
#' @examples 
#' params <- c(1.2, 0.75, 0.66, 2.78)
#' bgbb.PlotTransactionRateHeterogeneity(params)
#' params <- c(0.2, 1.5, 3.2, 6)
#' bgbb.PlotTransactionRateHeterogeneity(params)
#' @md
bgbb.PlotTransactionRateHeterogeneity <- function(params) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.PlotTransactionRateHeterogeneity")
    
    alpha <- params[1]
    beta <- params[2]
    x.axis.ticks <- 0.01 * 0:100
    heterogeneity <- dbeta(x.axis.ticks, alpha, beta)
    plot(x.axis.ticks, 
         heterogeneity, 
         type = "l", 
         xlab = "Transaction Rate", 
         ylab = "Density", 
         main = "Heterogeneity in Transaction Rate")
    rate.mean <- round(alpha/(alpha + beta), 4)
    rate.var <- round((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)), 4)
    mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
    mtext(mean.var.label, side = 3)
    return(heterogeneity)
}

#' BG/BB Plot Dropout Rate Heterogeneity
#' 
#' Plots and returns the estimated beta distribution of Theta (customers' propensities to drop out).
#' 
#' This returns the distribution of each customer's geometric parameter that determines their lifetime (using the BG/BB assumption that a customer's lifetime can be modeled with a geometric distribution).
#' 
#' @inheritParams bgbb.pmf
#' @return Distribution of customers' propensities to drop out.
#' @examples
#' params <- c(1.2, 0.75, 0.66, 2.78)
#' bgbb.PlotDropoutRateHeterogeneity(params)
#' params <- c(0.2, 1.5, 3.2, 6)
#' bgbb.PlotDropoutRateHeterogeneity(params)
#' @md
bgbb.PlotDropoutRateHeterogeneity <- function(params) {
    
    dc.check.model.params(printnames = c("alpha", "beta", "gamma", "delta"), 
                          params = params, 
                          func = "bgbb.PlotDropoutRateHeterogeneity")
    
    alpha <- params[3]
    beta <- params[4]
    x.axis.ticks <- 0.01 * 0:100
    heterogeneity <- dbeta(x.axis.ticks, alpha, beta)
    plot(x.axis.ticks, 
         heterogeneity, 
         type = "l", 
         xlab = "Dropout rate", 
         ylab = "Density", 
         main = "Heterogeneity in Dropout Rate")
    rate.mean <- round(alpha/(alpha + beta), 4)
    rate.var <- round((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)), 4)
    mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
    mtext(mean.var.label, side = 3)
    return(heterogeneity)
} 
