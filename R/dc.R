################################################################################ Functions for Manipulating Data

library(Matrix)

#' Check the inputs to functions that use this common pattern
#'
#' A bunch of functions whose names start with \code{pnbd} take a set of four
#' parameters as their first argument, and then a set of vectors or scalars such
#' as \code{x} or \code{T.cal} as their subsequent arguments. This function
#' started out as pnbd.InputCheck() and it was meant to run input checks for any
#' number of such subsequent vector arguments, as long as they all met the same
#' requirements as \code{x}, \code{t.x} and \code{T.cal} in
#' \code{\link{pnbd.LL}}: meaning, the length of the longest of these vectors is
#' a multiple of the lengths of all others, and all vectors are numeric and
#' positive.
#'
#' With an extra argument, \code{printnames}, pnbd.InputCheck() could also
#' accommodate input checks for functions whose names start with \code{bgbb},
#' \code{bgnbd}, and \code{spend} so it was basically useful everywhere. That's
#' when it became \code{dc.InputCheck()}. \code{params} can have any length as
#' long as that length is the same as the length of \code{printnames}, so
#' \code{dc.InputCheck()} can probably handle mixtures of distributions for
#' modeling BTYD behavior that are not yet implemented.
#'
#' By other arguments ... here we mean a bunch of named vectors that are used by
#' functions that call \code{dc.InputCheck}, such as x, t.x, T.cal, etc. The
#' standard rules for vector operations apply - if they are not of the same
#' length, shorter vectors will be recycled (start over at the first element)
#' until they are as long as the longest vector. Vector recycling is a good way
#' to get into trouble. Keep vectors to the same length and use single values
#' for parameters that are to be the same for all calculations. If one of these
#' parameters has a length greater than one, the output will be a vector of
#' probabilities.
#'
#' @param params If used by \code{pnbd.[...]} functions, Pareto/NBD parameters
#'   -- a vector with r, alpha, s, and beta, in that order. See
#'   \code{\link{pnbd.LL}}. If used by \code{bgnbd.[...]} functions, BG/NBD
#'   parameters -- a vector with r, alpha, a, and b, in that order. See
#'   \code{\link{bgnbd.LL}}. If used by \code{bgbb.[...]} functions, BG/BB
#'   parameters -- a vector with alpha, beta, gamma, and delta, in that order.
#'   See \code{\link{bgbb.LL}}. If used by \code{spend.[...]} functions, a
#'   vector of gamma-gamma parameters -- p, q, and gamma, in that order. See
#'   \code{\link{spend.LL}}.
#' @param func Function calling dc.InputCheck
#' @param printnames a string vector with the names of parameters to pass to
#'   \code{\link{dc.check.model.params}}
#' @param ... other arguments
#' @return If all is well, a data frame with everything you need in it, with
#'   nrow() equal to the length of the longest vector in \code{...}
#' @seealso \code{\link{pnbd.LL}}
#'   \code{\link{pnbd.ConditionalExpectedTransactions}}
dc.InputCheck <- function(params, 
                          func,
                          printnames = c("r", "alpha", "s", "beta"), 
                          ...) {
    inputs <- as.list(environment())
    vectors <- list(...)
    dc.check.model.params(printnames = inputs$printnames, 
                          params = inputs$params, 
                          func = inputs$func)
    max.length <- max(sapply(vectors, length))
    lapply(names(vectors), function(x) {
        if(max.length %% length(vectors[[x]])) 
            warning(paste("Maximum vector length not a multiple of the length of", 
                          x, sep = " "))
        if (any(vectors[[x]] < 0) || !is.numeric(vectors[[x]])) 
            stop(paste(x, 
                       "must be numeric and may not contain negative numbers.", 
                       sep = " "))
    })
    return(as.data.frame(lapply(vectors, 
                                rep, 
                                length.out = max.length)))
}

#' Compress Customer-by-Sufficient-Statistic (CBS) Matrix
#'
#' Combines all customers with the same combination of recency, frequency and
#' length of calibration period in the customer-by-sufficient-statistic matrix,
#' and adds a fourth column labelled "custs" (with the number of customers
#' belonging in each row).
#'
#' This function is meant to be used to speed up log-likelihood and 
#' parameter estimation functions in the Pareto/NBD (pnbd) set of functions. 
#' How much faster those function run depends on how similar customers are. 
#' You can used compressed CBS matrices in BG/NBD estimation too, but there
#' will be no speed gains there over using un-compressed CBS data.
#'
#' This function only takes columns "x", "t.x", and "T.cal" into account. All
#' other columns will be added together - for example, if you have a spend
#' column, the output's spend column will contain the total amount spent by all
#' customers with an identical recency, frequency, and time observed.
#'
#' @param cbs 	calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x"), recency ("t.x"), and total time
#'   observed ("T.cal"). Note that recency must be the time between the start of
#'   the calibration period and the customer's last transaction, not the time
#'   between the customer's last transaction and the end of the calibration
#'   period.
#' @param rounding  	the function tries to ensure that there are similar
#'   customers by rounding the customer-by-sufficient-statistic matrix first.
#'   This parameter determines how many decimal places are left in the data.
#'   Negative numbers are allowed; see the documentation for round in the base
#'   package. As of the time of writing, that documentation states: "Rounding to
#'   a negative number of digits means rounding to a power of ten, so for
#'   example round(x, digits = -2) rounds to the nearest hundred."
#' @return  A customer-by-sufficient-statistic matrix with an additional column
#'   "custs", which contains the number of customers with each combination of
#'   recency, frequency and length of calibration period.
#' @examples
#' # Create a sample customer-by-sufficient-statistic matrix:
#' set.seed(7)
#' x <- sample(1:4, 10, replace = TRUE)
#' t.x <- sample(1:4, 10, replace = TRUE)
#' T.cal <- rep(4, 10)
#' ave.spend <- sample(10:20, 10, replace = TRUE)
#' cbs <- cbind(x, t.x, T.cal, ave.spend)
#' cbs
#'
#' # If cbs is printed, you would note that the following
#' # sets of rows have the same x, t.x and T.cal:
#' # (1, 6, 8); (3, 9)
#'
#' dc.compress.cbs(cbs, 0)   # No rounding necessary
#'
#' # Note that all additional columns (in this case, ave.spend)
#' # are aggregated by sum.
dc.compress.cbs <- function(cbs, 
                            rounding = 3) {
    
    if (!("x" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a frequency column labelled \"x\"")
    if (!("t.x" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a recency column labelled \"t.x\"")
    if (!("T.cal" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a column for length of time observed labelled \"T.cal\"")
    
    orig.rows <- nrow(cbs)
    
    if (!("custs" %in% colnames(cbs))) {
        custs <- rep(1, nrow(cbs))
        cbs <- cbind(cbs, custs)
    }
    
    other.colnames <- colnames(cbs)[!(colnames(cbs) %in% c("x", "t.x", "T.cal"))]
    
    ## Round x, t.x and T.cal to the desired level
    cbs[, c("x", "t.x", "T.cal")] <- round(cbs[, c("x", "t.x", "T.cal")], rounding)
    
    ## Aggregate every column that is not x, t.x or T.cal by those columns. Do this by
    ## summing entries which have the same x, t.x and T.cal.
    cbs <- as.matrix(aggregate(cbs[, !(colnames(cbs) %in% c("x", "t.x", "T.cal"))], 
                               by = list(x = cbs[, "x"], t.x = cbs[, "t.x"], T.cal = cbs[, "T.cal"]), sum))
    
    colnames(cbs) <- c("x", "t.x", "T.cal", other.colnames)
    final.rows <- nrow(cbs)
    message("Data reduced from ", orig.rows, " rows to ", final.rows, " rows.")
    return(cbs)
}

#' Convert Event Log to CBS and CBT Matrices
#'
#' Uses an event log to return calibration period CBT and CBS, holdout period
#' CBT and CBS, and summary data for each customer (including times of first and
#' last transactions).
#'
#' This function automatically removes customers' first transactions, meaning
#' that the output matrices will only contain repeat transaction information.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. The "date" column must
#'   contain date objects, not character strings or factors.
#' @param per   interval of time for customer-by-sufficient-statistic matrix.
#'   May be "day", "week", "month", "quarter", or "year".
#' @param T.cal     R date object indicating when the calibration period ends.
#' @param T.tot T.tot	R date object indicating when holdout period ends.
#' @param merge.same.date   If TRUE, transactions from the same period count as
#'   a single transaction instead of counting as multiple transactions.
#' @param cohort.birth.per  Time interval used to filter the event log. Can be
#'   specified as a Date object or a vector of two Dates. If one date object is
#'   used, the birth period is from the minimum date in the dataset through the
#'   given date. If two dates are given, the birth period is set between
#'   (inclusive) the two dates.
#' @param dissipate.factor  integer indicating how much of the dataset to
#'   eliminate. If left as 1, none of the dataset is eliminated.
#'   (dissipate.factor-1)/(dissipate.factor) events will be removed from the
#'   event log. For example, if 2 is provided, 1/2 of the event log is
#'   eliminated, and if 10 is provided, 9/10 of the event log is eliminated.
#' @param statistic     Determines type of CBT returned: can be: "reach",
#'   "freq", "total.spend", or "average.spend." (note: spend requires $sales
#'   column in elog)
#' @return A list of items: - `$cal` list with CBS and CBT from the calibration
#'   period - `$holdout` list with CBS and CBT from holdout period -
#'   `$cust.data` data frame with each customer's first and last transaction
#'   details
#' @examples  
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' 
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#' 
#' data <- dc.ElogToCbsCbt(elog, per="week", T.cal=as.Date("1997-09-30"))
#' @md
dc.ElogToCbsCbt <- function(elog, 
                            per = "week", 
                            T.cal = max(elog$date), 
                            T.tot = max(elog$date), 
                            merge.same.date = TRUE, 
                            cohort.birth.per = T.cal, 
                            dissipate.factor = 1, 
                            statistic = "freq") {
    
    dc.WriteLine("Started making CBS and CBT from the ELOG...")
    
    elog <- dc.FilterCustByBirth(elog, cohort.birth.per)
    if (nrow(elog) == 0) 
        stop("error caused by customer birth filtering")
    
    elog <- elog[elog$date <= T.tot, ]
    if (nrow(elog) == 0) 
        stop("error caused by holdout period end date")
    
    elog <- dc.DissipateElog(elog, dissipate.factor)
    if (nrow(elog) == 0) 
        stop("error caused by event long dissipation")
    
    if (merge.same.date) {
        elog <- dc.MergeTransactionsOnSameDate(elog)
        if (nrow(elog) == 0) 
            stop("error caused by event log merging")
    }
    
    calibration.elog <- elog[elog$date <= T.cal, ]
    holdout.elog <- elog[elog$date > T.cal, ]
    
    split.elog.list <- dc.SplitUpElogForRepeatTrans(calibration.elog)
    
    repeat.transactions.elog <- split.elog.list$repeat.trans.elog
    cust.data <- split.elog.list$cust.data
    
    
    dc.WriteLine("Started Building CBS and CBT for calibration period...")
    cbt.cal <- dc.BuildCBTFromElog(calibration.elog, statistic)
    cbt.cal.rep.trans <- dc.BuildCBTFromElog(repeat.transactions.elog, statistic)
    cbt.cal <- dc.MergeCustomers(cbt.cal, cbt.cal.rep.trans)
    
    dates <- data.frame(cust.data$birth.per, cust.data$last.date, T.cal)
    
    cbs.cal <- dc.BuildCBSFromCBTAndDates(cbt.cal, dates, per, cbt.is.during.cal.period = TRUE)
    
    dc.WriteLine("Finished building CBS and CBT for calibration period.")
    
    cbt.holdout <- NULL
    cbs.holdout <- NULL
    if (nrow(holdout.elog) > 0) {
        dc.WriteLine("Started building CBS and CBT for holdout period...")
        cbt.holdout <- dc.BuildCBTFromElog(holdout.elog, statistic)
        
        dates <- c((T.cal + 1), T.tot)
        cbs.holdout <- dc.BuildCBSFromCBTAndDates(cbt.holdout, dates, per, cbt.is.during.cal.period = FALSE)
        cbt.holdout <- dc.MergeCustomers(cbt.cal, cbt.holdout)
        cbs.holdout <- dc.MergeCustomers(cbs.cal, cbs.holdout)
        dc.WriteLine("Finished building CBS and CBT for holdout.")
        dc.WriteLine("...Finished Making All CBS and CBT")
        return(list(cal = list(cbs = cbs.cal, cbt = cbt.cal), holdout = list(cbt = cbt.holdout, 
            cbs = cbs.holdout), cust.data = cust.data))
    }
    
    dc.WriteLine("...Finished Making All CBS and CBT")
    return(list(cal = list(cbs = cbs.cal, cbt = cbt.cal), holdout = list(cbt = cbt.holdout, 
        cbs = cbs.holdout), cust.data = cust.data))
}

#' Filter Customer by Birth
#'
#' Filters an event log, keeping all transactions made by customers who made
#' their first transactions in the given time interval.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. The date column must be
#'   formatted as Date objects.
#' @param cohort.birth.per  Time interval used to filter the event log. Can be
#'   specified as a Date object or a vector of two Dates. If one date object is
#'   used, the birth period is from the minimum date in the dataset through the
#'   given date. If two dates are given, the birth period is set between
#'   (inclusive) the two dates.
#' @return event log with only rows from customers who made their first
#'   transaction within the birth period.
#' @examples
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#'
#' # converting the date column to Date objects is
#' # necessary for this function.
#' elog$date <- as.Date(elog$date, "%Y%m%d")
#'
#' # starting date. Note that it must be a Date object.
#' start.date <- as.Date("1997-01-01")
#' # ending date. Note that it must be a Date object.
#' end.date <- as.Date("1997-01-31")
#'
#' # Filter the elog to include only customers who made their
#' # first transaction in January 1997
#' filtered.elog <- dc.FilterCustByBirth(elog, c(start.date, end.date))
dc.FilterCustByBirth <- function(elog, 
                                 cohort.birth.per) {
    L = length(cohort.birth.per)
    if (L > 2) {
        stop("Invalid cohort.birth.per argument")
    }
    if (L == 0) {
        return(elog)
    }
    if (L == 1) {
        start.date <- min(elog$date)
        end.date <- cohort.birth.per
    } else if (length(cohort.birth.per) == 2) {
        start.date <- min(cohort.birth.per)
        end.date <- max(cohort.birth.per)
    }
    cbt <- dc.CreateFreqCBT(elog)
    custs.first.transaction.indices <- dc.GetFirstPurchasePeriodsFromCBT(cbt)
    custs.first.transaction.dates <- as.Date(colnames(cbt)[custs.first.transaction.indices])
    custs.in.birth.period.indices <- which(custs.first.transaction.dates >= start.date & 
        custs.first.transaction.dates <= end.date)
    custs.in.birth.period <- rownames(cbt)[custs.in.birth.period.indices]
    elog <- elog[elog$cust %in% custs.in.birth.period, ]
    dc.WriteLine("Finished filtering out customers not in the birth period.")
    return(elog)
}

#' Dissipate Event Log
#'
#' Filters an event log, keeping a fraction of the original event log.
#'
#' @param elog event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction.
#' @param dissipate.factor integer indicating how much of the dataset to
#'   eliminate. It must be greater than 1 for the function to work.
#'   (dissipate.factor-1)/(dissipate.factor) events will be removed from the
#'   event log. For example, if 2 is provided, 1/2 of the event log is
#'   eliminated, and if 10 is provided, 9/10 of the event log is eliminated.
#' @return Reduced event log.
dc.DissipateElog <- function(elog, 
                             dissipate.factor) {
    if (dissipate.factor > 1) {
        x <- rep(FALSE, dissipate.factor)
        x[1] <- TRUE
        keptIndices <- rep(x, length.out = nrow(elog))
        elog <- elog[keptIndices, ]
        elog$cust <- factor(elog$cust)
        dc.WriteLine("Finished filtering out", dissipate.factor - 1, "of every", 
            dissipate.factor, "transactions.")
    } else {
        dc.WriteLine("No dissipation requested.")
    }
    return(elog)
}

#' Split Up Event Log for Repeat Transactions
#'
#' Turns an event log into a repeat transaction event log, removing customers'
#' first transactions. Also returns a data frame with information about
#' customers' first and last transactions.
#'
#' @param elog event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. The "date" column must
#'   contain date objects, not character strings or factors.
#' @return A named list: - `repeat.trans.elog`  an event log containing only
#'   repeat transactions - `cust.data`  data frame containing the first and last
#'   transaction information for each customer
dc.SplitUpElogForRepeatTrans <- function(elog) {
    dc.WriteLine("Started Creating Repeat Purchases")
    unique.custs <- unique(elog$cust)
    first.trans.indices <- rep(0, length(unique.custs))
    last.trans.indices <- rep(0, length(unique.custs))
    count <- 0
    for (cust in unique.custs) {
        count <- count + 1
        cust.indices <- which(elog$cust == cust)
        # Of this customer's transactions, find the index of the first one
        first.trans.indices[count] <- min(cust.indices[which(elog$date[cust.indices] == 
            min(elog$date[cust.indices]))])
        
        # Of this customer's transactions, find the index of the last one
        last.trans.indices[count] <- min(cust.indices[which(elog$date[cust.indices] == 
            max(elog$date[cust.indices]))])
    }
    repeat.trans.elog <- elog[-first.trans.indices, ]
    
    first.trans.data <- elog[first.trans.indices, ]
    last.trans.data <- elog[last.trans.indices, ]
    
    
    # [-1] is because we don't want to change the column name for custs
    names(first.trans.data)[-1] <- paste("first.", names(first.trans.data)[-1], sep = "")
    names(first.trans.data)[which(names(first.trans.data) == "first.date")] <- "birth.per"
    names(last.trans.data) <- paste("last.", names(last.trans.data), sep = "")
    
    # [-1] is because we don't want to include two custs columns
    cust.data <- data.frame(first.trans.data, last.trans.data[, -1])
    names(cust.data) <- c(names(first.trans.data), names(last.trans.data)[-1])
    
    dc.WriteLine("Finished Creating Repeat Purchases")
    return(list(repeat.trans.elog = repeat.trans.elog, cust.data = cust.data))
}

#' Build Customer-by-Time Matrix from Event Log
#' 
#' Creates a customer-by-time matrix from an event log.
#' 
#' @param elog 	event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction.. For the total spend and
#'   average spend matrices, the event log must have a "sales" column. If the
#'   dates are not formatted to be in the order year-month-day, the columns of
#'   the customer-by-time matrix may not be ordered chronologically if the
#'   "date" column does not consist of date objects (R will order them
#'   alphabetically). This will cause problems with other functions, so it is
#'   better to convert the date column to date objects before running this
#'   function.
#' @param statistic     either "freq", "reach", "total.spend", or
#'   "average.spend". This determines what type of customer-by-time matrix is
#'   returned.
#' @return Customer-by-time matrix.
dc.BuildCBTFromElog <- function(elog, 
                                statistic = "freq") {
    dc.WriteLine("Started Building CBT...")
    if (statistic == "freq") {
        return(dc.CreateFreqCBT(elog))
    } else if (statistic == "reach") {
        return(dc.CreateReachCBT(elog))
    } else if (statistic == "total.spend") {
        return(dc.CreateSpendCBT(elog))
    } else if (statistic == "average.spend") {
        return(dc.CreateSpendCBT(elog, is.avg.spend = TRUE))
    } else {
        stop("Invalid cbt build (var: statistic) specified.")
    }
}

#' Create Frequency Customer-by-Time Matrix
#' 
#' Creates a customer-by-time matrix with total number of transactions per time period.
#' 
#' @param elog 	event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. If the dates are not
#'   formatted to be in the order year-month-day, the columns of the
#'   customer-by-time matrix may not be ordered chronologically if the "date"
#'   column does not consist of date objects (R will order them alphabetically).
#'   This will cause problems with other functions, so it is better to convert
#'   the date column to date objects before running this function.
#' @return  Frequency customer-by-time matrix.
#' @examples 
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' 
#' # Given that the dates are in the order year-month-day,
#' # it is not strictly necessary to convert the date column
#' # to date formats. However, it is good practice:
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#' 
#' freq.cbt <- dc.CreateFreqCBT(elog)
dc.CreateFreqCBT <- function(elog) {
    # Factoring is so that when xtabs sorts customers, it does so in the original
    # order It doesn't matter that they're factors; rownames are stored as characters
    elog$cust <- factor(elog$cust, levels = unique(elog$cust))
    xt <- xtabs(~cust + date, data = elog)
    dc.WriteLine("...Completed Freq CBT")
    return(xt)
}

#' Create Reach Customer-by-Time Matrix
#'
#' Creates a customer-by-time matrix with 1's in periods that a customer made a
#' transaction and 0's otherwise.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. If the dates are not
#'   formatted to be in the order year-month-day, the columns of the
#'   customer-by-time matrix may not be ordered chronologically if the "date"
#'   column does not consist of date objects (R will order them alphabetically).
#'   This will cause problems with other functions, so it is better to convert
#'   the date column to date objects before running this function.
#' @return  Reach customer-by-time matrix.
#' @examples 
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' 
#' # Given that the dates are in the order year-month-day,
#' # it is not strictly necessary to convert the date column
#' # to date formats. However, it is good practice:
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#' 
#' reach.cbt <- dc.CreateReachCBT(elog)
dc.CreateReachCBT <- function(elog) {
    # Factoring is so that when xtabs sorts customers, it does so in the original
    # order It doesn't matter that they're factors; rownames are stored as characters
    elog$cust <- factor(elog$cust, levels = unique(elog$cust))
    xt <- xtabs(~cust + date, data = elog)
    xt[xt > 1] <- 1
    dc.WriteLine("...Completed Reach CBT")
    return(xt)
}

#' Create Spend Customer-by-Time Matrix
#'
#' Creates a customer-by-time matrix with spend per time period.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. If the dates are not
#'   formatted to be in the order year-month-day, the columns of the
#'   customer-by-time matrix may not be ordered chronologically if the "date"
#'   column does not consist of date objects (R will order them alphabetically).
#'   This will cause problems with other functions, so it is better to convert
#'   the date column to date objects before running this function.
#' @param is.avg.spend  if TRUE, return average spend customer-by-time matrix;
#'   else, return total spend customer-by-time matrix.
#' @return  Spend customer-by-time matrix.
#' @examples
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5);
#'
#' # Given that the dates are in the order year-month-day,
#' # it is not strictly necessary to convert the date column
#' # to date formats. However, it is good practice:
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#'
#' spend.cbt <- dc.CreateSpendCBT(elog)
dc.CreateSpendCBT <- function(elog, 
                              is.avg.spend = FALSE) {
    # Factoring is so that when xtabs sorts customers, it does so in the original
    # order It doesn't matter that they're factors; rownames are stored as characters
    elog$cust <- factor(elog$cust, levels = unique(elog$cust))
    sales.xt <- xtabs(sales ~ cust + date, data = elog)
    if (is.avg.spend) {
        suppressMessages(freq.cbt <- dc.CreateFreqCBT(elog))
        sales.xt <- sales.xt/freq.cbt
        # For the cases where there were no transactions
        sales.xt[which(!is.finite(sales.xt))] <- 0
    }
    dc.WriteLine("...Completed Spend CBT")
    return(sales.xt)
}

#' Make Recency-Frequency Matrix Skeleton
#'
#' Creates a matrix with all possible recency and frequency combinations.
#'
#' Makes the structure in which to input data for recency-frequency matrices.
#'
#' @param n.periods     number of transaction opportunities in the calibration
#'   period.
#' @return  Matrix with two columns: frequency ("x") and recency ("t.x"). All
#'   possible recency-frequency combinations in the calibration period are
#'   represented.
dc.MakeRFmatrixSkeleton <- function(n.periods) {
    ## note: to access the starting i'th t.x element use (i>0): i*(i-1)/2 + 2, ...
    ## this yields the sequence: 2, 3, 5, 8, ... there are n*(n+1)/2 + 1 elements in
    ## this table
    n <- n.periods
    rf.mx.skeleton <- matrix(0, n * (n + 1)/2 + 1, 2)
    colnames(rf.mx.skeleton) <- c("x", "t.x")
    for (ii in 1:n) {
        ith.t.index <- 2 + ii * (ii - 1)/2
        t.vector <- rep(ii, ii)
        x.vector <- c(1:ii)
        rf.mx.skeleton[ith.t.index:(ith.t.index + (ii - 1)), 1] <- x.vector
        rf.mx.skeleton[ith.t.index:(ith.t.index + (ii - 1)), 2] <- t.vector
    }
    return(rf.mx.skeleton)
}

#' Make Holdout Period Recency-Frequency Matrix
#'
#' Creates a recency-frequency matrix for the holdout period.
#'
#' @param holdout.cbt   holdout period frequency customer-by-time matrix. This
#'   is a matrix consisting of a row per customer and a column per time period.
#'   It should contain the number of transactions each customer made per time
#'   period.
#' @return  recency-frequency matrix for the holdout period, with three columns:
#'   frequency ("x.star"), recency ("t.x.star"), number of transaction
#'   opportunities in the holdout period ("n.star"), and the number of customers
#'   with each frequency-recency combination ("custs").
dc.MakeRFmatrixHoldout <- function(holdout.cbt) {
    
    holdout.length <- ncol(holdout.cbt)
    matrix.skeleton <- dc.MakeRFmatrixSkeleton(holdout.length)
    n.combinations <- nrow(matrix.skeleton)
    n.star <- rep(holdout.length, n.combinations)
    final.transactions <- dc.GetLastPurchasePeriodsFromCBT(holdout.cbt)
    custs <- rep(0, n.combinations)
    for (ii in 1:n.combinations) {
        custs.with.freq <- which(rowSums(holdout.cbt) == matrix.skeleton[ii, 1])
        custs.with.rec <- which(final.transactions == matrix.skeleton[ii, 2])
        custs[ii] <- length(intersect(custs.with.freq, custs.with.rec))
    }
    rf.holdout.matrix <- cbind(matrix.skeleton, n.star, custs)
    colnames(rf.holdout.matrix) <- c("x.star", "t.x.star", "n.star", "custs")
    return(rf.holdout.matrix)
}

#' Make Calibration Period Recency-Frequency Matrix
#'
#' Make a calibration period recency-frequency matrix.
#'
#' @param frequencies   vector which indicates the number of repeat transactions
#'   made by customers in the calibration period.
#' @param periods.of.final.purchases    a vector indicating in which period
#'   customers made their final purchases.
#' @param num.of.purchase.periods   	the number of transaction opportunities in
#'   the calibration period.
#' @param holdout.frequencies   an optional vector indicating the number of
#'   transactions made by customers in the holdout period.
#' @return  A matrix with all possible frequency-recency combinations, and the
#'   number of customers with each combination. It contains columns for
#'   frequency ("x"), recency ("t.x"), number of transaction opportunities in
#'   the calibration period ("n.cal"), number of customers with this combination
#'   of recency, frequency, and number of periods observed ("custs"), and
#'   optionally, number of transactions in the holdout period ("x.star").
#' @examples
#' elog <- dc.ReadLines(system.file("data/discreteSimElog.csv", package="BTYD"),1,2)
#' elog[,"date"] <- as.Date(elog[,"date"])
#'
#' cutoff.date <- as.Date("1977-01-01")
#' cbt <- dc.CreateReachCBT(elog)
#' cal.cbt <- cbt[,as.Date(colnames(cbt)) <= cutoff.date]
#' holdout.cbt <- cbt[,as.Date(colnames(cbt)) > cutoff.date]
#'
#' cal.start.dates.indices <- dc.GetFirstPurchasePeriodsFromCBT(cal.cbt)
#' cal.start.dates <- as.Date(colnames(cal.cbt)[cal.start.dates.indices])
#' cal.end.dates.indices <- dc.GetLastPurchasePeriodsFromCBT(cal.cbt)
#' cal.end.dates <- as.Date(colnames(cal.cbt)[cal.end.dates.indices])
#' T.cal.total <- rep(cutoff.date, nrow(cal.cbt))
#' cal.dates <- data.frame(cal.start.dates, cal.end.dates, T.cal.total)
#'
#' # Create calibration period customer-by-sufficient-statistic data frame,
#' # using years as the unit of time.
#' cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.cbt,
#'                                       cal.dates,
#'                                       per="year",
#'                                       cbt.is.during.cal.period=TRUE)
#'
#' holdout.start <- as.Date(colnames(holdout.cbt)[1])
#' holdout.end <- as.Date(tail(colnames(holdout.cbt),n=1))
#' # The (-1) below is to remove the effect of the birth period - we are only
#' # interested in repeat transactions in the calibration period.
#' frequencies <- (cal.cbs[,"x"] - 1)
#' periods.of.final.purchases <- cal.cbs[,"t.x"]
#' num.of.purchase.periods <- ncol(cal.cbt) - 1
#'
#' # Create a calibration period recency-frequency matrix
#' cal.rf.matrix <- dc.MakeRFmatrixCal(frequencies,
#'                                     periods.of.final.purchases,
#'                                     num.of.purchase.periods)
dc.MakeRFmatrixCal <- function(frequencies, 
                               periods.of.final.purchases, 
                               num.of.purchase.periods, 
                               holdout.frequencies = NULL) {
    
    if (!is.numeric(periods.of.final.purchases)) {
        stop("periods.of.final.purchases must be numeric")
    }
    if (length(periods.of.final.purchases) != length(frequencies)) {
        stop(paste("number of customers in frequencies is not equal", "to the last purchase period vector"))
    }
    ## initializes the data structures to later be filled in with counts
    rf.mx.skeleton <- dc.MakeRFmatrixSkeleton(num.of.purchase.periods)
    if (is.null(holdout.frequencies)) {
        RF.matrix <- cbind(rf.mx.skeleton, num.of.purchase.periods, 0)
        colnames(RF.matrix) <- c("x", "t.x", "n.cal", "custs")
    } else {
        RF.matrix <- cbind(rf.mx.skeleton, num.of.purchase.periods, 0, 0)
        colnames(RF.matrix) <- c("x", "t.x", "n.cal", "custs", "x.star")
    }
    
    
    ## create a matrix out of the frequencies & periods.of.final.purchases
    rf.n.custs <- cbind(frequencies, periods.of.final.purchases, holdout.frequencies)
    ## count all the pairs with zero for frequency and remove them
    zeroes.rf.subset <- which(rf.n.custs[, 1] == 0)  ##(which x == 0)
    RF.matrix[1, 4] <- length(zeroes.rf.subset)
    if (!is.null(holdout.frequencies)) {
        RF.matrix[1, 5] <- sum(holdout.frequencies[zeroes.rf.subset])
    }
    rf.n.custs <- rf.n.custs[-zeroes.rf.subset, ]
    
    ## sort the count data by both frequency and final purchase period
    rf.n.custs <- rf.n.custs[order(rf.n.custs[, 1], rf.n.custs[, 2]), ]
    
    ## formula: (x-1) + 1 + tx*(tx-1)/2 + 1 keep count of duplicates once different,
    ## use formula above to place count into the RF table.
    current.pair <- c(rf.n.custs[1, 1], rf.n.custs[1, 2])
    
    same.item.in.a.row.counter <- 1
    if (!is.null(holdout.frequencies)) {
        x.star.total <- rf.n.custs[1, 3]
    }
    num.count.points <- nrow(rf.n.custs)
    for (ii in 2:num.count.points) {
        last.pair <- current.pair
        current.pair <- c(rf.n.custs[ii, 1], rf.n.custs[ii, 2])
        if (identical(last.pair, current.pair)) {
            same.item.in.a.row.counter <- same.item.in.a.row.counter + 1
            if (!is.null(holdout.frequencies)) {
                x.star.total <- x.star.total + rf.n.custs[ii, 3]
            }
        } else {
            x <- last.pair[1]
            t.x <- last.pair[2]
            corresponding.rf.index <- (x - 1) + 1 + t.x * (t.x - 1)/2 + 1
            RF.matrix[corresponding.rf.index, 4] <- same.item.in.a.row.counter
            same.item.in.a.row.counter <- 1
            if (!is.null(holdout.frequencies)) {
                RF.matrix[corresponding.rf.index, 5] <- x.star.total
                x.star.total <- rf.n.custs[ii, 3]
            }
        }
        if (ii == num.count.points) {
            x <- current.pair[1]
            t.x <- current.pair[2]
            corresponding.rf.index <- (x - 1) + 1 + t.x * (t.x - 1)/2 + 1
            RF.matrix[corresponding.rf.index, 4] <- same.item.in.a.row.counter
            same.item.in.a.row.counter <- NULL
            if (!is.null(holdout.frequencies)) {
                RF.matrix[corresponding.rf.index, 5] <- x.star.total
                x.star.total = NULL
            }
        }
    }
    return(RF.matrix)
}

#' Build CBS matrix from CBT matrix
#'
#' Given a customer-by-time matrix, yields the resulting
#' customer-by-sufficient-statistic matrix.
#'
#' The customer-by-sufficient statistic matrix will contain the sum of the
#' statistic included in the customer-by-time matrix (see the cbt parameter),
#' the customer's last transaction date, and the total time period for which the
#' customer was observed.
#'
#' @param cbt 	customer-by-time matrix. This is a matrix consisting of a row per
#'   customer and a column per time period. It should contain numeric
#'   information about a customer's transactions in every time period - either
#'   the number of transactions in that time period (frequency), a 1 to indicate
#'   that at least 1 transaction occurred (reach), or the average/total amount
#'   spent in that time period.
#' @param dates     if cbt.is.during.cal.period is TRUE, then dates is a data
#'   frame with three columns: 1. the dates when customers made their first
#'   purchases 2. the dates when customers made their last purchases 3. the date
#'   of the end of the calibration period. if cbt.is.during.cal.period is FALSE,
#'   then dates is a vector with two elements: 1. the date of the beginning of
#'   the holdout period 2. the date of the end of the holdout period.
#' @param per   interval of time for customer-by-sufficient-statistic matrix.
#'   May be "day", "week", "month", "quarter", or "year".
#' @param cbt.is.during.cal.period  if TRUE, indicates the customer-by-time
#'   matrix is from the calibration period. If FALSE, indicates the
#'   customer-by-time matrix is from the holdout period.
#' @return Customer-by-sufficient-statistic matrix, with three columns:
#'   frequency("x"), recency("t.x") and total time observed("T.cal"). See
#'   details. Frequency is total transactions, not repeat transactions.
#' @examples
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#' 
#' # Transaction-flow models are about interpurchase times. Since we
#' # only know purchase times to the day, we merge all transaction on
#' # the same day. This example uses dc.MergeTransactionsOnSameDate
#' # to illustrate this; however, we could have simply used dc.CreateReachCBT
#' # instead of dc.CreateFreqCBT to obtain the same result.
#' merged.elog <- dc.MergeTransactionsOnSameDate(elog)
#' cutoff.date <- as.Date("1997-09-30")
#' freq.cbt <- dc.CreateFreqCBT(merged.elog)
#' cal.freq.cbt <- freq.cbt[,as.Date(colnames(freq.cbt)) <= cutoff.date]
#' holdout.freq.cbt <- freq.cbt[,as.Date(colnames(freq.cbt)) > cutoff.date]
#' 
#' cal.start.dates.indices <- dc.GetFirstPurchasePeriodsFromCBT(cal.freq.cbt)
#' cal.start.dates <- as.Date(colnames(cal.freq.cbt)[cal.start.dates.indices])
#' cal.end.dates.indices <- dc.GetLastPurchasePeriodsFromCBT(cal.freq.cbt)
#' cal.end.dates <- as.Date(colnames(cal.freq.cbt)[cal.end.dates.indices])
#' T.cal.total <- rep(cutoff.date, nrow(cal.freq.cbt))
#' cal.dates <- data.frame(cal.start.dates, 
#'                         cal.end.dates, 
#'                         T.cal.total)
#' 
#' # Create calibration period customer-by-sufficient-statistic data frame,
#' # using weeks as the unit of time.
#' cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.freq.cbt, 
#'                                       cal.dates,
#'                                       per="week", 
#'                                       cbt.is.during.cal.period=TRUE)
#' # Force the calibration period customer-by-sufficient-statistic to only contain
#' # repeat transactions (required by BG/BB and Pareto/NBD models)
#' cal.cbs[,"x"] <- cal.cbs[,"x"] - 1
#' 
#' holdout.start <- cutoff.date+1
#' holdout.end <- as.Date(colnames(holdout.freq.cbt)[ncol(holdout.freq.cbt)])
#' holdout.dates <- c(holdout.start, holdout.end)
#' 
#' # Create holdout period customer-by-sufficient-statistic data frame, using weeks
#' # as the unit of time.
#' holdout.cbs <- dc.BuildCBSFromCBTAndDates(holdout.freq.cbt, 
#'                                           holdout.dates,
#'                                           per="week", 
#'                                           cbt.is.during.cal.period=FALSE)                        
#' @md
dc.BuildCBSFromCBTAndDates <- function(cbt, 
                                       dates, 
                                       per, 
                                       cbt.is.during.cal.period = TRUE) {
    if (cbt.is.during.cal.period == TRUE) {
        dc.WriteLine("Started making calibration period CBS...")
        custs.first.dates <- dates[, 1]
        custs.last.dates <- dates[, 2]
        T.cal <- dates[, 3]
        if (length(custs.first.dates) != length(custs.last.dates)) {
            stop("Invalid dates (different lengths) in BuildCBSFromFreqCBTAndDates")
        }
        
        f <- rowSums(cbt)
        r <- as.numeric(difftime(custs.last.dates, custs.first.dates, units = "days"))
        T <- as.numeric(difftime(T.cal, custs.first.dates, units = "days"))
        x <- switch(per, day = 1, week = 7, month = 365/12, quarter = 365/4, year = 365)
        r = r/x
        T = T/x
        cbs = cbind(f, r, T)
        # cbs <- data.frame(f=f, r=r/x, T=T/x)
        rownames(cbs) <- rownames(cbt)
        colnames(cbs) <- c("x", "t.x", "T.cal")
    } else {
        ## cbt is during holdout period
        dc.WriteLine("Started making holdout period CBS...")
        date.begin.holdout.period <- dates[1]
        date.end.holdout.period <- dates[2]
        f <- rowSums(cbt)
        T <- as.numeric(difftime(date.end.holdout.period, date.begin.holdout.period, 
            units = "days")) + 1
        x <- switch(per, day = 1, week = 7, month = 365/12, quarter = 365/4, year = 365)
        T = T/x
        cbs = cbind(f, T)
        # cbs <- data.frame( f=f, T=T/x)
        rownames(cbs) <- rownames(cbt)
        colnames(cbs) <- c("x.star", "T.star")
    }
    
    dc.WriteLine("Finished building CBS.")
    return(cbs)
}

#' Merge Customers
#'
#' Takes two CBT or CBS matrices and ensures that the second one has the same
#' row names as the first.
#'
#' Care should be taken in using this function. It inserts zero values in all
#' rows that were not in the original holdout period data. This behavior does
#' not cause a problem if using CBT matrices, but will cause a problem if using
#' CBS matrices (for example, the output will report all customers with a
#' holdout period length of zero). However, this particular issue is easily
#' fixed (see examples) and should not cause problems.
#'
#' A work-around to avoid using this function is presented in the example for
#' [`dc.BuildCBSFromCBTAndDates`] - build the full CBT and only use the columns
#' applying to each particular time period to construct separate CBTs, and from
#' them, CBSs. That is a much cleaner and less error-prone method; however, on
#' occasion the data will not be available in event log format and you may not
#' be able to construct a CBT for both time periods together.
#'
#' @param data.correct  CBT or CBS with the correct customer IDs as row names.
#'   Usually from the calibration period.
#' @param data.to.correct   CBT or CBS which needs to be fixed (customer IDs
#'   inserted). Usually from the holdout period.
#' @return Updated holdout period CBT or CBS.
#' @examples
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#' cutoff.date <- as.Date("1997-09-30")
#' cal.elog <- elog[which(elog[,"date"] <= cutoff.date),]
#' holdout.elog <- elog[which(elog[,"date"] > cutoff.date),]
#'
#' # Create calibration period CBT from cal.elog
#' cal.reach.cbt <- dc.CreateReachCBT(cal.elog)
#' # Create holdout period CBT from holdout.elog
#' holdout.reach.cbt <- dc.CreateReachCBT(holdout.elog)
#'
#' # Note the difference:
#' nrow(cal.reach.cbt)            # 2357 customers
#' nrow(holdout.reach.cbt)        # 684 customers
#'
#' # Create a "fixed" holdout period CBT, with the same number
#' # of customers in the same order as the calibration period CBT
#' fixed.holdout.reach.cbt <- dc.MergeCustomers(cal.reach.cbt, holdout.reach.cbt)
#' nrow(fixed.holdout.reach.cbt)  # 2357 customers
#'
#' # You can verify that the above is correct by turning these into a CBS
#' # (see \code{\link{dc.BuildCBSFromCBTAndDates}} and using
#' # \code{\link{pnbd.PlotFreqVsConditionalExpectedFrequency}}, for example
#'
#' # Alternatively, we can fix the CBS, instead of the CBS:
#'
#' cal.start.dates.indices <- dc.GetFirstPurchasePeriodsFromCBT(cal.reach.cbt)
#' cal.start.dates <- as.Date(colnames(cal.reach.cbt)[cal.start.dates.indices])
#' cal.end.dates.indices <- dc.GetLastPurchasePeriodsFromCBT(cal.reach.cbt)
#' cal.end.dates <- as.Date(colnames(cal.reach.cbt)[cal.end.dates.indices])
#' T.cal.total <- rep(cutoff.date, nrow(cal.reach.cbt))
#' cal.dates <- data.frame(cal.start.dates, cal.end.dates, T.cal.total)
#'
#' # Create calibration period customer-by-sufficient-statistic data frame,
#' # using weeks as the unit of time.
#' cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.reach.cbt,
#'                                       cal.dates,
#'                                       per="week",
#'                                       cbt.is.during.cal.period=TRUE)
#'
#' # Force the calibration period customer-by-sufficient-statistic to only
#' # 	contain repeat transactions (required by BG/BB and Pareto/NBD models)
#' cal.cbs[,"x"] <- cal.cbs[,"x"] - 1
#'
#' holdout.start <- cutoff.date+1
#' holdout.end <- as.Date(colnames(fixed.holdout.reach.cbt)[ncol(fixed.holdout.reach.cbt)])
#' holdout.dates <- c(holdout.start, holdout.end)
#'
#' # Create holdout period customer-by-sufficient-statistic data frame,
#' # using weeks as the unit of time.
#' holdout.cbs <- dc.BuildCBSFromCBTAndDates(holdout.reach.cbt,
#'                                           holdout.dates,
#'                                           per="week",
#'                                           cbt.is.during.cal.period=FALSE)
#'
#' # Note the difference:
#' nrow(cal.cbs)            # 2357 customers
#' nrow(holdout.cbs)        # 684 customers
#'
#' # Create a "fixed" holdout period CBS, with the same number
#' # of customers in the same order as the calibration period CBS
#' fixed.holdout.cbs <- dc.MergeCustomers(cal.cbs, holdout.cbs)
#' nrow(fixed.holdout.cbs)  # 2357 customers
#'
#' # Furthermore, this function will assign a zero value to all fields
#' # that were not in the original holdout period CBS. Since T.star is the
#' # same for all customers in the holdout period, we should fix that:
#' fixed.holdout.cbs[,"T.star"] <- rep(max(fixed.holdout.cbs[,"T.star"]),nrow(fixed.holdout.cbs))
dc.MergeCustomers <- function(data.correct, 
                              data.to.correct) {
    
    ## Initialize a new data frame
    data.to.correct.new <- matrix(0, nrow = nrow(data.correct), ncol = ncol(data.to.correct))
    # data.to.correct.new <- data.frame(data.to.correct.new.size)
    orig.order <- 1:nrow(data.correct)
    orig.order <- orig.order[order(rownames(data.correct))]
    data.correct.ordered <- data.correct[order(rownames(data.correct)), ]
    ## obscure code: handles boundary case when data.correct has one column and
    ## coerces data.correct.ordered to be a vector
    if (is.null(nrow(data.correct.ordered))) {
        # data.correct.ordered <- data.frame(data.correct.ordered)
        rownames(data.correct.ordered) <- rownames(data.correct)[order(rownames(data.correct))]
        colnames(data.correct.ordered) <- colnames(data.correct)
    }
    
    data.to.correct <- data.to.correct[order(rownames(data.to.correct)), ]
    rownames(data.to.correct.new) <- rownames(data.correct.ordered)
    colnames(data.to.correct.new) <- colnames(data.to.correct)
    
    ## Initialize the two iterators ii.correct, ii.to.correct
    ii.correct <- 1
    ii.to.correct <- 1
    
    ## Grab the data to hold the stopping conditions
    max.correct.iterations <- nrow(data.correct.ordered)
    max.to.correct.iterations <- nrow(data.to.correct)
    
    ## Grab the lists of customers from the data frames and convert them to optimize
    ## the loop speed
    cust.list.correct <- rownames(data.correct.ordered)
    cust.list.to.correct <- rownames(data.to.correct)
    
    cust.correct.indices <- c()
    cust.to.correct.indices <- c()
    
    
    while (ii.correct <= max.correct.iterations & ii.to.correct <= max.to.correct.iterations) {
        cur.cust.correct <- cust.list.correct[ii.correct]
        cur.cust.to.correct <- cust.list.to.correct[ii.to.correct]
        if (cur.cust.correct < cur.cust.to.correct) {
            ii.correct <- ii.correct + 1
        } else if (cur.cust.correct > cur.cust.to.correct) {
            ii.to.correct <- ii.to.correct + 1
        } else if (cur.cust.correct == cur.cust.to.correct) {
            ## data.to.correct.new[ii.correct, ] = data.to.correct[ii.to.correct, ]
            cust.correct.indices <- c(cust.correct.indices, ii.correct)
            cust.to.correct.indices <- c(cust.to.correct.indices, ii.to.correct)
            
            ii.correct <- ii.correct + 1
            ii.to.correct <- ii.to.correct + 1
        } else {
            stop("Array checking error in MergeCustomers")
        }
    }
    data.to.correct.new[cust.correct.indices, ] <- data.to.correct
    data.to.correct.new <- data.to.correct.new[order(orig.order), ]
    return(data.to.correct.new)
}

#' Merge Transactions on Same Day
#'
#' Updates an event log; any transactions made by the same customer on the same
#' day are combined into one transaction.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction.
#' @return  Event log with transactions made by the same customer on the same
#'   day merged into one transaction.
dc.MergeTransactionsOnSameDate <- function(elog) {
    dc.WriteLine("Started merging same-date transactions...")
    elog <- cbind(elog, 1:nrow(elog) * (!duplicated(elog[, c("cust", "date")])))
    aggr.elog <- aggregate(elog[, !(colnames(elog) %in% c("cust", "date"))], by = list(cust = elog[, 
                                                                                                   "cust"], date = elog[, "date"]), sum)
    aggr.elog <- aggr.elog[order(aggr.elog[, ncol(aggr.elog)]), ][, -ncol(aggr.elog)]
    dc.WriteLine("... Finished merging same-date transactions.")
    return(aggr.elog)
}

#' Remove Time Between
#'
#' This function creates a new event log, with time in the middle removed. Used,
#' for example, in sports with off-seasons.
#'
#' The four date parameters must be in ascending order.
#'
#' @param elog  event log, which is a data frame with columns for customer ID
#'   ("cust"), date ("date"), and optionally other columns such as "sales". Each
#'   row represents an event, such as a transaction. The "date" column must
#'   consist of date objects, not character strings.
#' @param day1  date of beginning of first period. Must be a date object.
#' @param day2 	date of end of first period. Must be a date object.
#' @param day3  date of beginning of second period. Must be a date object.
#' @param day4  date of third period. Must be a date object.
#' @return list - `elog1` the event log with all elog$date entries between day1
#'   and day2 - `elog2` the event with all elog$date entries between day3 and
#'   day4 - `elog3` elog1 combined with elog2, with all dates from elog2 reduced
#'   by the time removed between elog1 and elog2
#' @examples
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' elog[,"date"] <- as.Date(elog[,"date"], "%Y%m%d")
#'
#' # Use the cdnow data to return a 6 month event log for January, February,
#' # March, October, November, December.
#' period.one.start <- as.Date("1997-01-01")
#' period.one.end <- as.Date("1997-03-31")
#' period.two.start <- as.Date("1997-10-01")
#' period.two.end <- as.Date("1997-12-31")
#' reduced.elog <- dc.RemoveTimeBetween(elog, period.one.start, period.one.end,
#'                                      period.two.start, period.two.end)
#'
#' # Note that the new elog will go up to June 30 at a maximum, since we
#' # are only using 6 months of data starting on January 1
#' max(reduced.elog$elog3$date)  # "1997-06-30"
#' @md
dc.RemoveTimeBetween <- function(elog, 
                                 day1, 
                                 day2, 
                                 day3, 
                                 day4) {
    if (day1 > day2 || day2 > day3 || day3 > day4) {
        stop("Days are not input in increasing order.")
    }
    elog1 <- elog[which(elog$date >= day1 & elog$date <= day2), ]
    elog2 <- elog[which(elog$date >= day3 & elog$date <= day4), ]
    time.between.periods <- as.numeric(day3 - day2)
    
    elog2timeErased <- elog2
    elog2timeErased$date <- elog2$date - time.between.periods
    elog3 = rbind(elog1, elog2timeErased)
    
    elogsToReturn = list()
    elogsToReturn$elog1 <- elog1
    elogsToReturn$elog2 <- elog2
    elogsToReturn$elog3 <- elog3
    return(elogsToReturn)
}

#' Get First Purchase Periods from Customer-by-Time Matrix
#'
#' Uses a customer-by-time matrix to return a vector containing the periods in
#' which customers made their first purchase.
#'
#' @param cbt   	customer-by-time matrix. This is a matrix consisting of a row
#'   per customer and a column per time period. It should contain numeric
#'   information about a customer's transactions in every time period - either
#'   the number of transactions in that time period (frequency), a 1 to indicate
#'   that at least 1 transaction occurred (reach), or the average/total amount
#'   spent in that time period.
#' @return a vector containing the indices of periods in which customers made
#'   their first transactions. To convert to actual dates (if your
#'   customer-by-time matrix has dates as column names), use
#'   colnames(cbt)\[RESULT\]
dc.GetFirstPurchasePeriodsFromCBT <- function(cbt) {
    cbt <- as.matrix(cbt)
    num.custs <- nrow(cbt)
    num.periods <- ncol(cbt)
    first.periods <- c(num.custs)
    
    ## loops through the customers and periods and locates the first purchase periods
    ## of each customer. Records them in first.periods
    for (ii in 1:num.custs) {
        curr.cust.transactions <- as.numeric(cbt[ii, ])
        transaction.index <- 1
        made.purchase <- FALSE
        while (made.purchase == FALSE & transaction.index <= num.periods) {
            if (curr.cust.transactions[transaction.index] > 0) {
                made.purchase <- TRUE
            } else {
                transaction.index <- transaction.index + 1
            }
        }
        if (made.purchase == FALSE) {
            first.periods[ii] <- 0
        } else {
            first.periods[ii] <- transaction.index
        }
    }
    return(first.periods)
}

#' Get Last Purchase Periods from Customer-by-Time Matrix
#'
#' Uses a customer-by-time matrix to return a vector containing the periods in
#' which customers made their last purchase.
#'
#' @inheritParams   dc.GetFirstPurchasePeriodsFromCBT
#' @return  a vector containing the indices of periods in which customers made
#'   their last transactions. To convert to actual dates (if your
#'   customer-by-time matrix has dates as column names), use
#'   colnames(cbt)\[RESULT\]
dc.GetLastPurchasePeriodsFromCBT <- function(cbt) {
    cbt <- as.matrix(cbt)
    num.custs <- nrow(cbt)
    num.periods <- ncol(cbt)
    last.periods <- c(num.custs)
    
    ## loops through the customers and periods and locates the first purchase periods
    ## of each customer. Records them in last.periods
    for (ii in 1:num.custs) {
        curr.cust.transactions <- as.numeric(cbt[ii, ])
        transaction.index <- num.periods
        made.purchase <- FALSE
        while (made.purchase == FALSE & transaction.index >= 1) {
            if (curr.cust.transactions[transaction.index] > 0) {
                made.purchase <- TRUE
            } else {
                transaction.index <- transaction.index - 1
            }
        }
        if (made.purchase == FALSE) {
            last.periods[ii] <- 0
        } else {
            last.periods[ii] <- transaction.index
        }
    }
    return(last.periods)
}

#' Write Line
#' 
#' Writes any number of arguments to the console.
#' 
#' The code is literally: cat(..., fill=TRUE); flush.console();
#' 
#' @param ... objects to print to the R console.
dc.WriteLine <- function(...) {
    message(...)
    flush.console()
}

#' Add Logs 
#' 
#' Given log(a) and log(b), returns log(a + b)
#' 
#' @param loga 	first number in log space.
#' @param logb first number in log space.
addLogs <- function(loga, 
                    logb) {
    return(logb + log(exp(loga - logb) + 1))
}

#' Subtract Logs 
#' 
#' Given log(a) and log(b), returns log(a - b)
#' 
#' @inheritParams addLogs
subLogs <- function(loga, 
                    logb) {
    return(logb + log(exp(loga - logb) - 1))
}

#' Plot Log-Likelihood Contours
#'
#' Creates a set of contour plots, such that there is a contour plot for every
#' pair of parameters varying.
#'
#' For each contour plot, the non-varying parameters are kept constant at the
#' predicted values.
#'
#' The contour will extend out by (n.divs * zoom.percent) in both directions and
#' both dimensions from the estimated parameter values. The exception is if
#' allow.neg.params is FALSE. In this case, the contour plot will end at zero if
#' it would have extended into negative parameter values.
#'
#' The estimated parameter values will be indicated by the intersection of two
#' red lines.
#' 
#' @inheritParams dc.PlotLogLikelihoodContour
#' @param multiple.screens  if TRUE, plots each contour plot on a separate R graphics window.
#' @seealso [`dc.PlotLogLikelihoodContour`]
#' @examples
#' # **Example for BG/BB model:
#' data(donationsSummary)
#' rf.matrix <- donationsSummary$rf.matrix
#' 
#' # starting-point parameters
#' bgbb.startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' bgbb.est.params <- bgbb.EstimateParameters(rf.matrix, bgbb.startingparams)
#' 
#' # set up parameter names for a more descriptive result
#' bgbb.param.names <- c("alpha", "beta", "gamma", "delta")
#' 
#' # plot-log likelihood contours:
#' dc.PlotLogLikelihoodContours(bgbb.rf.matrix.LL,
#'                              bgbb.est.params, 
#'                              rf.matrix = rf.matrix, 
#'                              n.divs = 5, 
#'                              num.contour.lines = 8, 
#'                              zoom.percent = 0.3, 
#'                              allow.neg.params = FALSE, 
#'                              param.names = bgbb.param.names)
#' 
#' # **Example for Pareto/NBD model:
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' 
#' # Speed up calculations:
#' cbs <- dc.compress.cbs(cbs)
#' 
#' # parameters estimated using pnbd.EstimateParameters
#' pnbd.est.params <- cdnowSummary$est.params
#' 
#' # set up parameter names for a more descriptive result
#' pnbd.param.names <- c("r", "alpha", "s", "beta")
#' 
#' # plot log-likelihood contours:
#' dc.PlotLogLikelihoodContours(pnbd.cbs.LL, 
#'                              pnbd.est.params, 
#'                              cal.cbs = cbs, 
#'                              hardie = TRUE,
#'                              n.divs = 5, 
#'                              num.contour.lines = 15, 
#'                              zoom.percent = 0.3, 
#'                              allow.neg.params = FALSE, 
#'                              param.names = pnbd.param.names)
#' 
#' # **Example for BG/NBD model:
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' 
#' # parameters estimated using bgnbd.EstimateParameters
#' bgnbd.est.params <- cdnowSummary$est.params
#' 
#' # set up parameter names for a more descriptive result
#' bgnbd.param.names <- c("r", "alpha", "s", "beta")
#' 
#' # plot log-likelihood contours:
#' dc.PlotLogLikelihoodContours(bgnbd.cbs.LL, 
#'                              bgnbd.est.params, 
#'                              cal.cbs = cbs, 
#'                              n.divs = 5, 
#'                              num.contour.lines = 15, 
#'                              zoom.percent = 0.3, 
#'                              allow.neg.params = FALSE, 
#'                              param.names = bgnbd.param.names)
#' @md
dc.PlotLogLikelihoodContours <- function(loglikelihood.fcn, 
                                         predicted.params, 
                                         ..., 
                                         n.divs = 2, 
                                         multiple.screens = FALSE, 
                                         num.contour.lines = 10, 
                                         zoom.percent = 0.9, 
                                         allow.neg.params = FALSE, 
                                         param.names = c("param 1", "param 2", "param 3", "param 4")) {
    permutations <- combn(length(predicted.params), 2)
    num.permutations <- ncol(permutations)
    contour.plots <- list()
    
    if (multiple.screens == FALSE) {
        dev.new()
        plot.window.num.cols <- ceiling(num.permutations/2)
        plot.window.num.rows <- 2
        par(mfrow = c(plot.window.num.rows, plot.window.num.cols))
    }
    
    for (jj in 1:num.permutations) {
        vary.or.fix.param <- rep("fix", 4)
        vary.or.fix.param[permutations[, jj]] <- "vary"
        contour.plots[[jj]] <- dc.PlotLogLikelihoodContour(loglikelihood.fcn, vary.or.fix.param, 
            predicted.params, ..., n.divs = n.divs, new.dev = multiple.screens, num.contour.lines = num.contour.lines, 
            zoom.percent = zoom.percent, allow.neg.params = allow.neg.params, param.names = param.names)
    }
    
    if (multiple.screens == FALSE) {
        par(mfrow = c(1, 1))
    }
    
    return(contour.plots)
}

#' Plot Log-Likelihood Contour
#'
#' Makes a contour plot of a loglikelihood function that varies over two
#' designated parameters, centered around a set of previously estimated
#' parameters.
#'
#' The contour plot will have the first parameter labelled "vary" on the x-axis,
#' and the second parameter labelled "vary" on the y-axis. It will extend out by
#' (n.divs * zoom.percent) in both directions and both dimensions from the
#' estimated parameter values. The exception is if allow.neg.params is FALSE. In
#' this case, the contour plot will end at zero if it would have extended into
#' negative parameter values.
#'
#' The estimated parameter values will be indicated by the intersection of two
#' red lines.
#'
#' @param loglikelihood.fcn     log-likelihood function to plot.
#' @param vary.or.fix.param     a vector of strings containing either "vary" or
#'   "fix". The parameters in the same indices as "vary" will be plotted while
#'   the other parameters will remain fixed at the estimated values. See
#'   details.
#' @param predicted.params  estimated parameters.
#' @param ...   all additional arguments required by the log-likelihood
#'   function. For example, [`bgbb.rf.matrix.LL`] requires rf.matrix;
#'   [`pnbd.cbs.LL`] requires cal.cbs and hardie (defaults to TRUE); and [`bgnbd.cbs.LL`] requires
#'   cal.cbs.
#' @param n.divs 	integer representing how fine-grained the contour plot is. A
#'   higher value will produce a higher resolution plot with smoother contour
#'   lines, but will take longer to plot. n.divs also affects the boundaries of
#'   the contour plot; see details.
#' @param new.dev   if TRUE, makes a new window for each contour plot.
#' @param num.contour.lines     number of contour lines to plot in the window.
#' @param zoom.percent  determines boundaries of contour plot. See details.
#' @param allow.neg.params  if FALSE, the contour plot will not include negative
#'   values (see details). This should be set to false for the BG/BB and
#'   Pareto/NBD models.
#' @param param.names   a vector containing parameter names.
#' @seealso [`dc.PlotLogLikelihoodContours`]
#' @examples 
#' # **Examples for BG/BB model:
#' data(donationsSummary)
#' rf.matrix <- donationsSummary$rf.matrix
#' 
#' # starting-point parameters
#' bgbb.startingparams <- c(1, 1, 0.5, 3)
#' # estimated parameters
#' bgbb.est.params <- bgbb.EstimateParameters(rf.matrix, bgbb.startingparams)
#' 
#' # set up parameter names for a more descriptive result
#' bgbb.param.names <- c("alpha", "beta", "gamma", "delta")
#' 
#' # plot a log-likelihood contour of alpha and beta, the unobserved
#' # parameters for the beta-Bernoulli transaction process of the BG/BB.
#' # Note that allow.neg.params has been set to false as BG/BB parameters
#' # cannot be negative:
#' dc.PlotLogLikelihoodContour(bgbb.rf.matrix.LL, 
#'                             c("vary", "vary", "fix", "fix"), 
#'                             bgbb.est.params, 
#'                             rf.matrix = rf.matrix, 
#'                             n.divs = 15, 
#'                             num.contour.lines = 15, 
#'                             zoom.percent = 0.2, 
#'                             allow.neg.params = FALSE, 
#'                             param.names = bgbb.param.names)
#' 
#' # plot a log-likelihood contour of gamma and delta, the unobserved
#' # parameters for the beta-geometric dropout process of the BG/BB.
#' # Note that allow.neg.params has been set to false as BG/BB parameters
#' # cannot be negative:
#' dc.PlotLogLikelihoodContour(bgbb.rf.matrix.LL, 
#'                             c("fix", "fix", "vary", "vary"), 
#'                             bgbb.est.params, 
#'                             rf.matrix = rf.matrix, 
#'                             n.divs = 15, 
#'                             num.contour.lines = 15, 
#'                             zoom.percent = 0.2, 
#'                             allow.neg.params = FALSE, 
#'                             param.names = bgbb.param.names)
#' 
#' # **Example for Pareto/NBD model:
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' 
#' # Speed up calculations:
#' cbs <- dc.compress.cbs(cbs)
#' 
#' # parameters estimated using pnbd.EstimateParameters
#' pnbd.est.params <- cdnowSummary$est.params
#' 
#' # set up parameter names for a more descriptive result
#' pnbd.param.names <- c("r", "alpha", "s", "beta")
#' 
#' # plot a log-likelihood contour of r and s, the shape parameters
#' # of the transaction and dropout process models (respectively).
#' # Note that allow.neg.params has been set to false as Pareto/NBD
#' # parameters cannot be negative:
#' dc.PlotLogLikelihoodContour(pnbd.cbs.LL, 
#'                             c("vary", "fix", "vary", "fix"), 
#'                             pnbd.est.params, 
#'                             cal.cbs = cbs, 
#'                             hardie = TRUE,
#'                             n.divs = 20, 
#'                             num.contour.lines = 20, 
#'                             zoom.percent = 0.1, 
#'                             allow.neg.params = FALSE, 
#'                             param.names = pnbd.param.names)
#' 
#' # **Example for BG/NBD model:
#' data(cdnowSummary)
#' cbs <- cdnowSummary$cbs
#' 
#' # parameters estimated using bgnbd.EstimateParameters
#' bgnbd.est.params <- cdnowSummary$est.params
#' 
#' # set up parameter names for a more descriptive result
#' bgnbd.param.names <- c("r", "alpha", "s", "beta")
#' 
#' # plot a log-likelihood contour of r and s, the shape parameters
#' # of the transaction and dropout process models (respectively).
#' # Note that allow.neg.params has been set to false as BG/NBD
#' # parameters cannot be negative:
#' dc.PlotLogLikelihoodContour(bgnbd.cbs.LL, 
#'                             c("vary", "fix", "vary", "fix"), 
#'                             bgnbd.est.params, 
#'                             cal.cbs = cbs, 
#'                             n.divs = 20, 
#'                             num.contour.lines = 20, 
#'                             zoom.percent = 0.1, 
#'                             allow.neg.params = FALSE, 
#'                             param.names = bgnbd.param.names)
#' @md
dc.PlotLogLikelihoodContour <- function(loglikelihood.fcn, 
                                        vary.or.fix.param, 
                                        predicted.params, 
                                        ..., 
                                        n.divs = 3, 
                                        new.dev = FALSE, 
                                        num.contour.lines = 10, 
                                        zoom.percent = 0.9, 
                                        allow.neg.params = FALSE, 
                                        param.names = c("param 1", "param 2", "param 3", "param 4")) {
    if (new.dev) {
        dev.new()
    }
    idx.par.vary <- which(vary.or.fix.param == "vary")
    
    if (length(idx.par.vary) != 2) {
        stop("vary.or.fix.param must have exactly two elements: \"vary\" ")
    }
    
    values.par.vary <- predicted.params[idx.par.vary]
    v1 <- values.par.vary[1]
    v2 <- values.par.vary[2]
    par1.ticks <- c(v1 - (n.divs:1) * zoom.percent, v1, v1 + (1:n.divs) * zoom.percent)
    par2.ticks <- c(v2 - (n.divs:1) * zoom.percent, v2, v2 + (1:n.divs) * zoom.percent)
    
    param.names.vary <- param.names[idx.par.vary]
    
    if (!allow.neg.params) {
        par1.ticks <- par1.ticks[par1.ticks > 0]
        par2.ticks <- par2.ticks[par2.ticks > 0]
    }
    n.par1.ticks = length(par1.ticks)
    n.par2.ticks = length(par2.ticks)
    
    ll <- sapply(0:(n.par1.ticks * n.par2.ticks - 1), function(e) {
        i <- (e%%n.par1.ticks) + 1
        j <- (e%/%n.par1.ticks) + 1
        current.params <- predicted.params
        current.params[idx.par.vary] <- c(par1.ticks[i], par2.ticks[j])
        loglikelihood.fcn(current.params, ...)
    })
    
    loglikelihood.contours <- matrix(ll, nrow = n.par1.ticks, ncol = n.par2.ticks)
    
    if (FALSE) {
        for (ii in 1:n.par1.ticks) {
            for (jj in 1:n.par2.ticks) {
                current.params <- predicted.params
                current.params[idx.par.vary] <- c(par1.ticks[ii], par2.ticks[jj])
                loglikelihood.contours[ii, jj] <- loglikelihood.fcn(current.params, 
                  ...)
                ## cat('finished', (ii-1)*2*n.divs+jj, 'of', 4*n.divs*n.divs, fill=TRUE)
            }
        }
    }
    contour.plot <- contour(x = par1.ticks, y = par2.ticks, z = loglikelihood.contours, 
        nlevels = num.contour.lines)
    # label.varying.params <- paste(idx.par.vary, collapse=', ')
    
    contour.plot.main.label <- paste("Log-likelihood contour of", param.names.vary[1], 
        "and", param.names.vary[2])
    abline(v = values.par.vary[1], h = values.par.vary[2], col = "red")
    
    title(main = contour.plot.main.label, 
          xlab = param.names.vary[1], 
          ylab = param.names.vary[2])
}

#' Read Lines
#' 
#' Given a .csv file that throws errors when read in by the usual read.csv and read.table methods, 
#' loops through the file line-by-line and picks out the customer, date, and sales (optional) 
#' transaction data to return an event log.
#' 
#' Once this function has been run, you may need to convert the date column to Date objects for 
#' the event log to work with other functions in this library. See the as.Date function in the 
#' `base` R package for more details.
#' 
#' @param csv.filename The name of the comma-delimited file to be read. This file must contain headers.
#' @param cust.idx The index of the customer ID column in the comma-delimited file.
#' @param date.idx The index of the date column in the comma-delimited file.
#' @param sales.idx The index of the sales column in the comma-delimited file.
#' 
#' @examples 
#' # Create event log from file "cdnowElog.csv", which has
#' # customer IDs in the second column, dates in the third column, and
#' # sales numbers in the fifth column.
#' elog <- dc.ReadLines(system.file("data/cdnowElog.csv", package="BTYD"),2,3,5)
#' 
#' # convert date column to date objects, as required by some other functions
#' elog$date <- as.Date(elog$date, "$Y%m%d")
#' @md
dc.ReadLines <- function(csv.filename, 
                         cust.idx, 
                         date.idx, 
                         sales.idx = -1) {
    dc.WriteLine("Started reading file. Progress:")
    elog.file <- file(csv.filename, open = "r")
    elog.lines <- readLines(elog.file)
    n.lines <- length(elog.lines) - 1
    cust <- rep("", n.lines)
    date <- rep("", n.lines)
    if (sales.idx != -1) {
        sales <- rep(0, n.lines)
    }
    
    for (ii in 2:(n.lines + 1)) {
        ## splitting each line by commas
        split.string <- strsplit(elog.lines[ii], ",")
        ## assigning the comma delimited values to our vector
        this.cust <- split.string[[1]][cust.idx]
        this.date <- split.string[[1]][date.idx]
        if (is.na(this.cust) | is.na(this.date)) {
            next
        }
        cust[ii - 1] <- this.cust
        date[ii - 1] <- this.date
        if (sales.idx != -1) {
            sales[ii - 1] <- split.string[[1]][sales.idx]
        }
        ## Progress bar:
        if (ii%%1000 == 0) {
            dc.WriteLine(ii, "/", n.lines)
        }
    }
    
    elog <- cbind(cust, date)
    elog.colnames <- c("cust", "date")
    
    if (sales.idx != -1) {
        elog <- cbind(elog, sales)
        elog.colnames <- c(elog.colnames, "sales")
    }
    
    elog <- data.frame(elog, stringsAsFactors = FALSE)
    colnames(elog) <- elog.colnames
    
    if (sales.idx != -1) {
        elog$sales <- as.numeric(elog$sales)
    }
    close(elog.file)
    dc.WriteLine("File successfully read.")
    return(elog)
}

#' Check model params
#' 
#' Check model parameters for correctness.
#' 
#' @param printnames Names to print parameter errors.
#' @param params Model parameters.
#' @param func Function calling dc.check.model.params.
#' @return Stops the program if there is something wrong with the parameters.
dc.check.model.params <- function(printnames, 
                                  params, 
                                  func) {
    if (length(params) != length(printnames)) {
        stop("Error in ", func, ": Incorrect number of parameters; there should be ", 
            length(printnames), ".", call. = FALSE)
    }
    if (!is.numeric(params)) {
        stop("Error in ", func, ": parameters must be numeric, but are of class ", 
            class(params), call. = FALSE)
    }
    if (any(params < 0)) {
        stop("Error in ", func, ": All parameters must be positive. Negative parameters: ", 
            paste(printnames[params < 0], collapse = ", "), call. = FALSE)
    }
}

#' Cumulative to Incremental
#' 
#' Converts a vector of cumulative transactions to a vector of incremental transactions.
#' 
#' @param cu A vector containing cumulative transactions over time.
#' @return Vector of incremental transactions.
dc.CumulativeToIncremental <- function(cu) {
    inc <- cu - c(0, cu)[-(length(cu) + 1)]
    return(inc)
} 
