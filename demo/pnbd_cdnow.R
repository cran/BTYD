## Authors: Lukasz Dziurzynski, Edward Wadsworth

data(cdnowSummary)

## Get the calibration period customer-by-sufficient-statistic matrix from the cdnow data:
cbs <- cdnowSummary$cbs

## Estimate parameters for the Pareto/NBD model from the CBS:
par.start <- c(0.5, 1, 0.5, 1)
params <- pnbd.EstimateParameters(cbs, par.start)
params

## Check log-likelihood of the params:
pnbd.cbs.LL(params, cbs)

## Plot the comparison of actual and expected calibration period frequencies:
pnbd.PlotFrequencyInCalibration(params, cbs, censor=7, plotZero=TRUE)

T.star <- 39                # Length of holdout period
x.star <- cbs[,"x.star"]    # Transactions made by each customer in the holdout period

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period frequencies:
pnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star, cbs, x.star, censor=7)

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period recencies:
pnbd.PlotRecVsConditionalExpectedFrequency(params, cbs, T.star, x.star)

cum.trans <- cdnowSummary$cu.tracking                           # cumulative weekly transactions
inc.trans <- pnbd.CumulativeToIncremental(cum.trans)     # incremental weekly transactions
T.cal <- max(cbs[,"T.cal"])                                     # length of calibration period
T.tot <- 78                                                     # end of holdout period
time.periods <- 78                                              # total number of time periods

## Plot the comparison of actual and expected total cumulative transactions across
## both the calibration and holdout periods:
pnbd.PlotTrackingCum(params, cbs, T.cal, T.tot, cum.trans, time.periods)

## Plot the comparison of actual and expected total incremental transactions across
## both the calibration and holdout periods:
pnbd.PlotTrackingInc(params, cbs, T.cal, T.tot, inc.trans, time.periods)
