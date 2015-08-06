source("setup.R")

params <- bgnbd.EstimateParameters(cal.cbs)
params
LL <- bgnbd.cbs.LL(params, cal.cbs)
p.matrix <- c(params, LL);
for (i in 1:2){
  params <- bgnbd.EstimateParameters(cal.cbs, params);
  LL <- bgnbd.cbs.LL(params, cal.cbs);
  p.matrix.row <- c(params, LL);
  p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("r", "alpha", "a", "b", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;

bgnbd.PlotTransactionRateHeterogeneity(params)
bgnbd.PlotDropoutRateHeterogeneity(params)
bgnbd.Expectation(params, t=52);

# Expected behavior from a particular customer
custName <- sample(cal.cbs[,1],1)
custName

cal.cbs[custName,]

x<-cal.cbs[custName,"x"]
t.x <-cal.cbs[custName,"t.x"]
T.cal <- cal.cbs[custName,"T.cal"]

bgnbd.ConditionalExpectedTransactions(params, T.star = 52, x, t.x, T.cal)
bgnbd.PAlive(params, x, t.x, T.cal)


for (i in seq(10, 25, 5)){
  cond.expectation <- bgnbd.ConditionalExpectedTransactions(
    params, T.star = 52, x = i,
    t.x = 20, T.cal = 39)
  cat ("x:",i,"\t Expectation:",cond.expectation, fill = TRUE)
}

censor <- 7
bgnbd.PlotFrequencyInCalibration(params, cal.cbs, censor)

# Verify in Holdout Period
x.star <- hold.cbs[,"x.star"]
comp <- bgnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star=52, cal.cbs, x.star, censor)
rownames(comp) <- c("act", "exp", "bin")
comp

bgnbd.PlotRecVsConditionalExpectedFrequency(params, cal.cbs, T.star=52, x.star)

# Plot Actual V/s Expected Transactions on a weekly basis
tot.cbt <- dc.CreateFreqCBT(elog)
head(tot.cbt)

# ...Completed Freq CBT
d.track.data <- rep(0, 7 * 105)
origin <- as.Date("2013-01-01")
for (i in colnames(tot.cbt)){
  date.index <- difftime(as.Date(i), origin) + 1;
  d.track.data[date.index] <- sum(tot.cbt[,i]);
}
w.track.data <- rep(0, 105)
for (j in 1:105){
  w.track.data[j] <- sum(d.track.data[(j*7-6):(j*7)])
}

T.cal <- cal.cbs[,"T.cal"]
T.tot <- 105
n.periods.final <- 105
inc.tracking <- bgnbd.PlotTrackingInc(params, T.cal,
                                     T.tot, w.track.data,
                                     n.periods.final)
inc.tracking[,20:25]

cum.tracking.data <- cumsum(w.track.data)
cum.tracking <- bgnbd.PlotTrackingCum(params, T.cal,
                                     T.tot, cum.tracking.data,
                                     n.periods.final)
cum.tracking[,20:25]



