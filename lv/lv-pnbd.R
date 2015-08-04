source("setup.R")

#cal.cbs <- as.data.frame(cal.cbs)
params <- pnbd.EstimateParameters(cal.cbs)
params

LL <- pnbd.cbs.LL(params, cal.cbs);
LL

p.matrix <- c(params, LL);
for (i in 1:2){
  params <- pnbd.EstimateParameters(cal.cbs, params);
  LL <- pnbd.cbs.LL(params, cal.cbs);
  p.matrix.row <- c(params, LL);
  p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("r", "alpha", "s", "beta", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;


pnbd.PlotTransactionRateHeterogeneity(params)
pnbd.PlotDropoutRateHeterogeneity(params)

# Individual level estimations

pnbd.Expectation(params, t=52)

# Expected behavior from a particular customer
custName <- sample(cal.cbs[,1],1)
custName

cal.cbs[custName,]

x<-cal.cbs[custName,"x"]
t.x <-cal.cbs[custName,"t.x"]
T.cal <- cal.cbs[custName,"T.cal"]

pnbd.ConditionalExpectedTransactions(params, T.star = 52, x, t.x, T.cal)
pnbd.PAlive(params, x, t.x, T.cal)

# To visualize the distribution of P(Alive) across customers:
p.alives <- pnbd.PAlive(params, cal.cbs[,"x"], cal.cbs[,"t.x"], cal.cbs[,"T.cal"])

ggplot(as.data.frame(p.alives),aes(x=p.alives))+
  geom_histogram(colour="grey",fill="orange")+
  ylab("Number of Customers")+
  xlab("Probability Customer is 'Live'")+
  theme_minimal()

# Goodness of Fit
censor <- 10
pnbd.PlotFrequencyInCalibration(params, cal.cbs, censor)

# Verify in Holdout Period
x.star <- hold.cbs[,"x.star"]
comp <- pnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star=52, cal.cbs, x.star, censor)
rownames(comp) <- c("act", "exp", "bin")
comp

pnbd.PlotRecVsConditionalExpectedFrequency(params, cal.cbs, T.star=52, x.star)

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
inc.tracking <- pnbd.PlotTrackingInc(params, T.cal,
                                     T.tot, w.track.data,
                                     n.periods.final)
inc.tracking[,20:25]

cum.tracking.data <- cumsum(w.track.data)
cum.tracking <- pnbd.PlotTrackingCum(params, T.cal,
                                     T.tot, cum.tracking.data,
                                     n.periods.final)
cum.tracking[,20:25]




