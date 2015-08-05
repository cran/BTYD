# Import the data and create the customer-by-time-matrix
# Get ready to model

##Install packages
toInstallCandidates <- c("ggplot2", "reshape2", "plyr", "lubridate", "devtools", "gsl")

# check if pkgs are already present

toInstall <- toInstallCandidates[!toInstallCandidates%in%library()$results[,1]] 
if(length(toInstall)!=0)
{install.packages(toInstall, repos = "http://cran.r-project.org")}

# load pkgs
lapply(toInstallCandidates, library, character.only = TRUE)


# Cleanup existing installation of BTYD
if ("package:BTYD" %in% search()) { detach("package:BTYD", unload=TRUE) }
if ("BTYD" %in% rownames(installed.packages())) { remove.packages("BTYD") }

#Change local to 1 if using the local version of BTYD instead of CRAN version
local <-1

if (local){
  devtools::install_github("jamespaul007/BTYD", ref="pnbd-fix")
} else{
  install.packages("BTYD")
}
library(BTYD)

elogFile <- "TxData.csv"
elog <- dc.ReadLines(elogFile, cust.idx = 1,
                     date.idx = 2, sales.idx = 3)
head(elog)
elog$date <- as.Date(elog$date, "%m/%d/%Y");
head(elog)
summary(elog)

# Optional - Sample the data
#class(elog)
#elog <- elog[sample(nrow(elog), 3000), ]
#str(elog)

max(elog$date);
min(elog$date);
qplot(date, sales, data = elog)
xtabs( sales ~ cust, data =elog)

# Approach #1
# Merge transactions on the same day
elog <- dc.MergeTransactionsOnSameDate(elog)

# Split into calibration and holdout period
end.of.cal.period <- as.Date("2013-12-31")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]

# Information about the customer
split.data <- dc.SplitUpElogForRepeatTrans(elog.cal)
clean.elog <- split.data$repeat.trans.elog

## Create customer-by-time matrix
freq.cbt <- dc.CreateFreqCBT(clean.elog)
tot.cbt <- dc.CreateFreqCBT(elog)
cal.cbt <- dc.MergeCustomers(tot.cbt, freq.cbt)
birth.periods <- split.data$cust.data$birth.per
last.dates <- split.data$cust.data$last.date
cal.cbs.dates <- data.frame(birth.periods, last.dates,
                            end.of.cal.period)
cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.cbt, cal.cbs.dates,
                                      per="week")
# Approach #2
# Alternative Method - ElogToCbsCbt

T.cal <- as.Date("2013-12-31")
simData <- dc.ElogToCbsCbt(elog, per="week", T.cal)
cal.cbs <- simData$cal$cbs
cal.cbt <- simData$cal$cbt

# for Holdout period
hold.cbs <- simData$holdout$cbs
head(hold.cbs)
hold.cbt <- simData$holdout$cbt

