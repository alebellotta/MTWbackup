require(quantmod)
require(PerformanceAnalytics)

#get DJIA since 1896 from St. Louis Fed Fred
getSymbols("DJIA",src="FRED")

#do monthly to shorten the lengthy calculation
#on my old computer
DJIA <- to.monthly(DJIA)[,4]
retDJIA <- ROC(DJIA,n=1,type="discrete")

#paper HURST EXPONENT AND FINANCIAL MARKET PREDICTABILITY by
#Bo Qian and Khaled Rasheed use 1024 days or 4 years
#or 208 weeks assuming 52 weeks/year
#I first tried HurstIndex from PerformanceAnalytics
hurst.48 <- apply.rolling(retDJIA, FUN="HurstIndex", width=48)
chart.TimeSeries(hurst.48,colorset=c("cadetblue"), main = "Hurst Index (48 Month) for DJIA 1896-2011")
#my results are very different from the research paper whether I use
#daily, weekly, or monthly
#so I must be doing something wrong
hurst.avg <- apply.rolling(retDJIA, FUN="HurstIndex", width=12)
hurst.avg <- runMean(hurst.avg,n=48)
hurst <- merge(hurst.48,hurst.avg)
colnames(hurst) <- c("Hurst 48 Month","Hurst 4y Avg of 12 Month")
chart.TimeSeries(hurst,colorset=c("cadetblue","darkolivegreen3"),
	#legend.loc="bottomleft",
	main = "Hurst Index Comparison for DJIA
	1896-2011")
abline(h=0.5) #this represents no trend or mean reversion


#now let's try a different package
require(FGN)
#get DJIA since 1896 from St. Louis Fed Fred
getSymbols("DJIA",src="FRED")
#will do daily; takes about 2 hours
#on my old computer
DJIA <- DJIA
retDJIA <- ROC(DJIA,n=1)
#about 2 hours for a result on my old computer
hurstK <- apply.rolling(retDJIA, FUN="HurstK", width=1024)
chart.TimeSeries(hurstK,colorset=c("cadetblue"),
	main = "HurstK Calculation of Hurst Exponent (1024 days) for DJIA
	1896-2011")


#######system-building time
#do monthly just to test more quickly
DJIA <- to.monthly(DJIA)[,4]
index(DJIA) <- as.Date(index(DJIA))
retDJIA<-ROC(DJIA,n=1,type="discrete")
index(retDJIA) <- as.Date(index(retDJIA))
hurstKmonthly <- apply.rolling(retDJIA, FUN="HurstK", width = 12)
colnames(hurstKmonthly) <- "HurstK.monthly"
index(hurstKmonthly) <- as.Date(index(hurstKmonthly))
serialcorr <- runCor(cbind(coredata(retDJIA)),cbind(index(retDJIA)),n=12)
serialcorr <- as.xts(serialcorr,order.by=index(retDJIA))
autoreg <- runCor(retDJIA,lag(retDJIA,k=1),n=12)
colnames(serialcorr) <- "SerialCorrelation.monthly"
colnames(autoreg) <- "AutoRegression.monthly"
#check for correlation of potential signals
chart.Correlation(merge(hurstKmonthly,serialcorr,autoreg))
#try a sum to enter in strong trends
signal <- hurstKmonthly+serialcorr+autoreg
colnames(signal) <- "HurstCorrelationSum"
chart.TimeSeries(signal)
signal <- lag(signal,k=1)
retSys <- ifelse(signal > 0.1, 1, 0) * retDJIA
charts.PerformanceSummary(retSys,ylog=TRUE)
#ok performance but let's see if we can enter
#only strong trends up and reduce the drawdown
signalUpTrend <- runMean(hurstKmonthly+serialcorr+autoreg,n=6) + (DJIA/runMean(DJIA,n=12)-1)*10
chart.TimeSeries(signalUpTrend)
signalUpTrend <- lag(signalUpTrend,k=1)
retSys <- merge(retSys,ifelse(signalUpTrend > 1, 1, 0) * retDJIA,retDJIA)
colnames(retSys) <- c("DJIA Hurst System","DJIA HurstUp System",
	"DJIA")
charts.PerformanceSummary(retSys,ylog=TRUE,cex.legend=1.25,
	colorset=c("cadetblue","darkolivegreen3","gray70"))

#now let's take it out of sample to see how it works
getSymbols("^N225",from="1980-01-01",to=format(Sys.Date(),"%Y-%m-%d"))
N225 <- to.monthly(N225)[,4]
index(N225) <- as.Date(index(N225))
retN225<-ROC(N225,n=1,type="discrete")
index(retN225) <- as.Date(index(retN225))
hurstKmonthly <- apply.rolling(retN225, FUN="HurstK", width = 12)
colnames(hurstKmonthly) <- "HurstK.monthly"
index(hurstKmonthly) <- as.Date(index(hurstKmonthly))
serialcorr <- runCor(cbind(coredata(retN225)),cbind(index(retN225)),n=12)
serialcorr <- as.xts(serialcorr,order.by=index(retN225))
autoreg <- runCor(retN225,lag(retN225,k=1),n=12)
colnames(serialcorr) <- "SerialCorrelation.monthly"
colnames(autoreg) <- "AutoRegression.monthly"
signalUpTrend <- runMean(hurstKmonthly+serialcorr+autoreg,n=6) + (N225/runMean(N225,n=12)-1)*10
chart.TimeSeries(signalUpTrend)
signalUpTrend <- lag(signalUpTrend,k=1)
retSys <- merge(ifelse(signalUpTrend > 1, 1, 0) * retN225,retN225)
colnames(retSys) <- c("Nikkei 225 HurstUp System","Nikkei 225")
charts.PerformanceSummary(retSys,ylog=TRUE,cex.legend=1.25,
	colorset=c("cadetblue","darkolivegreen3"))
###########################