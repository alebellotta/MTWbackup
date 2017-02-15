####  JUMP COUNTER IN A TIME SERIES ####
#Data
ticker = "AAPL"
stockprice = as.numeric(get.hist.quote(instrument = ticker, quote = "Close"))

jcount = NULL
jcount[1] = 0
for (i in 2:length(stockprice)){
  if (stockprice[i] > stockprice[i-1]*1.025 & stockprice[i] < stockprice[i-1]*1.05){
    jcount[i] = 1
  } else if (stockprice[i] >= stockprice[i-1]*1.5 & stockprice[i] < stockprice[i-1]*1.075 ){
    jcount[i] = 2
  } else if (stockprice[i] >= stockprice[i-1]*1.075 & stockprice[i] < stockprice[i-1]*1.10 ){
    jcount[i] = 3
  } else if (stockprice[i] >= stockprice[i-1]*1.10 & stockprice[i] < stockprice[i-1]*1.125 ){
    jcount[i] = 4
  } else if (stockprice[i] <= stockprice[i-1]*0.975 & stockprice[i] > stockprice[i-1]*0.95){
    jcount[i] = -1
  } else if (stockprice[i] <= stockprice[i-1]*0.95 & stockprice[i] > stockprice[i-1]*0.925){
    jcount[i] = -2
  } else if (stockprice[i] <= stockprice[i-1]*0.925 & stockprice[i] > stockprice[i-1]*0.875){
    jcount[i] = -3
  } else if (stockprice[i] <= stockprice[i-1]*0.875 & stockprice[i] > stockprice[i-1]*0.85){
    jcount[i] = -4
  } else if (stockprice[i] >= stockprice[i-1]*1.125 & stockprice[i] < stockprice[i-1]*1.15){
    jcount[i] = 5
  } else if (stockprice[i] >= stockprice[i-1]*1.15 & stockprice[i] < stockprice[i-1]*1.175 ){
    jcount[i] = 6
  } else if (stockprice[i] <= stockprice[i-1]*0.85 & stockprice[i] > stockprice[i-1]*0.825){
    jcount[i] = -5
  } else if (stockprice[i] <= stockprice[i-1]*0.825 & stockprice[i] > stockprice[i-1]*0.8){
    jcount[i] = -6
  } else {
    jcount[i] = 0
  }
}

plot(sort(jcount))
muJ = mean(jcount)
logjcount = pmax(log(jcount), 0)
logjcount[is.na(logjcount)] = 0
vJ = var(logjcount)
