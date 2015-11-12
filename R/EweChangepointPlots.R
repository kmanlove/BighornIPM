EweChangepointsPlot <- function(compd.data){
  aso <- subset(compd.data, Pop == "Asotin", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  aso <- aso[complete.cases(aso), ]
  bb <- subset(compd.data, Pop == "BlackButte", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  bb <- bb[complete.cases(bb), ]
  bc <- subset(compd.data, Pop == "BigCanyon", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  bc <- bc[complete.cases(bc), ]
  im <- subset(compd.data, Pop == "Imnaha", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  im <- im[complete.cases(im), ]
  lo <- subset(compd.data, Pop == "Lostine", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  lo <- lo[complete.cases(lo), ]
  los <- subset(compd.data, Pop == "Lostine", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  los <- los[complete.cases(los), ]
  lhc <- subset(compd.data, Pop == "LowerHells", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  lhc <- lhc[complete.cases(lhc), ]
  mv <- subset(compd.data, Pop == "MountainView" & year >= 1994, select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  mv <- mv[complete.cases(mv), ]
  mu <- subset(compd.data, Pop == "MuirCreek", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  mu <- mu[complete.cases(mu), ]
  my <- subset(compd.data, Pop == "MyersCreek", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  my <- my[complete.cases(my), ]
  rb <- subset(compd.data, Pop == "Redbird" & year >= 1992, select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  rb <- rb[complete.cases(rb), ]
  uhc <- subset(compd.data, Pop == "UHCOR")
  uhc <- subset(compd.data, Pop == "UHCOR", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  uhc <- uhc[complete.cases(uhc), ]
  we <- subset(compd.data, Pop == "Wenaha", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est"))
  we <- we[complete.cases(we), ]
  
  # time and time.2 are both square matrices: n.pop x n.year
  aso.years <- aso$year - 2011
  aso.counts <- aso$Ewes
  aso.r_est <- aso$r_est
  bb.years <- bb$year - 1995
  bb.counts <- bb$Ewes
  bb.r_est <- bb$r_est
  bc.years <- bc$year - 2000
  bc.counts <- bc$Ewes
  bc.r_est <- bc$r_est
  im.years <- im$year - 2000
  im.counts <- im$Ewes
  im.r_est <- im$r_est
  los.years <- los$year - 1986
  los.counts <- los$Ewes
  los.r_est <- los$r_est
  lhc.years <- lhc$year - 1995
  lhc.counts <- lhc$Ewes
  lhc.r_est <- lhc$r_est
  mv.years <- mv$year - 1995
  mv.counts <- mv$Ewes
  mv.r_est <- mv$r_est
  mu.years <- mu$year - 2000
  mu.counts <- mu$Ewes
  mu.r_est <- mu$r_est
  my.years <- my$year - 2002
  my.counts <- my$Ewes
  my.r_est <- my$r_est
  rb.years <- rb$year - 1995
  rb.counts <- rb$Ewes
  rb.r_est <- rb$r_est
  uhc.years <- uhc$year - 2002
  uhc.counts <- uhc$Ewes
  uhc.r_est <- uhc$r_est
  we.years <- we$year - 1995
  we.counts <- we$Ewes
  we.r_est <- we$r_est
  
   par(oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 2))
   layout.mat <- rbind(c(1, 1, 2), c(1, 1, 3))
   layout(layout.mat)
   plot(aso.counts ~ aso.years, xlab = "Years from first observed pneumonia", 
        ylab = "Ewe count", type = "l", ylim = c(0, 100), xlim = c(-20, 20), 
        col = rgb(.2, .2, .2, alpha = .8))
   lines(bb.counts ~ bb.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(bc.counts ~ bc.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(im.counts ~ im.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(lhc.counts ~ lhc.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(los.counts ~ los.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(mv.counts ~ mv.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(mu.counts ~ mu.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(my.counts ~ my.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(rb.counts ~ rb.years, col = rgb(.2, .2, .2, alpha = .8))
   lines(we.counts ~ we.years, col = rgb(.2, .2, .2, alpha = .8))
   abline(v = 0, lty = 2, col = "grey60")
   polygon(x = c(-3.14, -3.14, -1.01, -1.01), 
           y = c(0, 100, 100, 0), 
           col = rgb(1, .1, .1, alpha = .2), 
           border = rgb(1, .1, .1, alpha = .2))
   polygon(x = c(6.74, 6.74, 8.59, 8.59), 
           y = c(0, 100, 100, 0), 
           col = rgb(1, .1, .1, alpha = .2), 
           border = rgb(1, .1, .1, alpha = .2))
   
   plot(im.counts[im.years >= 0] ~ im.years[im.years >= 0], type = "l", 
        ylim = c(0, 100), xlim = c(0, 15),
        xlab = "Years post-invasion", ylab = "Ewe count",
        main = "Second changepoint", col = "darkgreen")
   text(x = 3, y = 90, "Imnaha", cex = .8, col = "darkgreen")
   lines(los.counts[los.years >= 0] ~ los.years[los.years >= 0], col = "darkblue")
   text(x = 3, y = 15, "Lostine", cex = .8, col = "darkblue")
   lines(rb.counts[rb.years >= 0] ~ rb.years[rb.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   text(x = 6, y = 50, "Redbird", cex = .8, col = "black")
   polygon(x = c(6.74, 6.74, 8.59, 8.59), 
           y = c(0, 100, 100, 0), 
           col = rgb(1, .1, .1, alpha = .2), 
           border = rgb(1, .1, .1, alpha = .2))
   
   plot(aso.counts[aso.years >= 0] ~ aso.years[aso.years >= 0], 
        type = "l", ylim = c(0, 100), xlim = c(0, 15),
        xlab = "Years post-invasion", ylab = "Ewe count",
        main = "No second changepoint", col = rgb(.2, .2, .2, alpha = .8))
   lines(bb.counts[bb.years >= 0] ~ bb.years[bb.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(bc.counts[bc.years >= 0] ~ bc.years[bc.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(lhc.counts[los.years >= 0] ~ lhc.years[lhc.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(mv.counts[mv.years >= 0] ~ mv.years[mv.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(mu.counts[mu.years >= 0] ~ mu.years[mu.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(my.counts[my.years >= 0] ~ my.years[my.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   lines(we.counts[we.years >= 0] ~ we.years[we.years >= 0], col = rgb(.2, .2, .2, alpha = .8))
   
 }