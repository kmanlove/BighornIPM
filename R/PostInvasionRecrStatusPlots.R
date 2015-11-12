PostInvasionRecrStatusPlots <- function(compd.data){
   aso <- subset(compd.data, Pop == "Asotin", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   aso <- aso[complete.cases(aso[ ,1:3]), ]
   bb <- subset(compd.data, Pop == "BlackButte", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   bb <- bb[complete.cases(bb[ , 1:3]), ]
   bc <- subset(compd.data, Pop == "BigCanyon", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   bc <- bc[complete.cases(bc[ , 1:3]), ]
   im <- subset(compd.data, Pop == "Imnaha", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   im <- im[complete.cases(im[ , 1:3]), ]
   los <- subset(compd.data, Pop == "Lostine", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   los <- los[complete.cases(los[ , 1:3]), ]
   lhc <- subset(compd.data, Pop == "LowerHells", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   lhc <- lhc[complete.cases(lhc[ , 1:3]), ]
   mv <- subset(compd.data, Pop == "MountainView" & year >= 1994, select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   mv <- mv[complete.cases(mv[ , 1:3]), ]
   mu <- subset(compd.data, Pop == "MuirCreek", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   mu <- mu[complete.cases(mu[ , 1:3]), ]
   my <- subset(compd.data, Pop == "MyersCreek", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   my <- my[complete.cases(my[ , 1:3]), ]
   rb <- subset(compd.data, Pop == "Redbird" & year >= 1992, select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   rb <- rb[complete.cases(rb[ , 1:3]), ]
   uhc <- subset(compd.data, Pop == "UHCOR", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   uhc <- uhc[complete.cases(uhc[ , 1:3]), ]
   we <- subset(compd.data, Pop == "Wenaha", select = c("year", "PopEst", "Ewes", "Lambs", "SumLambSurv", "r_est", "Recr", "CLASS", "CLASS_SUSP"))
   we <- we[complete.cases(we[ , 1:3]), ]

   aso.years <- aso$year - 2011
   aso.recr <- aso$Recr
   aso.counts <- aso$Ewes
   aso.r_est <- aso$r_est
   bb.years <- bb$year - 1995
   bb.counts <- bb$Ewes
   bb.recr <- bb$Recr
   bb.r_est <- bb$r_est
   bc.years <- bc$year - 2000
   bc.counts <- bc$Ewes
   bc.r_est <- bc$r_est
   bc.recr <- bc$Recr
   im.years <- im$year - 2000
   im.counts <- im$Ewes
   im.recr <- im$Recr
   im.r_est <- im$r_est
   los.years <- los$year - 1986
   los.counts <- los$Ewes
   los.r_est <- los$r_est
   los.recr <- los$Recr
   lhc.years <- lhc$year - 1995
   lhc.counts <- lhc$Ewes
   lhc.r_est <- lhc$r_est
   lhc.recr <- lhc$Recr
   mv.years <- mv$year - 1995
   mv.counts <- mv$Ewes
   mv.r_est <- mv$r_est
   mv.recr <- mv$Recr
   mu.years <- mu$year - 2000
   mu.counts <- mu$Ewes
   mu.r_est <- mu$r_est
   mu.recr <- mu$Recr
   my.years <- my$year - 2002
   my.counts <- my$Ewes
   my.r_est <- my$r_est
   my.recr <- my$Recr
   rb.years <- rb$year - 1995
   rb.counts <- rb$Ewes
   rb.r_est <- rb$r_est
   rb.recr <- rb$Recr
   uhc.years <- uhc$year - 2002
   uhc.counts <- uhc$Ewes
   uhc.r_est <- uhc$r_est
   uhc.recr <- uhc$Recr
   we.years <- we$year - 1995
   we.counts <- we$Ewes
   we.r_est <- we$r_est
   we.recr <- we$Recr
   
   full.r_ests <- c(
     aso.r_est,
     bb.r_est,
     bc.r_est,
     im.r_est,
     lhc.r_est,
     los.r_est, 
     mv.r_est,
     mu.r_est,
     my.r_est,
     rb.r_est,
     uhc.r_est,
     we.r_est
   )
   
   full.counts <- c(
     aso.counts,
     bb.counts,
     bc.counts,
     im.counts,
     lhc.counts,
     los.counts, 
     mv.counts,
     mu.counts,
     my.counts,
     rb.counts,
     uhc.counts,
     we.counts
   )
   
   time <- c(
     aso.years,
     bb.years,
     bc.years,
     im.years,
     lhc.years,
     los.years,
     mv.years,
     mu.years,
     my.years,
     rb.years,
     uhc.years,
     we.years
   )
   
   full.recr <- c(
     as.numeric(as.character(aso.recr)),
     as.numeric(as.character(bb.recr)),
     as.numeric(as.character(bc.recr)),
     as.numeric(as.character(im.recr)),
     as.numeric(as.character(lhc.recr)),
     as.numeric(as.character(los.recr)),
     as.numeric(as.character(mv.recr)),
     as.numeric(as.character(mu.recr)),
     as.numeric(as.character(my.recr)),
     as.numeric(as.character(rb.recr)),
     as.numeric(as.character(uhc.recr)),
     as.numeric(as.character(we.recr))
   )
   
   pop.ind <- c(
     rep("Aso", length(aso.years)),
     rep("BB", length(bb.years)),
     rep("BC", length(bc.years)),
     rep("Imn", length(im.years)),
     rep("LHC", length(lhc.years)),
     rep("Los", length(los.years)),
     rep("MV", length(mv.years)),
     rep("Mu", length(mu.years)),
     rep("My", length(my.years)),
     rep("RB", length(rb.years)),
     rep("UHC", length(uhc.years)),
     rep("Wen", length(we.years))
   )
   
   class.ind <- c(
     as.character(aso$CLASS),
     as.character(bb$CLASS),
     as.character(bc$CLASS),
     as.character(im$CLASS),
     as.character(lhc$CLASS),
     as.character(los$CLASS),
     as.character(mv$CLASS),
     as.character(mu$CLASS),
     as.character(my$CLASS),
     as.character(rb$CLASS),
     as.character(uhc$CLASS),
     as.character(we$CLASS)
   )
   
   bio.yr <- c(aso$year,
               bb$year,
               bc$year,
               im$year,
               lhc$year,
               los$year,
               mv$year,
               mu$year,
               my$year,
               rb$year,
               uhc$year,
               we$year)
     
   full.data <- as.data.frame(cbind(time, full.recr, pop.ind, class.ind, bio.yr))
   names(full.data) <- c("years.post.invasion", "full.recr", "pop.ind", "class.ind", "bio.yr")
   full.data$years.post.invasion <- as.numeric(as.character(full.data$years.post.invasion))
   full.data <- subset(full.data, !(pop.ind %in% c("UHC", "LHC")) & 
                         is.na(full.recr) == F &
                         years.post.invasion <= 19)
#   write.csv(full.data, "./Data/RecrRelativeToInvasion_151109.csv")
   
   post.invasion <- subset(full.data, years.post.invasion >= 0 & years.post.invasion <= 18)
   pre.invasion <- subset(full.data, years.post.invasion <= -1)
   pre.invasion.quants <- quantile(as.numeric(as.character(pre.invasion$full.recr)), c(0.25, 0.75))
   
   N.tab <- table(post.invasion$time)
   #svg("./Plots/RecrVsYearsPostInvasion_09Nov2015.svg", width = 6, height = 4)
   #dev.off()
   
   post.invasion.sm <- subset(post.invasion, is.na(class.ind) == F & class.ind != "" & !(pop.ind %in% c("UHC", "LHC")))
   post.invasion.sm$class.ind <- factor(post.invasion.sm$class.ind)
   post.invasion.class.tab <- table(post.invasion.sm$class.ind, factor(post.invasion.sm$years.post.invasion))
   prop.post.i.class <- prop.table(post.invasion.class.tab[c(2, 1, 4, 5, 3), c(1:2, 12:19, 3:11)], margin = 2)
    
   # Fig 1c
   #svg("./Plots/HealthClassYearsPostInvasion_09Nov2015.svg", width = 6, height = 4)
#   par(oma = c(0, 0, 0, 0))
   par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(5, 5, 1, 1))
   barplot(prop.post.i.class, col = c(rgb(1, .1, .1, alpha = 1), 
                                      rgb(1, .1, .1, alpha = .75), 
                                      rgb(1, .1, .1, alpha = .5),  
                                      rgb(1, .1, .1, alpha = .25), 
                                      "lightblue"),
           xlab = "Years post-invasion", ylab = "Proportion per health class", 
           legend = T, ylim = c(0, 1.25), args.legend = list("top", ncol = 3, bty = "n", cex = .7))
      
      # Fig 1d
      plot(as.numeric(as.character(full.recr)) ~ factor(years.post.invasion), 
          data = post.invasion, 
          col = rgb(1, 0, 0, alpha = .4), 
          border = rgb(1, .4, .4, alpha = .7),
          xlab = "Years post-invasion",
          ylab = "Recruitment",
          staplewex = .5)
      abline(h = .20, lty = 1, col = "black", lwd = 3)
      lines(lowess(as.numeric(as.character(post.invasion$full.recr)) ~ as.numeric(as.character(post.invasion$years.post.invasion)))
            , lty = 1, lwd = 3, col = "red")
      text(x = as.numeric(as.character(names(N.tab))) + 1, y = 70, paste("(", N.tab, ")", sep = ""),
          cex = .7, col = "grey40")
      polygon(x = c(0, 0, 20, 20), y = c(pre.invasion.quants, rev(pre.invasion.quants)), col = rgb(0, 0, 1, alpha = .2),
          border = rgb(0, 0, 1, alpha = .4))

#dev.off()
 }