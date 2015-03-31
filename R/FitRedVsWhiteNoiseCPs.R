FitRedVsWhiteNoiseCPs <- function(RedVsWhiteOut, reps)
{
  struct.alph.05 <- RedVsWhiteOut[[1]]
  struct.alph.2 <- RedVsWhiteOut[[2]]
  struct.alph1 <- RedVsWhiteOut[[3]]
  struct.alph.05.cp <-struct.alph.2.cp <- struct.alph1.cp <- vector("list", reps)
  struct.alph.05.lm <-struct.alph.2.lm <- struct.alph1.lm <- vector("list", reps)
  diff.struct.alph.05 <- diff.struct.alph.2 <- diff.struct.alph1 <- rep(NA, reps)
  
  unstruct.alph.05 <- RedVsWhiteOut[[4]]
  unstruct.alph.2 <- RedVsWhiteOut[[5]]
  unstruct.alph1 <- RedVsWhiteOut[[6]]
  unstruct.alph.05.cp <-unstruct.alph.2.cp <- unstruct.alph1.cp <- vector("list", reps)
  unstruct.alph.05.lm <-unstruct.alph.2.lm <- unstruct.alph1.lm <- vector("list", reps)
  diff.unstruct.alph.05 <- diff.unstruct.alph.2 <- diff.unstruct.alph1 <- rep(NA, reps)

  spillover.alph.05 <- RedVsWhiteOut[[4]]
  spillover.alph.2 <- RedVsWhiteOut[[5]]
  spillover.alph1 <- RedVsWhiteOut[[6]]
  spillover.alph.05.cp <-spillover.alph.2.cp <- spillover.alph1.cp <- vector("list", reps)
  spillover.alph.05.lm <-spillover.alph.2.lm <- spillover.alph1.lm <- vector("list", reps)
  diff.spillover.alph.05 <- diff.spillover.alph.2 <- diff.spillover.alph1 <- rep(NA, reps)
  
  for(i in 1:reps){
    struct.alph.05.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = struct.alph.05[i, 1:59], CI = T, bootstrap.samples = 500)
    struct.alph.05.lm[[i]] <- lm(struct.alph.05[i, 1:59] ~ seq(1:59))
    diff.struct.alph.05[i] <- AIC(struct.alph.05.lm[[i]]) - AIC(struct.alph.05.cp[[i]]) 

    struct.alph.2.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = struct.alph.2[i, 1:59], CI = T, bootstrap.samples = 500)
    struct.alph.2.lm[[i]] <- lm(struct.alph.2[i, 1:59] ~ seq(1:59))
    diff.struct.alph.2[i] <- AIC(struct.alph.2.lm[[i]]) - AIC(struct.alph.2.cp[[i]]) 

    struct.alph1.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = struct.alph1[i, 1:59], CI = T, bootstrap.samples = 500)
    struct.alph1.lm[[i]] <- lm(struct.alph1[i, 1:59] ~ seq(1:59))
    diff.struct.alph1[i] <- AIC(struct.alph1.lm[[i]]) - AIC(struct.alph1.cp[[i]]) 
    
    unstruct.alph.05.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = unstruct.alph.05[i, 1:59], CI = T, bootstrap.samples = 500)
    unstruct.alph.05.lm[[i]] <- lm(unstruct.alph.05[i, 1:59] ~ seq(1:59))
    diff.unstruct.alph.05[i] <- AIC(unstruct.alph.05.lm[[i]]) - AIC(unstruct.alph.05.cp[[i]]) 

    unstruct.alph.2.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = unstruct.alph.2[i, 1:59], CI = T, bootstrap.samples = 500)
    unstruct.alph.2.lm[[i]] <- lm(unstruct.alph.2[i, 1:59] ~ seq(1:59))
    diff.unstruct.alph.2[i] <- AIC(unstruct.alph.2.lm[[i]]) - AIC(unstruct.alph.2.cp[[i]]) 

    unstruct.alph1.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = unstruct.alph1[i, 1:59], CI = T, bootstrap.samples = 500)
    unstruct.alph1.lm[[i]] <- lm(unstruct.alph1[i, 1:59] ~ seq(1:59))
    diff.unstruct.alph1[i] <- AIC(unstruct.alph1.lm[[i]]) - AIC(unstruct.alph1.cp[[i]]) 
    
    spillover.alph.05.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = spillover.alph.05[i, 1:59], CI = T, bootstrap.samples = 500)
    spillover.alph.05.lm[[i]] <- lm(spillover.alph.05[i, 1:59] ~ seq(1:59))
    diff.spillover.alph.05[i] <- AIC(spillover.alph.05.lm[[i]]) - AIC(spillover.alph.05.cp[[i]]) 
    
    spillover.alph.2.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = spillover.alph.2[i, 1:59], CI = T, bootstrap.samples = 500)
    spillover.alph.2.lm[[i]] <- lm(spillover.alph.2[i, 1:59] ~ seq(1:59))
    diff.spillover.alph.2[i] <- AIC(spillover.alph.2.lm[[i]]) - AIC(spillover.alph.2.cp[[i]]) 
    
    spillover.alph1.cp[[i]] <- piecewise.linear(x = as.numeric(seq(1,59)), y = spillover.alph1[i, 1:59], CI = T, bootstrap.samples = 500)
    spillover.alph1.lm[[i]] <- lm(spillover.alph1[i, 1:59] ~ seq(1:59))
    diff.spillover.alph1[i] <- AIC(spillover.alph1.lm[[i]]) - AIC(spillover.alph1.cp[[i]]) 

    print(i)
  }

  breaks.in <- hist(c(diff.struct.alph.05, diff.struct.alph.2, diff.struct.alph1, 
                      diff.unstruct.alph.05, diff.unstruct.alph.2, diff.unstruct.alph1,
                      diff.spillover.alph.05, diff.spillover.alph.2, diff.spillover.alph1
                      ), breaks = 15)$breaks
    
  par(mfrow = c(3, 3))
  hist(diff.struct.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.unstruct.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.spillover.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.struct.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.unstruct.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.spillover.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.struct.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.unstruct.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  hist(diff.spillover.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
 
  par(mfrow = c(3, 3))

  plot(y = c(1, 1), x = c(struct.alph.05.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = struct.alph.05.cp[[i]]$intervals[1 ,1], x1 = struct.alph.05.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(unstruct.alph.05.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = unstruct.alph.05.cp[[i]]$intervals[1 ,1], x1 = unstruct.alph.05.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(spillover.alph.05.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = spillover.alph.05.cp[[i]]$intervals[1 ,1], x1 = spillover.alph.05.cp[[i]]$intervals[2 ,1])
  }
  
  plot(y = c(1, 1), x = c(struct.alph.2.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = struct.alph.2.cp[[i]]$intervals[1 ,1], x1 = struct.alph.2.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(unstruct.alph.2.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = unstruct.alph.2.cp[[i]]$intervals[1 ,1], x1 = unstruct.alph.2.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(spillover.alph.2.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = spillover.alph.2.cp[[i]]$intervals[1 ,1], x1 = spillover.alph.2.cp[[i]]$intervals[2 ,1])
  }
  
  plot(y = c(1, 1), x = c(struct.alph1.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = struct.alph1.cp[[i]]$intervals[1 ,1], x1 = struct.alph1.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(unstruct.alph1.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = unstruct.alph1.cp[[i]]$intervals[1 ,1], x1 = unstruct.alph1.cp[[i]]$intervals[2 ,1])
  }
  plot(y = c(1, 1), x = c(spillover.alph1.cp[[1]]$intervals[ ,1]), xlab = "Changepoint interval", ylab = "", type = "l", xlim = c(0, 60), ylim = c(0, (reps + 1)))
  for(i in 2:reps){
    segments(y0 = i, y1 = i, x0 = spillover.alph1.cp[[i]]$intervals[1 ,1], x1 = spillover.alph1.cp[[i]]$intervals[2 ,1])
  }
  
#   hist(diff.struct.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.unstruct.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.spillover.alph.05, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.struct.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.unstruct.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.spillover.alph.2, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.struct.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.unstruct.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
#   hist(diff.spillover.alph1, col = "grey70", main = "", breaks = breaks.in, xlab = "difference in AIC", xlim = c(0, 200))
  return(list(struct.alph.05.cp =struct.alph.05.cp,
                struct.alph.2.cp = struct.alph.2.cp,
                struct.alph1.cp = struct.alph1.cp,
                struct.alph.05.lm = struct.alph.05.lm,
                struct.alph.2.lm = struct.alph.2.lm,
                struct.alph1.lm = struct.alph1.lm,
                diff.struct.alph.05 = diff.struct.alph.05,
                diff.struct.alph.2 = diff.struct.alph.2,
                diff.struct.alph1 = diff.struct.alph1,
                unstruct.alph.05.cp = unstruct.alph.05.cp,
                unstruct.alph.2.cp = unstruct.alph.2.cp,
                unstruct.alph1.cp = unstruct.alph1.cp,
                unstruct.alph.05.lm = unstruct.alph.05.lm,
                unstruct.alph.2.lm = unstruct.alph.2.lm,
                unstruct.alph1.lm = unstruct.alph1.lm,
                diff.unstruct.alph.05 = diff.unstruct.alph.05,
                diff.unstruct.alph.2 = diff.unstruct.alph.2,
                diff.unstruct.alph1 = diff.unstruct.alph1,
                spillover.alph.05.cp = spillover.alph.05.cp,
                spillover.alph.2.cp = spillover.alph.2.cp,
                spillover.alph1.cp =  spillover.alph1.cp,
                spillover.alph.05.lm = spillover.alph.05.lm,
                spillover.alph.2.lm = spillover.alph.2.lm,
                spillover.alph1.lm = spillover.alph1.lm ,
                diff.spillover.alph.05 = diff.spillover.alph.05,
                diff.spillover.alph.2 = diff.spillover.alph.2 
  ))
  
}