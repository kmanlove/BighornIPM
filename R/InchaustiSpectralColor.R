InchaustiSpectralColor <- function(data, PopEst){
  # VARIANCE GROWTH EXPONENT
  # Adapted from Inchausti and Halley (Science 2001) supp info
  # 1. log-transform timeseries (N + 1)
  # 2. for each window of size k (3 < k < n), calculate
    # variance (V_k)for all subsequences of size k in the timeseries
    # ALLOW FOR PARTIALLY OVERLAPPING WINDOWS.
  # 3. exponent is regression coef of log(V_k) ~ log(k)
  # for a stationary process, gamma (the slope) will converge to 0 as n gets large
  
  
  # SPECTRAL REDNESS EXPONENT
  # 1. apply discrete fourier transform to (full) timeseries
  #   (use appropriate windowing function)
  require(stats)
  if(PopEst == F){
  log.Ewes <- log(data$Ewes[1:(dim(data)[1] - 1)] - data$NoFemRel[1:(dim(data)[1] - 1)] + data$NoFemRem[1:(dim(data)[1] - 1)] + 1)
  ewes.diff <- diff(log.Ewes)[-1]
  m <- floor((dim(data)[1] - 1) / 2)
  ewes.spectrum <- spectrum(ewes.diff, method = "pgram")
  ewes.freq <- ewes.spectrum$freq
  ewes.density <- ewes.spectrum$spec
#  ewes.dft <- fft(ewes.diff)
#  ewes.dft2 <- ewes.dft[2 : (m + 1)] 
  plot(log(ewes.density) ~ log(ewes.freq))
  # 2. note that spectrum values will generally be complex numbers
#  # 3. Square abs val of each value in spectrum
#  sqrd.abs <- abs(fft.out) ^ 2
  # 4. spectral redness exponent is the slope of the least squares regression
  #    line on this plot on the log-log scale
#  years.inv <- ((1 / seq(1:(dim(data)[1] - 1))))
#  plot(sqrd.abs ~ years.inv, log = "xy")
  loglog.fit <- lm(log(ewes.density) ~ log(ewes.freq))
  # 5. Reddened spectrum is associated with a line of negative slope
  redness.exp <- -1 * summary(loglog.fit)$coefficients[2, 1]
  
  return(list(ewes.spectrum = ewes.spectrum,
              loglog.fit = loglog.fit, 
              redness.exp = redness.exp))
  } 

  else if(PopEst == T)
    {
    log.PopEst <- log(data$PopEst[1:dim(data)[1] - 1] + 1)
    m <- floor((dim(data)[1] - 1) / 2)
    popest.diff <- diff(log.PopEst)
    popest.spectrum <- spectrum(popest.diff, method = "pgram")
    popest.freq <- popest.spectrum$freq
    popest.density <- popest.spectrum$spec
    plot(log(popest.density) ~ log(popest.freq))
    loglog.fit <- lm(log(popest.density) ~ log(popest.freq))
    redness.exp <- -1 * summary(loglog.fit)$coefficients[2, 1]
    
    return(list(popest.spectrum = popest.spectrum,
                loglog.fit = loglog.fit, 
                redness.exp = redness.exp))
  }

}