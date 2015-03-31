InchaustiSpectralColor <- function(data){
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
  fft.out <- fft(data$Ewes[1:(dim(data)[1] - 1)])
  
  # 2. note that spectrum values will generally be complex numbers
  # 3. Square abs val of each value in spectrum
  sqrd.abs <- abs(fft.out) ^ 2
  # 4. spectral redness exponent is the slope of the least squares regression
  #    line on this plot on the log-log scale
  plot(sqrd.abs ~ seq(1:(dim(data)[1] - 1))^2, log = "xy")
  loglog.fit <- lm(log(sqrd.abs) ~ seq(1:(dim(data)[1] - 1))^2)
  # 5. Reddened spectrum is associated with a line of negative slope
  redness.exp <- summary(loglog.fit)$coefficients[2, 1]
  
  return(list(fft.out = fft.out,
              loglog.fit = loglog.fit, 
              redness.exp = redness.exp))
}