PlotPostVitalRates <- function(beta.posts.ad.surv, beta.posts.wean, write, write.path)
{
  # Function to plot posteriors on all beta-estimates from IPM
  #
  # Args:
  # beta.posts.ad.surv = 18 x 5 matrix of posterior quantiles for all ad survival betas
  # beta.posts.wean = 18 x 5 matrix of posterior quantiles for all lamb weaning betas
  # write = logical for whether plot should be written out (T == write plot to file)
  # write.path = path to directory where file will be written if write == T
  #
  # Returns:
  # 
  # Posterior density plot for all IPM-estimated vital rates

  # plot betas out
  plot.cols <- c("white", "black", "grey60")
  if(write == T){
#  svg("./Plots/FigsV1/MoviDefFigs_05Jan2015/Figure2_Posteriors_MoviDef_05Jan2015.svg", width = 3, height = 2, pointsize = 8)
    svg(paste(write.path), width = 3, height = 2, pointsize = 8)
  } 

  par(mfrow = c(1, 2), las = 1, oma = c(4, 6, 2, 2), mar = c(3, 5, 0, 0))
  plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of survival", xlab = "", xaxt = "n")
  for(i in 3:15){
    segments(x0 = i, x1 = i, y0 = (exp(beta.posts.adsurv[i, 1])) / (1 + exp(beta.posts.adsurv[i, 1])), y1 = exp(beta.posts.adsurv[i, 5]) / (1 + exp(beta.posts.adsurv[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3), lwd = 2)
    segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.adsurv[i, 3])) / (1 + exp(beta.posts.adsurv[i, 3])), y1 = exp(beta.posts.adsurv[i, 3]) / (1 + exp(beta.posts.adsurv[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3), lwd = 2)
  }
  axis(side = 1, at = c(4.5, 7.5, 10.5, 13.5), labels = c("2", "3-7", "8-13", "> 13"))
  mtext(side = 1, line = 2, outer = F, "Ewe age")
  
  plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of producing a recruit", xlab = "", xaxt = "n")
  for(i in 3:15){
    segments(x0 = i, x1 = i, y0 = (exp(beta.posts.wean[i, 1])) / (1 + exp(beta.posts.wean[i, 1])), y1 = (exp(beta.posts.wean[i, 5])) / (1 + exp(beta.posts.wean[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3), lwd = 2)
    segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.wean[i, 3])) / (1 + exp(beta.posts.wean[i, 3])), y1 = exp(beta.posts.wean[i, 3]) / (1 + exp(beta.posts.wean[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3), lwd = 2)
  }
  axis(side = 1, at = c(4.5, 7.5, 10.5, 13.5), labels = c("2", "3-7", "8-13", "> 13"))
  mtext(side = 1, line = 2, outer = F, "Ewe age")
  #dev.off()
  
}