UpdateLeslieFun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, 
                            joint.posterior.coda, posterior.names, intro.cost)
{
  # function to update Leslie matrix parameters with new draws from posterior of 
  # designated environmental state
  #
  # Args:
  # 
  # current.state = character vector with label of current state
  # sex.ratio = sex ratio of lambs (usually .1)
  # samples.to.draw = numerical range of chain steps that are outside the burn-in period in joint.posterior.coda
  # tot.chains = total number of Markov chains in joint.posterior.coda
  # joint.posterior.coda = posterior coda object output from specified IPM
  # posterior.names = list of names corresponding to each vector in the coda object
  # intro.cost = increase in overall mortality rate associated with disease introduction
  #
  # Returns:
  #
  # leslie = new leslie matrix
  
  chain <- ceiling(runif(1, 0, tot.chains))
  post.draw <- as.data.frame(t(joint.posterior.coda[[chain]][sample(x = samples.to.draw, size = 1), ]))
  names(post.draw) <- posterior.names
  
  post.draw$beta.wean.1.2 <- exp(post.draw$beta.wean.1.2) / (1 + exp(post.draw$beta.wean.1.2))
  post.draw$beta.wean.1.3 <- exp(post.draw$beta.wean.1.3) / (1 + exp(post.draw$beta.wean.1.3))
  post.draw$beta.wean.1.4 <- exp(post.draw$beta.wean.1.4) / (1 + exp(post.draw$beta.wean.1.4))
  post.draw$beta.wean.1.5 <- exp(post.draw$beta.wean.1.5) / (1 + exp(post.draw$beta.wean.1.5))
  post.draw$beta.wean.1.6 <- exp(post.draw$beta.wean.1.6) / (1 + exp(post.draw$beta.wean.1.6))
  post.draw$beta.adsurv.1.2 <- exp(post.draw$beta.adsurv.1.2) / (1 + exp(post.draw$beta.adsurv.1.2))
  post.draw$beta.adsurv.1.3 <- exp(post.draw$beta.adsurv.1.3) / (1 + exp(post.draw$beta.adsurv.1.3))
  post.draw$beta.adsurv.1.4 <- exp(post.draw$beta.adsurv.1.4) / (1 + exp(post.draw$beta.adsurv.1.4))
  post.draw$beta.adsurv.1.5 <- exp(post.draw$beta.adsurv.1.5) / (1 + exp(post.draw$beta.adsurv.1.5))
  post.draw$beta.adsurv.1.6 <- exp(post.draw$beta.adsurv.1.6) / (1 + exp(post.draw$beta.adsurv.1.6))
  
  post.draw$beta.wean.2.2 <- exp(post.draw$beta.wean.2.2) / (1 + exp(post.draw$beta.wean.2.2))
  post.draw$beta.wean.2.3 <- exp(post.draw$beta.wean.2.3) / (1 + exp(post.draw$beta.wean.2.3))
  post.draw$beta.wean.2.4 <- exp(post.draw$beta.wean.2.4) / (1 + exp(post.draw$beta.wean.2.4))
  post.draw$beta.wean.2.5 <- exp(post.draw$beta.wean.2.5) / (1 + exp(post.draw$beta.wean.2.5))
  post.draw$beta.adsurv.2.2 <- exp(post.draw$beta.adsurv.2.2) / (1 + exp(post.draw$beta.adsurv.2.2))
  post.draw$beta.adsurv.2.3 <- exp(post.draw$beta.adsurv.2.3) / (1 + exp(post.draw$beta.adsurv.2.3))
  post.draw$beta.adsurv.2.4 <- exp(post.draw$beta.adsurv.2.4) / (1 + exp(post.draw$beta.adsurv.2.4))
  post.draw$beta.adsurv.2.5 <- exp(post.draw$beta.adsurv.2.5) / (1 + exp(post.draw$beta.adsurv.2.5))
  
  if(current.state == "healthy"){
    repros <- c(0, rep((post.draw$beta.wean.1.2 * sex.ratio), 1), rep((post.draw$beta.wean.1.3  * sex.ratio), 5), rep((post.draw$beta.wean.1.4  * sex.ratio), 6), rep((post.draw$beta.wean.1.5  * sex.ratio), 6))    
    survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5))    
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, dim(diag(survs))[1])))
  }
  else if(current.state == "infected"){
    repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
    survs <- c(1, rep(post.draw$beta.adsurv.2.2, 1), rep(post.draw$beta.adsurv.2.3, 5), rep(post.draw$beta.adsurv.2.4, 6), rep(post.draw$beta.adsurv.2.5, 5))    
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  } 
  else if(current.state == "spillover"){
    repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
    survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5)) * (1 - intro.cost)
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs)))) 
  }
  return(leslie)
}
