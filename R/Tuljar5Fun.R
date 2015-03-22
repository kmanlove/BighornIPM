TuljarFun5 <- function(Tuljar1.out, Tuljar2.out, Tuljar3.out, Tuljar4.out)
{
  g1 <- (1 / Tuljar2.out$lambda.naught) * cbind(Tuljar3.out$Cx.he, Tuljar3.out$Cx.sp, Tuljar3.out$Cx.inf)
  g2 <- (1 / Tuljar2.out$lambda.naught) * rbind(Tuljar1.out$w.naught[1] * t(Tuljar3.out$Cx.he), 
                                                Tuljar1.out$w.naught[2] * t(Tuljar3.out$Cx.sp), 
                                                Tuljar1.out$w.naught[3] * t(Tuljar3.out$Cx.inf))
  
  return(list(g1 = g1, g2 = g2))
}
