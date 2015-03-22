TuljarFun3 <- function(Tuljar2.out){
  D <- Tuljar2.out$b / Tuljar2.out$lambda.naught - Tuljar2.out$u.naught %*% t(Tuljar2.out$v.naught)
  
  Cx.he <- Tuljar2.out$L.healthy - Tuljar2.out$b
  Cx.sp <- Tuljar2.out$L.sp - Tuljar2.out$b
  Cx.inf <- Tuljar2.out$L.infected - Tuljar2.out$b
  
  return(list(D = D, Cx.he = Cx.he, Cx.sp = Cx.sp, Cx.inf = Cx.inf))
}