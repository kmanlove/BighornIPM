pupilplot <- function (wf, cp = NULL, col = topo.colors(256), addContours = FALSE, 
                       cscale = TRUE, ...) 
{
  if (cscale) {
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
    par(las = 1)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar) 
    thelist <- list(...)  
    findz <- which(names(thelist) == 'zlim')  
    if (length(findz) > 0 ) {   
      zlim <- thelist$zlim  
    }else{  
      zlim <- range(wf, finite = TRUE) #the original line  
    } 
    # end of my hack  
    levels <- seq(zlim[1], zlim[2], length = length(col))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col,  density = NA)
    axis(4)
    box()
    mar <- mar.orig
    mar[4] <- 0
    par(mar = mar)
  }
  if (is.null(cp)) {
    axis1 <- 1:nrow(wf)
    axis2 <- 1:ncol(wf)
  }
  else {
    axis1 <- ((1:nrow(wf)) - cp$xc)/cp$rx
    axis2 <- ((1:ncol(wf)) - cp$yc)/cp$ry
  }
  image(axis1, axis2, wf, col = col, asp = 1, xlab = "Years to Reintroduction", ylab = "Years to Fade-out",  ...)
  if (addContours) 
    contour(axis1, axis2, wf, add = TRUE)
}