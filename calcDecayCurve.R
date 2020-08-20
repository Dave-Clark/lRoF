calcDecayCurve <- function(x){
  xCoord <- seq(min(x$data[, 1]), max(x$data[, 1]),  0.01)
  pred <- 1 - (1 - x$a.intercept) * exp(-x$b.slope * xCoord)
  decayFit <- data.table(x = xCoord, y = pred)
  return(decayFit)
}
