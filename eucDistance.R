eucDistance <- function(x, Y){
  
  n <- nrow(Y)
  
  out <- rep(0, n)
  
  for (i in 1:n){
    
    out[i] <- sqrt(sum( (x - Y[i,]) ^ 2 ))
    
  }
  
  
  return(out)
  
}
