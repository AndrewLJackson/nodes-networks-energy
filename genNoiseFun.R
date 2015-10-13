genNoiseFun <- function(d, s, N, tmin, tmax){
 
  # Ornstein-Uhlenbeck process
  
  y <- sde.sim(X0=0, drift=d, sigma=s, T=1, N=N, M=1)
  
  x <- seq(tmin, tmax, length = N+1)
  
  
  # linear interpolation of the process for passing to ode
  yy <- approxfun(x, y) 
  
  return(yy)
  
  
}





