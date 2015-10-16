energyStabilityMetrics <- function(euc.dist){
  
  # Resistance
  max.idx <- which.max(euc.dist)
  max.deflect <- euc.dist[max.idx]
  max.deflect.time <- times[max.idx]
  
  resistance <- c(deflection = max.deflect, time = max.deflect.time)
  
  # Resilience
  
  idx <- which(euc.dist <= max.deflect / exp(1) & 
                 times > max.deflect.time)
  
  decay <- times[idx[1]]
  
  idx <- which(euc.dist <= quantile(euc.dist, 0.05, na.rm=T) & 
                 times > max.deflect.time)
  
  res.quant <- times[idx[1]]
  
  idx <- which(euc.dist <= 10^-4 & 
                 times > max.deflect.time)
  
  res.tol <- times[idx[1]]
  
  resilience <- c(decay = decay, quant = res.quant, tol = res.tol)
  
  out <- list()
  out$resistance <- resistance
  out$resilience <- resilience
  
  return(out)

}
