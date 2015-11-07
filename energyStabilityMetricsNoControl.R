energyStabilityMetricsNoControl <- function(euc.dist, times, t.disturb){
  
  
  # pre-perturbation behaviour
  q.95 <- quantile(euc.dist[times < t.disturb], 0.95)
  
  # fit lowess to post-perturbation behaviour
  smo <- lowess(times >= t.disturb, euc.dist[times >= t.disturb])
  
  
  # Resistance
  max.idx <- which.max(euc.dist)
  max.deflect <- euc.dist[max.idx]
  max.deflect.time <- times[max.idx]
  
  resistance <- c(deflection = max.deflect, time = max.deflect.time)
  
  # Resilience
  
  idx <- which(smo$y <= max.deflect / exp(1))
  
  decay <- smo$x[idx[1]]
  
  idx <- which(smo$y <= q.95 )
  
  res.quant <- smo$x[idx[1]]
  
  resilience <- c(decay = decay, quant = res.quant)
  
  out <- list()
  out$resistance <- resistance
  out$resilience <- resilience
  
  return(out)

}
