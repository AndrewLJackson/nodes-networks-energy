energyStabilityMetricsNoControl <- function(euc.dist, times, t.disturb){
  
  
  # pre-perturbation behaviour
  q.95 <- quantile(euc.dist[times < t.disturb], 0.95)
  
  # ----------------------------------------------------------------------------
  # Resistance
  
  # find the max deflection that occurs after the perturbation
  max.idx <- which.max(euc.dist & times > t.disturb)
  max.deflect <- euc.dist[max.idx]
  max.deflect.time <- times[max.idx]
  
  resistance <- c(deflection = max.deflect, time = max.deflect.time)
  
  # ----------------------------------------------------------------------------
  # Resilience
  
  idx <- which(euc.dist <= max.deflect / exp(1) & 
                 times > max.deflect.time)
  
  decay <- times[idx[1]]
  
  idx <- which(euc.dist <= q.95 & 
                 times > max.deflect.time)
  
  res.quant <- times[idx[1]]
  
  resilience <- c(decay = decay, quant = res.quant)
  
  out <- list()
  out$resistance <- resistance
  out$resilience <- resilience
  out$q.95 <- q.95
  
  return(out)

}

# --------------------------------------------------------------



## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y))
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste0(prefix, txt)
#   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor * r)
# }
# 
# pairs(perturbed[,2:(n.nodes+1)], lower.panel = panel.smooth, upper.panel = panel.cor)

