energyFlow <- function(t, G, Pars) {
  with(Pars, {

    # e is the noise term applied only to node==1
    e[1] <- noiseFun(t)
      
    dG <- (a %*% G) + f + e #+ noiseFun(t)
    
    return(list(dG))
  })
}

# Include a perturbation event function which can be evaluated at 
# defined time points in the call to ode()
node.fall <- function(t, G, pars){
  with (as.list(pars),{
    G[node.hit] <- G[node.hit] * prop.fall
    return(as.vector(G))
  })
  
}
