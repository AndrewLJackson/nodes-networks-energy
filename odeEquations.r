energyFlow <- function(t, G, Pars) {
  with(as.list(Pars), {

    dG <- (a %*% G) + f
    
    return(list(dG))
  })
}

# Include a perturbation event function which can be evaluated at 
# defined time points in the call to ode()
node.fall <- function(t, G, pars){
  with (as.list(pars),{
    G[node.hit] <- G[node.hit] * prop.fall
    return(G)
  })
  
}
