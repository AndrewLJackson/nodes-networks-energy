#' Evaluates a set of energy-node-network ODEs using same noise function.
#' 
#' A series OU noise functions are specified using the \code{d} and \copde{s}
#' expresssion parameters. The random number is set to the same seed each time
#' the noise series is generated, thereby creating a set of paired experiments
#' where the effect of noise on the stability of the system can be evaluated in
#' isolation. The system is evaluated for all pairs of \code{d} and \copde{s}
#' provided, with a perturbation applied to node[1] at the specified time point,
#' and a control series is also evaluated with no such perturbation.
#' 
#' @param d.list A list containing a set of expressions specifying the drift 
#'   fuction for passing to \code{\link{sde:sde.sim}}
#' @param s.list A list containing a set of expressions specifying the sigma 
#'   fuction for passing to \code{\link{sde:sde.sim}}
#' @param pars the set of parameters for passing to \code{\link{ode}}
#' @param times the set of time points at which the ode is to be evaluated at:
#'   for passing to \code{\link{ode}}
#'   
#' @return A list of matrices for each evaluation of the system, along with
#'   their corresponding noise series, and random seed.
#'   


evalPairedSims <- function(d.list, s.list, pars, times){
  
  # get the current system time for the seed to use for all subsequent noise
  # series in this repeated experiment.
  seed <- Sys.time
  
  # loop over d.list
  for (i in 1:length(d.list)){
    
    # loop over s.list
    for (j in 1:length(s.list)){
      
      # re-set the seed
      set.seed(seed)
      
      # generate the noise funcition
      ee <- genNoiseFun(d.list[[i]], s.list[[j]], N = 10^3, 
                        min(times), max(times))
      
      
      
      
    }
    
  }
  
  set.seed(seed)
  
  # Generate the noise function
  
  
  
}
