# script to loop over replicates of calls to system perturbation

#-------------------------------------------------------------------------------
# Initiation
rm(list=ls())  # clear memory
graphics.off() # close open graphics windows from previous runs
#set.seed(1)   # for debugging

# This script requires installation of these packages, along with their
# dependencies
# in the R command window, type
# install.packages(c('deSolve', 'diagram'))
# to pull them down from the cran repository.
library('deSolve')
library('diagram')
library('sde')


# source relevant files
source('generateRandomTransitionMatrix.r')
source('odeEquations.r')
source('genNoiseFun.r')
source('eucDistance.R')
source('energyStabilityMetricsNoControl.r')
source('evalNoControlSims.R')


#-------------------------------------------------------------------------------
# global parameters

# replicates
reps <- 30

# set of sigma values for noise generation
s.list <- list(expression(0.1),
               expression(0.5),
               expression(2),
               expression(5)
               )

d.list <- list( expression(0 - 100 * x) )

# specify the times at which we want to evaluate the system
times <- seq(1, 500, by = 0.25)

event.t <- 200

# events currently hardcoded to occur at t=100, so here i 
# get it to evaluate the system at t=100+error
times <- sort(c(times, event.t + .Machine$double.eps ))

# Define parameters for the Transition Matrix

# number of nodes in a network to create
n.nodes <- 8

# probability of a connection
p.connect <- 0.6

# upper and lower bounds for values in A
a.upper <- 0.2
a.lower <- 0.1

# Each node loses energy. This simulates energy loss during transfer,
# and also potentially acts as a feedback of energy to the environment,
# but that would need some re-coding to achieve.
e.upper <- 0.15 # energy loss lower bound
e.lower <- 0.05 # energy loss upper bound

# Energy input to the system is defined  by the vector f
f <- double(n.nodes)
# in this simple model, only node[1] acquires energy from the environment
f[1] <- 1 


# ------------------------------------------------------------------------------
# loop over replicates
# ------------------------------------------------------------------------------

# prep results matrices.... lazily here
experiments <- length(d.list)*length(s.list) 

resistance <- matrix(NA, reps*experiments, 2)
resilience <- matrix(NA, reps*experiments, 2)
sigma <- matrix(NA, reps*experiments, 1)

# a counter
ct <- 1
ct.rep <- 0

for (i in 1:reps){
  
  a <- generateRandomTransitionMatrix(n.nodes, p.connect, upper.tri = FALSE,
                                      a.upper, a.lower, e.upper, e.lower)
  
  # check that the system is stable
  if (!all(eigen(a)$values <= 0)){warning("One of the eigen values is positive")}
  
  # solve for equilibrium
  G.star <- -solve(a) %*% f
  
  # define the parameters for passing to the energy.flow() function
  # The following parameters control the perturbation:
  #     node.hit = integer identifying which node will be proportionally changed
  #     prop.fall = proportional change to node.hit. If == 1, then no change
  #     event.t = time point to lose node.
  # Notes: set even.t = Inf to not apply any perturbation.
  # If removing a node, then you should also use node.hit = node.loss.id and 
  # prop.fall = 0 to set its state variable value to zero, in addition to removing
  # its connections. You can turn off a proportional change in any node by 
  # setting prop.fail = 1.
  pars  <- list(a = a, f = f,
                node.hit = 1, prop.fall = 0.1, event.t = event.t,
                e = double(n.nodes))
  
  out <- evalNoControlSims(d.list, s.list, pars, times, G.star)
  
  # store results
  resistance[ct:(ct+experiments-1), ] <- out$resistance
  resilience[ct:(ct+experiments-1), ] <- out$resilience
  sigma[ct:(ct+experiments-1), ] <- unlist(s.list)
  
  # update counter
  ct <- ct + experiments
  ct.rep <- ct.rep + 1
  print(ct.rep)
}

# ------------------------------------------------------------------------------
# plot results
# ------------------------------------------------------------------------------

par(mfrow = c(2,2))
boxplot(resistance[,1] ~ sigma, main = "Max Deflection", ylab = "euclidean distance")
boxplot(resilience[,1] - event.t ~ sigma, main = "resilience (decay)", ylab = "time")
boxplot(resilience[,2] - event.t ~ sigma, main = "resilience (quantile)", ylab = "time")
boxplot(log(resilience[,1] - event.t) ~ sigma, main = "log resilience (decay)", ylab = "time")

# ------------------------------------------------------------------------------
save(resistance, resilience, d.list, s.list, 
     file = "perturbation.rdata")



