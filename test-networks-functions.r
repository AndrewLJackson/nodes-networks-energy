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

# source relevant files
source('generateRandomTransitionMatrix.r')
source('odeEquations.r')


#-------------------------------------------------------------------------------
# Define parameters for the Transition Matrix

# number of nodes in a network to create
n.nodes <- 10

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

a <- generateRandomTransitionMatrix(n.nodes, p.connect, upper.tri = FALSE,
                                    a.upper, a.lower, e.upper, e.lower)

#-------------------------------------------------------------------------------
# Define other parameters of the system
#-------------------------------------------------------------------------------

# Energy input to the system is defined  by the vector f
f <- double(n.nodes)
# in this simple model, only node1 acquires energy from the environment
f[1] <- 1 

# determine noise to apply to node == 1


# solve for equilibrium
G.star <- -solve(a) %*% f



#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------


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
pars  <- c(a = a, f = f,
           node.hit = 1, prop.fall = 0.1, event.t = 100)

# initial conditions of the system at equilibrium
yini  <- c(G = G.star)

# specify the times at which we want to evaluate the system
times <- seq(1, 200, by = 1)

# events currently hardcoded to occur at t=100, so here i 
# get it to evaluate the system at t=100+error
times <- sort(c(times, 100 + .Machine$double.eps ))

#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------

# Simulate the system over time with a call to ode()
out   <- ode(yini, times, energyFlow, pars,
             events = list(func = node.fall, time = pars["event.t"]))



#-------------------------------------------------------------------------------
# Plot the results, and the network structure
#-------------------------------------------------------------------------------


# Plot the energy per node over time
#dev.new(height = 5, width = 5)
matplot(out[,1], out[,2:(n.nodes+1)], type="l", main = "ODE model", 
        xlab = "time", ylab = "energy in each node", 
        lwd = c(2, rep(1, n.nodes-1)), bty = "L",
        xlim=c(0, max(times) * 1.3))
# Add a grey vertical line to indicate the time at which the perturbation event
# is to be applied (if there is one specified)
abline(v = pars["event.t"], col="grey")
text((max(times)*1.05 + 2*(1:n.nodes)),
     out[nrow(out), 2:(n.nodes+1)], 1:n.nodes, cex = 0.75)

# add the analytically derived equilibrium points
points(rep(max(times), length(G.star)), G.star, pch = 19)







