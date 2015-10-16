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


#-------------------------------------------------------------------------------
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

a <- generateRandomTransitionMatrix(n.nodes, p.connect, upper.tri = FALSE,
                                    a.upper, a.lower, e.upper, e.lower)

#-------------------------------------------------------------------------------
# Define other parameters of the system
#-------------------------------------------------------------------------------

# Energy input to the system is defined  by the vector f
f <- double(n.nodes)
# in this simple model, only node[1] acquires energy from the environment
f[1] <- 1 


# check that the system is stable
if (!all(eigen(a)$values <= 0)){warning("One of the eigen values is positive")}

# solve for equilibrium
G.star <- -solve(a) %*% f



#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------




# initial conditions of the system at equilibrium
yini  <- c(G = G.star)

# specify the times at which we want to evaluate the system
times <- seq(1, 200, by = 0.25)

# events currently hardcoded to occur at t=100, so here i 
# get it to evaluate the system at t=100+error
times <- sort(c(times, 100 + .Machine$double.eps ))


# AJ - Adding noise
# Either: 1) use something like tuneR::noise() to generate a time series, then
# approxfun() within the ode defined function to estimated noise at time=t; or
# 2) check out the book that goes with package sde and use the BM or OU 
# processes from that. I am leaning towards the latter.

d <- expression(0 - 20 * x)

s <- expression(0.1) 

ee <- genNoiseFun(d, s, N=10^3, min(times), max(times))



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
           node.hit = 1, prop.fall = 0.1, event.t = 100,
           e = double(n.nodes),
           noiseFun = ee)



#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------

# Simulate the system over time with a call to ode()
perturbed <- ode(yini, times, energyFlow, pars,
             events = list(func = node.fall, time = pars["event.t"]))

control   <- ode(yini, times, energyFlow, pars)


#-------------------------------------------------------------------------------
# Calculate stability metrics
#-------------------------------------------------------------------------------

# euclidean distance

difference <- perturbed[,2:(n.nodes+1)] - control[,2:(n.nodes+1)]

euc.dist <- sqrt(rowSums( difference^ 2))

plot(times, euc.dist, type = "l")

max.idx <- which.max(euc.dist)
max.deflect <- euc.dist[max.idx]
max.deflect.time <- times[max.idx]

tol <- 

#-------------------------------------------------------------------------------
# Plot the results, and the network structure
#-------------------------------------------------------------------------------


# Plot the energy per node over time
#dev.new(height = 5, width = 5)
matplot(perturbed[,1], perturbed[,2:(n.nodes+1)], type="l",
        main = "Perturbed System Dynamics", 
        xlab = "time", ylab = "energy in each node", 
        lwd = c(2, rep(1, n.nodes-1)), bty = "L",
        xlim=c(0, max(times) * 1.3))
# Add a grey vertical line to indicate the time at which the perturbation event
# is to be applied (if there is one specified)
abline(v = pars["event.t"], col="grey")
text((max(times)*1.05 + 2*(1:n.nodes)),
     perturbed[nrow(perturbed), 2:(n.nodes+1)], 1:n.nodes, cex = 0.75)

# add the analytically derived equilibrium points
points(rep(min(times), length(G.star)), G.star, pch = 19)
points(rep(max(times), length(G.star)), G.star, pch = 19)

# ------------------------------------------------------------------------------
# Plot perturbed - control

matplot(perturbed[,1], perturbed[,2:(n.nodes+1)] - control[,2:(n.nodes+1)],
        type="l", main = "Perturbed - Control", 
        xlab = "time", ylab = "energy in each node", 
        lwd = c(2, rep(1, n.nodes-1)), bty = "L",
        xlim=c(90, 140))
# Add a grey vertical line to indicate the time at which the perturbation event
# is to be applied (if there is one specified)
abline(v = pars["event.t"], col="grey")
#text((max(times)*1.05 + 2*(1:n.nodes)),
#     perturbed[nrow(out), 2:(n.nodes+1)], 1:n.nodes, cex = 0.75)

#-------------------------------------------------------------------------------
# Use pkg 'diagram' to visualise the network
# I intend to tidy up this to make it more visually appealing, and
# accurate regarding the energy gains and losses to and from the environment.

a.cnx <- a
diag(a.cnx) <- 0
a.self <- diag(a)


#dev.new()
pp <- plotmat(round(a.cnx, digits = 2),
              curve = 0.1,
              lwd = 1, box.lwd = 2, cex.txt = F,
              box.type = "square", box.prop = 0.8,  box.size = 0.05,
              arr.type = "triangle", arr.length = 0.3,
              arr.width = 0.3 / 2,
              arr.pos = 0.7,
              arr.lwd = a.cnx, 
              shadow.size = 0, prefix = "",
              relsize = 0.75,
              self.lwd = 1,
              self.shiftx = 0.07, self.shifty = 0.07, 
              main = "Energy network", 
              box.col = c('red', rep('white',n.nodes-1)))


# add red arrows representing gain from the environment
dx <- 0.06
dy <- 0.06

for ( i in 1:n.nodes){
  if (f[i]) {
    straightarrow(from = c(pp$rect[i,"xright"] + dx, 
                           mean(pp$rect[i,c("ytop","ybot")])),
                  to =   c(pp$rect[i,"xright"],
                           mean(pp$rect[i,c("ytop","ybot")])),
                  lcol = "red", arr.col = "red", lty = 1,
                  arr.pos = 0, arr.type = "triangle"
    )
  }
}


# add blue arrows representing loss to the environment
for ( i in 1:n.nodes){
  if (a.self[i]) {
    straightarrow(from = c(pp$rect[i,"xleft"], 
                           pp$rect[i,"ytop"]),
                  to =   c(pp$rect[i,"xleft"] - dx,
                           pp$rect[i,"ytop"]  + dy),
                  lcol = "blue", arr.col = "blue", lty = 3,
                  arr.pos = 1, arr.type = "triangle"
    )
  }
}






