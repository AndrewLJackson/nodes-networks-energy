# About: create random networks of nodes, among which energy flows 
# according to system of simple linear ODEs.
# Author: Andrew Jackson
# Date created: 22/11/14
# Based on: http://books.google.ie/books?id=hZVMPYaK6GgC&pg=PA104&lpg=PA104&dq=compartment+model+energy+flow&source=bl&ots=_yZyDalndi&sig=SVdQEUG6B031AEyoiYmb9XIZqbE&hl=en&sa=X&ei=hpxwVMwyhrE88O2BuAs&ved=0CC8Q6AEwAQ#v=onepage&q=compartment%20model%20energy%20flow&f=false



#-------------------------------------------------------------------------------
# Initiation
rm(list=ls())  # clear memory
graphics.off() # close open graphics windows from previous runs
set.seed(1)   # for debugging

# This script requires installation of these packages, along with their
# dependencies
# in the R command window, type
# install.packages(c('deSolve', 'diagram'))
# to pull them down from the cran repository.
library('deSolve')
library('diagram')

#-------------------------------------------------------------------------------
# Network set up
# I am toying with turning this portion of code into a function

# number of nodes in a network to create
n.nodes <- 10

# specify the binary interaction matrix B
# this matrix deterimes which nodes are connected.
# This is defined here symmetically, so that an amount of energy g
# is removed from one node, and added to another. But see later, we 
# can specify energy losses which essentially makes these asymmetric,
# with an amount g lost from one node, and a corresponding amount qg
# gained by the other node with 0 <= g <= 1


# probability of a connection
p.connect <- 0.6

B <- matrix(rbinom(n.nodes^2, 1, p.connect), ncol = n.nodes, nrow = n.nodes)


# only fill in the lower triangle of B
# nodes in the lower triangle consume energy from the 
# node specified by each column. That is, if a[2,1] <-1
# then node2 gains one unit from node1 and so node1 has to
# have a negative value associated in its diagonal.
B[upper.tri(B)] <- 0
B[diag(B)] <- 0

# specify the weighted interaction matrix A
# this matrix determines the rates of energy flow between nodes.
# We will let energy production stay constant at 1, and scale the others
# below this amount.

# upper and lower bounds for values in A
a.ub <- 0.2
a.lb <- 0.1

aa <- matrix(runif(n.nodes^2, a.lb, a.ub), 
	          ncol = n.nodes, nrow = n.nodes)

# Create the transition matrix a for each valid connection specified by 
# the binary connection matrix B
a <- aa * B


# Each node loses energy. This simulates energy loss during transfer,
# and also potentially acts as a feedback of energy to the environment,
# but that would need some re-coding to achieve.
e.ub <- 0.15 # energy loss lower bound
e.lb <- 0.05 # energy loss upper bound
e <- runif(n.nodes, e.lb, e.ub)


# calculate the diagonal values representing sums of losses to connected nodes
# in addition to loss to the environment
diag(a) <- -colSums(a) - e

# Energy input to the system is defined  by the vector f
f <- double(n.nodes)
# in this simple model, only node1 acquires energy from the environment
f[1] <- 1 
 
# set initial conditions, i.e. Energy per node.
# HEre i am going to just seed the network with 1 unit of energy
# in the energy node#1
G.0 <- matrix(0, nrow = n.nodes, ncol = 1)
G.0[1] <- 1


# solve for equilibrium
G.star <- -solve(a) %*% f

#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------


energy.flow <- function(t, G, Pars) {
  with(as.list(Pars), {

    if (t >= event.t) {
      a[,node.loss.id] <- 0
      a[node.loss.id] <- 0
      }
    dG <- (a %*% G) + f #+ c(rnorm(0,0.1), rep(0, nrow(G-1)))

    return(list(dG))
  })
}

# define the parameters for passing to the energy.flow() function
# The following parameters control the perturbation:
#     node.hit = integer identifying which node will be proportionally changed
#     prop.fall = proportional change to node.hit. If == 1, then no change
#     node.loss.id = integer identifying which node to remove entirely
#     event.t = time point to lose node.
# Notes: set even.t = Inf to not apply any perturbation.
# If removing a node, then you should also use node.hit = node.loss.id and 
# prop.fall = 0 to set its state variable value to zero, in addition to removing
# its connections. You can turn off a proportional change in any node by 
# setting prop.fail = 1.
pars  <- c(a = a, f = f,
           node.hit = 1, prop.fall = 0.1,
           node.loss.id = NA, event.t = 100)

# initial conditions of the system
yini  <- c(G = G.0)

# specify the times at which we want to evaluate the system
times <- seq(1, 200, by = 1)

# events currently hardcoded to occur at t=100, so here i 
# get it to evaluate the system at t=100+error
times <- sort(c(times, 100 + .Machine$double.eps ))


# Include a perturbation event function which can be evaluated at 
# defined time points in the call to ode()
node.fall <- function(t, G, pars){
  with (as.list(pars),{
    G[node.hit] <- G[node.hit] * prop.fall
    return(G)
  })

}


#*******************************************************************************
# This is where the system is actually evaluated.
# Simulate the system over time with a call to ode()
out   <- ode(yini, times, energy.flow, pars,
             events = list(func = node.fall, time = pars["event.t"]))
#*******************************************************************************

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












