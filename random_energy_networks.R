# About: create random networks of nodes, among which energy flows 
# according to simple linear ODEs.
# Author: Andrew Jackson
# Date created: 22/11/14
# Based on: http://books.google.ie/books?id=hZVMPYaK6GgC&pg=PA104&lpg=PA104&dq=compartment+model+energy+flow&source=bl&ots=_yZyDalndi&sig=SVdQEUG6B031AEyoiYmb9XIZqbE&hl=en&sa=X&ei=hpxwVMwyhrE88O2BuAs&ved=0CC8Q6AEwAQ#v=onepage&q=compartment%20model%20energy%20flow&f=false



#-------------------------------------------------------------------------------
# Initiation
rm(list=ls())
graphics.off()
set.seed(1)

library('deSolve')
library('diagram')

#-------------------------------------------------------------------------------
# Network set up

n.nodes <- 5

# specify the binary interaction matrix B
# this matrix deterimes which nodes are connected.
# This is defined here symmetically, so that an amount of energy g
# is removed from one node, and added to another. But see later, we 
# can specify energy losses which essentially makes these asymmetric,
# with an amount g lost from one node, and a corresponding amount qg
# gained by the other node with 0 <= g <= 1


# probability of a connection
p.connect <- 0.7

B <- matrix(rbinom(n.nodes^2, 1, p.connect), ncol = n.nodes, nrow = n.nodes)


# only fill in the lower triangle of B
# nodes in the lower triangle consume energy from the 
# node specified by each column. That is, if a[2,1] <-1
# then node2 gains one unit from node1 and so node1 has to
# have a negative value associated in its diagonal.
B[upper.tri(B)] <- 0
B[diag(B)] <- 0

#B[2,1] <- 1 # debug
#B[2,2] <- 0 # debug


# specify the weighted interaction matrix A
# this matrix determines the rates of energy flow between nodes.
# We will let energy production stay constant at 1, and scale the others
# below this amount.

# upper and lower bounds for values in A
a.ub <- 0.2
a.lb <- 0.1

aa <- matrix(runif(n.nodes^2, a.lb, a.ub), 
	          ncol = n.nodes, nrow = n.nodes)

a <- aa * B



# specify some energy losses to the environment
e.ub <- 0.15
e.lb <- 0.05
e <- runif(n.nodes, e.lb, e.ub)


# calculate the diagonal values representing sums of losses to nodes
# in addition to loss to the environment
diag(a) <- -colSums(a) - e


# Energy input to the system is defined  by the vector f
f <- double(n.nodes)
f[1] <- 1
 
# set initial conditions, i.e. Energy per node.
# HEre i am going to just seed the network with 1 unit of energy
# in the energy node#1
G.0 <- matrix(0, nrow = n.nodes, ncol = 1)
G.0[1] <- 1



#-------------------------------------------------------------------------------
# Evaluate the network as a set of coupled ODEs
#-------------------------------------------------------------------------------

## =======================================================================
## Example1: Predator-Prey Lotka-Volterra model (with logistic prey)
## =======================================================================

energy.flow <- function(Time, G, Pars) {
  with(as.list(c(Pars)), {

  
    dG <- (a %*% G) + f


    return(list(dG))
  })
}

pars  <- c(a = a, f = f)
yini  <- c(G = G.0)
times <- seq(0, 100, by = 1)
out   <- ode(yini, times, energy.flow, pars)
summary(out)


dev.new(height = 5, width = 5)
matplot(out[,1], out[,2:(n.nodes+1)], type="l", main = "ODE model", 
	     xlab = "time", ylab = "energy in each node")

dev.new()
pp <- plotmat(a, curve = 0,
               lwd = 1, box.lwd = 2, cex.txt = 0.8,
               box.type = "square", box.prop = 0.5, arr.type = "triangle",
               arr.pos = 0.4, shadow.size = 0.01, prefix = "f",
               main = "Energy network")














