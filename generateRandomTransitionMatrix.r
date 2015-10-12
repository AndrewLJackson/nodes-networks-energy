generateRandomTransitionMatrix <- function(n.nodes, p.connect, 
                                           upper.tri = FALSE, a.upper, a.lower,
                                           e.upper, e.lower) 
  
{
  
  # Random connections between nodes
  B <- matrix(rbinom(n.nodes^2, 1, p.connect), ncol = n.nodes, nrow = n.nodes)
  
  # only fill in the lower triangle of B
  # nodes in the lower triangle consume energy from the 
  # node specified by each column. That is, if a[2,1] <-1
  # then node2 gains one unit from node1 and so node1 has to
  # have a negative value associated in its diagonal.
  if (!upper.tri){B[upper.tri(B)] <- 0}
  
  # nodes are not connected to themselves
  B[diag(B)] <- 0
  
  # specify the weighted interaction matrix A
  # this matrix determines the rates of energy flow between nodes.
  # We will let energy production stay constant at 1, and scale the others
  # below this amount.
  
  aa <- matrix(runif(n.nodes^2, a.lower, a.upper), 
               ncol = n.nodes, nrow = n.nodes)
  
  # Create the transition matrix a for each valid connection specified by 
  # the binary connection matrix B
  a <- aa * B
  
  
  # Each node loses energy. This simulates energy loss during transfer,
  # and also potentially acts as a feedback of energy to the environment,
  # but that would need some re-coding to achieve.
  e <- runif(n.nodes, e.lower, e.upper)
  
  
  # calculate the diagonal values representing sums of losses to connected nodes
  # in addition to loss to the environment
  diag(a) <- -colSums(a) - e
  
  return(a)
  
  
}
