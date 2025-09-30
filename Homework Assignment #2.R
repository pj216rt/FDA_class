#hw2
library(tidyverse)

#problem 2
#program to compute and display geodesic path between any two points
#looking at equation 2.8 in Anuj's book

#function to build test functions
f_t <- function(t, a, b){
  part1 <- beta(a,b)
}

problem2 <- function(p, q, n_steps = 1000){
  #need to rescale the points
  p <- p / sqrt(sum(p^2))
  q <- q / sqrt(sum(q^2))
  
  #for unit vectors a,b in Hilbert space, \theta = arccos({x,y})
  angle <- acos(sum(p * q))
  
  #create matrix to hold the path
  #each row will be the path at t = n_step_i, and the length will
  #correspond to the length of the vector
  path <- matrix(NA, nrow = n_steps, ncol = length(p))
  
  #computing path, assuming that the points are not antipodal, since then any path will do.
  #I think we only compute the angle once?
  for (i in 1:n_steps) {
    #compute time, starting at 0
    t <- (i-1) / (n_steps-1)
    
    #the first term
    term1 <- sin(angle*(1-t))*p
    term2 <- sin(t*angle)*q
    
    #path at time t
    combined <- (term1+term2)/sin(angle)
    path[i, ] <- combined
  }
  
  #return the path
  return(path)
}
