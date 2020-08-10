# ______________________________________________________________________________
# Compute lower bound for low-dimensional subspace matrix based on Menon theorem
# ______________________________________________________________________________

# input:
# ------
# n = sample size
# epsi = distance preservation factor
# prob = 1 - probability of distance preservation

# output:
# -------
# k = Menon lower bound


lbMenon <- function(n, epsi, prob){
  
  beta <- -logb(prob, base=n)
  k <- (16+8*beta)/(epsi^2) * log(n)
  
  return(list(k=k))
}