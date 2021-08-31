###############################################################################
### This computes the average waiting time for Class-1 and Class-2 customers in 
#   an M/M/1 Delayed APQ.
###############################################################################

#########
## Function to compute the average waiting time
## K is the truncation size of the infinite sum
Avg_MM1 <- function(mu, lam1, lam2, b, d, K){
  
  # constants
  a <- b*d
  lam <- lam1 + lam2
  rho1 <- lam1 / mu
  rho2 <- lam2 / mu
  rho <- rho1 + rho2
  lam1A <- lam1 * (1-b)
  rho1A <- lam1A / mu 
  
  # transition matrix
  nu <- mu + lam1
  q <- mu / nu
  p <- lam1 / nu
  r <- p + q*rho^2
  delta <- (lam1 * b) / (mu - lam1A) / (mu - lam1)
  
  # numerical params
  x <- matrix(0L, nrow=K, ncol=K)
  
  # setup first 2 rows
  x[1,1] <- r - p
  x[2,1] <- q*r*p
  x[2,2] <- r^2 - p^2
  x[3,1] <- q * x[2,2]
  x[3,2] <- p * x[2,1] + q * r^3
  x[3,3] <- r^3 - p^3
  
  # compute lower triangle x[k,ell]
  # for (row in 3:K){
  #   x[row,1] <- q * x[row-1,2]
  #   for (col in 2:(K-1)){
  #     x[row,col] <- p * x[row-1,col-1] + q * x[row-1,col+1]
  #   }
  #   x[row,row] <- r^row - p^row
  # }
  for (row in 4:K){
    x[row,1] <- q * x[row-1,2]
    for (col in 2:(row-2)){
      x[row,col] <- p * x[row-1,col-1] + q * x[row-1,col+1]
    }
    x[row,row-1] <- p * x[row-1,row-2] + q * r^(row-1)
    x[row,row] <- r^row - p^row
  }
  
  # multiply with j
  J <- matrix(1:K,nrow=K,ncol=1)
  xK <- x %*% J
  if (xK[K] > 1e-10) print('ERROR: x_l did not converge')
  
  # compute poisson probabilities
  pK <- dpois(1:K, nu*d)
  
  # combine terms
  T1 <- sum(pK * xK)
  T2 <- rho * exp(nu*d*(r-1)) * (1/(1-rho) + r*nu*d) 
  diff <- delta * ((1-rho)*T1 + T2)
  
  # compute expectations
  EW_NP2 <- lam / (mu^2 * (1-rho) * (1-rho1))
  EW_NP1 <- lam / (mu^2 * (1-rho1))
  EW_APQ2 <- lam / (mu^2 * (1-rho) * (1-rho1A))
  EW_APQ1 <- lam / mu^2 / (1-rho) - rho2 * EW_APQ2 * (1-b)
  EW_D2 <- EW_NP2 - diff
  EW_D1 <- (rho^2 / (mu - lam) - rho2 * EW_D2) / rho1
  EW_FCFS <- (rho / mu) / (1-rho)
  
  return(data.frame('Non-Preemptive Class 1'=EW_NP1,
                    'Non-Preemptive Class 2'=EW_NP2,
                    'APQ Class 1'=EW_APQ1,
                    'APQ Class 2'=EW_APQ2,
                    'Delayed APQ Class 1'=EW_D1,
                    'Delayed APQ Class 2'=EW_D2,
                    'First Come First Served'=EW_FCFS))
}

## Vectorize the arguments to allow for efficient plotting
Avg_MM1.lam1 <- Vectorize(Avg_MM1, vectorize.args = 'lam1')
Avg_MM1.b <- Vectorize(Avg_MM1, vectorize.args = 'b')
Avg_MM1.d <- Vectorize(Avg_MM1, vectorize.args = 'd')




