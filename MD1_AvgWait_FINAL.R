###############################################################################
### This computes the average waiting time for Class-1 and Class-2 customers in 
#   an M/D/1 Delayed APQ.
### Prepared by Blair Bilodeau on December 15, 2019
###############################################################################

#########
## Function to compute the average waiting time
Avg_MD1 <- function(lam1, lam2, b, d){
  
  # constants
  mu <- 1 # fixed
  ell <- d
  lam <- lam1 + lam2
  rho1 <- lam1 / mu
  rho2 <- lam2 / mu
  rho <- rho1 + rho2
  lam1A <- lam1 * (1-b)
  rho1A <- lam1A / mu
  
  delta <- (lam1*b) / (mu - lam1A) / (mu - lam1)
  
  # function to get pi_i vals
  # only works up to about i=20
  pi.i.exact <- function(i){
    if (i == 0) {
      val <- 1
    } else if (i == 1){
      val <- exp(rho)-1
    } else{
      k_vec <- 1:(i-1)
      k_rho_vec <- k_vec * rho
      ik_vec <- i - k_vec
      fac_vec <- unlist(lapply(ik_vec, factorial))
      summand_vec <- (-1)^k_vec * exp(k_rho_vec) * (k_rho_vec)^(ik_vec - 1) * (i - k_vec*(1-rho)) / fac_vec
      val <- exp(i*rho) + (-1)^i * sum(summand_vec)
    }
    return((1-rho)*val)
  } 
  pi.i.exact.vec <- Vectorize(pi.i.exact)
  
  # compute pi vec
  # can change the number of options you allow
  samp <- pi.i.exact.vec(0:20) 
  
  # can change these limits, work for (0.05,0.95)
  decay.rate <- 1/uniroot(function(sigma) exp(rho*sigma) - sigma*exp(rho), c(1.001,1000), tol=1e-10)$root 
  
  # find samples decay rate
  samp.decay.rate <- samp[-1] / head(samp,-1)
  
  # check how it's performing (okay, only worry about rho > 0.5)
  # this tolerance can get higher for each rho
  samp.performance <- abs(decay.rate - samp.decay.rate) < 1e-08
  
  # good indices
  samp.indices <- which(samp.performance)
  
  if (length(samp.indices) == 0) print("Bad pi estimation!")
  
  # last index
  samp.I <- max(samp.indices) # last index
  samp.prelim <- samp[1:samp.I]
  prelim.sum <- sum(samp.prelim)
  
  x <- function(eps) log(1 - (1-decay.rate) * (1 - eps - prelim.sum) / (decay.rate * samp[samp.I])) / log(decay.rate)
  pi.vec <- c(samp.prelim, samp[samp.I] * (decay.rate ^ (1:ceiling(x(1e-8)))))
  I <- length(pi.vec)
  
  ###
  ### Term 1
  ###
  
  # ell is fixed (d= ell / mu)
  # j is the varying param
  # using 1/15! approx 0 to 1e-13 precision
  
  fj_1 <- function(j){
    sum0 <- 0
    if (2 <= (j+d)){
      for (k in 2:(j+d)){
        sum1 <- 0
        for (a in 0:(j+d-k)){
          t1 <- (-1)^a * d^{j+d-k-a} / factorial(j+d-k-a) / factorial(a)
          sum2 <- 0
          for (n in max(0,k-I+1):(k-1)){
            t2 <- pi.vec[k-n+1] * lam1^(j+d+n-k) / factorial(n) / (n+a+1)
            sum2 <- sum2 + t2
          }
          sum1 <- sum1 + t1 * sum2
        } 
        sum0 <- sum0 + sum1
      }
    }
    exp(-lam1*d) * sum0
  }
  fj_1.vec <- Vectorize(fj_1)
  
  ###
  ### Term 2
  ###
  
  fj_2 <- function(j){
    if (2 > d){
      0
    } else {
      sum0 <- 0
      for (k in 2:d){
        sum1 <- 0
        for (m in k:d){
          t1 <- (m-1)^(m-k) * (k-1) / factorial(m-k) / (m-1)
          sum2 <- 0
          for (a in 0:(j+d-m)){
            t2 <- (-1)^a * (d-m+1)^{j+d-m-a} / factorial(j+d-m-a) / factorial(a)
            sum3 <- 0
            for (n in max(0,k-I+1):(k-1)){
              t3 <- pi.vec[k-n+1] * lam1^(j+d+n-k) / factorial(n) / (n+a+1)
              sum3 <- sum3 + t3
            }
            sum2 <- sum2 + t2 * sum3
          } 
          sum1 <- sum1 + t1 * sum2 
        }
        sum0 <- sum0 + sum1
      }
      exp(-lam1*d) * sum0
    }
  }
  fj_2.vec <- Vectorize(fj_2)
  
  ###
  ### Term 1 (r)
  ###
  
  fj_1r <- function(j){
    sum0 <- 0
    if (2 <= (j+d)){
      for (k in 2:(j+d)){
        sum1 <- 0
        for (a in 0:(j+d-k)){
          t1 <- (-1)^a * d^{j+d-k-a} / factorial(j+d-k-a) / factorial(a)
          sum2 <- 0
          for (n in max(0,k-I+1):(k-1)){
            t2 <- pi.vec[k-n+1] * lam1^(j+d+n-k) / factorial(n) / (n+a+2)
            sum2 <- sum2 + t2
          }
          sum1 <- sum1 + t1 * sum2
        } 
        sum0 <- sum0 + sum1
      }
    }
    exp(-lam1*d) * sum0
  }
  fj_1r.vec <- Vectorize(fj_1r)
  
  ###
  ### Term 2 (r)
  ###
  
  fj_2r <- function(j){
    if (2 > d){
      0
    } else {
      sum0 <- 0
      for (k in 2:d){
        sum1 <- 0
        for (m in k:d){
          t1 <- (m-1)^(m-k) * (k-1) / factorial(m-k) / (m-1)
          sum2 <- 0
          for (a in 0:(j+d-m)){
            t2 <- (-1)^a * (d-m+1)^{j+d-m-a} / factorial(j+d-m-a) / factorial(a)
            sum3 <- 0
            for (n in max(0,k-I+1):(k-1)){
              t3 <- pi.vec[k-n+1] * lam1^(j+d+n-k) / factorial(n) / (n+a+2)
              sum3 <- sum3 + t3
            }
            sum2 <- sum2 + t2 * sum3
          } 
          sum1 <- sum1 + t1 * sum2 
        }
        sum0 <- sum0 + sum1
      }
      exp(-lam1*d) * sum0
    }
  }
  fj_2r.vec <- Vectorize(fj_2r)
  
  ###
  ### Putting them together
  ###
  
  if (d==0){
    diff <- delta * 0.5 * (rho + rho^2 / (1-rho))
    
  } else {
    j_summand_1 <- function(j) (j-1) * fj_1.vec(j) + fj_1r.vec(j)
    j_summand_2 <- function(j) (j-1) * fj_2.vec(j) + fj_2r.vec(j)
    #j <- 1:(I-d-1) # want this to be 50 but for speed go with 15
    j <- 1:I
    diff <- delta * (sum(j_summand_1(j)) - sum(j_summand_2(j)))
  }
  
  # compute expectations
  EW_FCFS <- rho / (2 * mu * (1-rho))
  EW_NP2 <- lam / (2 * mu^2 * (1-rho) * (1-rho1))
  EW_NP1 <- lam / (2 * mu^2 * (1-rho1))
  EW_APQ2 <- lam / (mu^2 * (1-rho) * (1-rho1A)) # wrong
  EW_APQ1 <- lam / mu^2 / (1-rho) - rho2 * EW_APQ2 * (1-b) # wrong
  EW_D2 <- EW_NP2 - diff
  EW_D1 <- (rho * EW_FCFS - rho2 * EW_D2) / rho1
  
  return(data.frame('Non-Preemptive Class 1'=EW_NP1,
                    'Non-Preemptive Class 2'=EW_NP2,
                    'APQ Class 1'=EW_APQ1,
                    'APQ Class 2'=EW_APQ2,
                    'Delayed APQ Class 1'=EW_D1,
                    'Delayed APQ Class 2'=EW_D2,
                    'First Come First Served'=EW_FCFS))
}

## Vectorize the arguments to allow for efficient plotting
Avg_MD1.lam1 <- Vectorize(Avg_MD1, vectorize.args = 'lam1')
Avg_MD1.b <- Vectorize(Avg_MD1, vectorize.args = 'b')
Avg_MD1.d <- Vectorize(Avg_MD1, vectorize.args = 'd')
