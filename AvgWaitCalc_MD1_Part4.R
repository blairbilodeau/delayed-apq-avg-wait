
#################
# Testing Phase 2
#################

mu <- 1 # fixed
lam1 <- 0.6
lam2 <- 0.3
b <- 1

ell <- 3
d <- ell / mu
lam <- lam1 + lam2
rho1 <- lam1 / mu
rho2 <- lam2 / mu
rho <- rho1 + rho2
lam1A <- lam1 * (1-b)
rho1A <- lam1A / mu

delta <- (lam1*b) / (mu - lam1A) / (mu - lam1)

#################
# Stationary Probabilities
#################

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

decay.rate <- function(rho){
  # these limits are maybe not great, but work for rho in (0.05,0.95)
  1/uniroot(function(sigma) exp(rho*sigma) - sigma*exp(rho), c(1.001,1000), tol=1e-10)$root
}
decay.rate.vec <- Vectorize(decay.rate)
decay.rate.vec(seq(0.05,0.95,0.05))

### compute pi vec

# can change the number of options you allow
samp <- pi.i.exact.vec(0:20) 

# can change these limits, work for (0.05,0.95)
decay.rate <- 1/uniroot(function(sigma) exp(rho*sigma) - sigma*exp(rho), c(1.001,1000), tol=1e-10)$root 

# find samples decay rate
samp.decay.rate <- samp[-1] / head(samp,-1)

# check how it's performing (okay, only worry about rho > 0.5)
samp.performance <- abs(decay.rate - samp.decay.rate) < 1e-07

# good indices
samp.indices <- which(samp.performance)

if (length(samp.indices) == 0) print("Bad pi estimation!")

# last index
samp.I <- max(samp.indices) # last index
samp.prelim <- samp[1:samp.I]
plot(0:(samp.I-1),samp.prelim)
prelim.sum <- sum(samp.prelim)

x <- function(eps) log(1 - (1-decay.rate) * (1 - eps - prelim.sum) / (decay.rate * samp[samp.I])) / log(decay.rate)
plot(seq(1e-5,1-prelim.sum,2e-5), x(seq(1e-5,1-prelim.sum,2e-5)),type='l')

pi.vec <- c(samp.prelim, samp[samp.I] * (decay.rate ^ (1:ceiling(x(1e-8)))))
I <- length(pi.vec)

plot(pi.vec)
abline(v=samp.I)

#################
# Term 1
#################

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

#################
# Term 2
#################

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

#################
# Term 1 (r)
#################

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
fj_1r.vec <- Vectorize(fj_1)

#################
# Term 2 (r)
#################

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
fj_2r.vec <- Vectorize(fj_2)

#################
# Putting them together
#################
if (d==0){
  diff <- delta * 0.5 * (rho + rho^2 / (1-rho))
  
} else {
  j_summand_1 <- function(j) (j-1) * fj_1.vec(j) + fj_1r.vec(j)
  j_summand_2 <- function(j) (j-1) * fj_2.vec(j) + fj_2r.vec(j)
  #j <- 1:(I-d-1) # want this to be 50 but for speed go with 15
  j <- 1:15
  diff <- delta * (sum(j_summand_1(j)) - sum(j_summand_2(j))) / rho
}

# compute expectations
EW_FCFS <- rho / (2 * mu * (1-rho))
EW_NP2 <- lam / (2 * mu^2 * (1-rho) * (1-rho1))
EW_NP1 <- lam / (2 * mu^2 * (1-rho1))
EW_APQ2 <- lam / (mu^2 * (1-rho) * (1-rho1A)) # wrong
EW_APQ1 <- lam / mu^2 / (1-rho) - rho2 * EW_APQ2 * (1-b) # wrong
EW_D2 <- EW_NP2 - diff
EW_D1 <- (rho * EW_FCFS - rho2 * EW_D2) / rho1

#################
# Implementation
#################

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

Avg_MD1.lam1 <- Vectorize(Avg_MD1, vectorize.args = 'lam1')
Avg_MD1.b <- Vectorize(Avg_MD1, vectorize.args = 'b')
Avg_MD1.d <- Vectorize(Avg_MD1, vectorize.args = 'd')

#################
# Plots
#################

Avg_MD1(lam1=0.6, lam2=0.3, b=0.5, d=2)

test <- Avg_MD1.lam1(lam1=seq(0,0.65,0.01), lam2=0.3, b=0.5, d=2)

plot(seq(0,0.65,0.01), test['Non.Preemptive.Class.1',], type='l', col='red')
lines(seq(0,0.65,0.01), test['Non.Preemptive.Class.2',], type='l', col='blue')
plot(seq(0,0.65,0.01), test['Delayed.APQ.Class.1',], type='l', col='green')
lines(seq(0,0.65,0.01), test['Delayed.APQ.Class.2',], type='l', col='purple')


#################
# d Values
#################

# d 3
d_seq <- c(c(0,1),seq(2,30,2))
d_seq <- seq(0,15,1)

Avg_MD1(lam1=0.5, lam2=0.3, b=1, d=0)

adj_d1 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.2, d=d_seq)
adj_d2 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.4, d=d_seq)
adj_d3 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.6, d=d_seq)
adj_d4 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.8, d=d_seq)
adj_d5 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=1, d=d_seq)

plot(d_seq, adj_d1['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("M/D/1 Class 1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')
points(d_seq, adj_d1['Delayed.APQ.Class.1',], pch=20, col='darkblue')
#lines(d_seq, adj_d1['Delayed.APQ.Class.1',], col='darkblue')
#points(0, adj_d1['Delayed.APQ.Class.1',1], pch=20, col='darkblue')
points(d_seq, adj_d2['Delayed.APQ.Class.1',], pch=20, col='blue') 
#lines(d_seq, adj_d2['Delayed.APQ.Class.1',], col='blue') 
#points(0, adj_d2['Delayed.APQ.Class.1',1], pch=20, col='blue')
points(d_seq, adj_d3['Delayed.APQ.Class.1',], pch=20, col='darkmagenta') 
#lines(d_seq, adj_d3['Delayed.APQ.Class.1',], col='darkmagenta') 
#points(0, adj_d3['Delayed.APQ.Class.1',1], pch=20, col='darkmagenta')
points(d_seq, adj_d4['Delayed.APQ.Class.1',], pch=20, col='darkgreen')
#lines(d_seq, adj_d4['Delayed.APQ.Class.1',], col='darkgreen')
#points(0, adj_d4['Delayed.APQ.Class.1',1], pch=20, col='darkgreen')
points(d_seq, adj_d5['Delayed.APQ.Class.1',], pch=20, col='darkorange') 
#lines(d_seq, adj_d5['Delayed.APQ.Class.1',], col='darkorange') 
#points(0, adj_d5['Delayed.APQ.Class.1',1], pch=20, col='darkorange')
lines(d_seq, adj_d1['First.Come.First.Served',], type='l', col='red') 
axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=1.0', 'Delayed APQ, b=0.8', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.2', 'Non-Preemptive'), 
       lty = c(1,NA,NA,NA,NA,NA,1),
       pch = c(NA,20,20,20,20,20,NA), 
       col=c('red', 'darkorange', 'darkgreen', 'darkmagenta', 'blue', 'darkblue', 'black'))

plot(d_seq, adj_d1['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("M/D/1 Class 2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,9), xaxt='n')
points(d_seq, adj_d1['Delayed.APQ.Class.2',], pch=20, col='darkblue')
#points(0, adj_d1['Delayed.APQ.Class.2',1], pch=20, col='darkblue')
points(d_seq, adj_d2['Delayed.APQ.Class.2',], pch=20, col='blue')
#points(0, adj_d2['Delayed.APQ.Class.2',1], pch=20, col='blue')
points(d_seq, adj_d3['Delayed.APQ.Class.2',], pch=20, col='darkmagenta')
#points(0, adj_d3['Delayed.APQ.Class.2',1], pch=20, col='darkmagenta')
points(d_seq, adj_d4['Delayed.APQ.Class.2',], pch=20, col='darkgreen')
#points(0, adj_d4['Delayed.APQ.Class.2',1], pch=20, col='darkgreen')
points(d_seq, adj_d5['Delayed.APQ.Class.2',], pch=20, col='darkorange')
#points(0, adj_d5['Delayed.APQ.Class.2',1], pch=20, col='darkorange')
lines(d_seq, adj_d1['First.Come.First.Served',], type='l', col='red')
axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=1.0', 'Delayed APQ, b=0.8', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.2', 'Non-Preemptive'), 
       lty = c(1,NA,NA,NA,NA,NA,1),
       pch = c(NA,20,20,20,20,20,NA), 
       col=c('red', 'darkorange', 'darkgreen', 'darkmagenta', 'blue', 'darkblue', 'black'))


#################
# b Values
#################

# b 1
b_seq <- seq(0,1,0.02)
adj_b0 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=0)
adj_b1 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=1)
adj_b2 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=2)
adj_b3 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=3)
adj_b4 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=4)

plot(b_seq, adj_b1['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("M/D/1 Class 1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='b', ylab='Expected Waiting Time',
     ylim=c(0,5))
lines(b_seq, adj_b1['Delayed.APQ.Class.1',], type='l', col='darkorange')
lines(b_seq, adj_b2['Delayed.APQ.Class.1',], type='l', col='darkgreen')
lines(b_seq, adj_b3['Delayed.APQ.Class.1',], type='l', col='blue')
lines(b_seq, adj_b4['Delayed.APQ.Class.1',], type='l', col='darkblue')
lines(b_seq, adj_b0['Delayed.APQ.Class.1',], type='l', col='red')
legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=1', 'Delayed APQ, d=2', 'Delayed APQ, d=3', 'Delayed APQ, d=4', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'blue', 'darkblue', 'black'))

plot(b_seq, adj_b1['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("M/D/1 Class 2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='b', ylab='Expected Waiting Time',
     ylim=c(2,5))
lines(b_seq, adj_b1['Delayed.APQ.Class.2',], type='l', col='darkorange')
lines(b_seq, adj_b2['Delayed.APQ.Class.2',], type='l', col='darkgreen')
lines(b_seq, adj_b3['Delayed.APQ.Class.2',], type='l', col='blue')
lines(b_seq, adj_b4['Delayed.APQ.Class.2',], type='l', col='darkblue')
lines(b_seq, adj_b0['Delayed.APQ.Class.2',], type='l', col='red')
legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=1', 'Delayed APQ, d=2', 'Delayed APQ, d=3', 'Delayed APQ, d=4', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'blue', 'darkblue', 'black'))
# legend('topleft', 
#        c('Non-Preemptive', 'Delayed APQ, d=4', 'Delayed APQ, d=3', 'Delayed APQ, d=2', 'Delayed APQ, d=1', 'APQ'), 
#        lty=rep(1,6), 
#        col=c('black', 'darkblue', 'blue', 'darkgreen', 'darkorange', 'red'))