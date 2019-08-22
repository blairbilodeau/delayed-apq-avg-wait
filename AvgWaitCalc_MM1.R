#################
# Implementation
#################

Avg_MM1 <- function(mu, lam1, lam2, b, d, k){

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
  x <- matrix(0L, nrow=k, ncol=k)
  
  # setup first 2 rows
  x[1,1] <- r - p
  x[2,1] <- q*r*p
  x[2,2] <- r^2 - p^2
  
  # compute lower triangle
  for (row in 3:k){
    x[row,1] <- q * x[row-1,2]
    for (col in 2:(k-1)){
      x[row,col] <- p * x[row-1,col-1] + q * x[row-1,col+1]
    }
    x[row,row] <- r^k - p^k
  }
  
  # test for convergence
  xk <- rowSums(x)
  # plot(1:k, xk)
  if (xk[k] > 1e-10) print('ERROR: x_l did not converge')
  
  # compute poisson probabilities
  pk <- dpois(1:k, nu*d)
  
  # combine terms
  T1 <- sum(pk * xk)
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

Avg_MM1.lam1 <- Vectorize(Avg_MM1, vectorize.args = 'lam1')
Avg_MM1.b <- Vectorize(Avg_MM1, vectorize.args = 'b')
Avg_MM1.d <- Vectorize(Avg_MM1, vectorize.args = 'd')

#################
# Plots
#################

Avg_MM1(mu=1, lam1=0.5, lam2=0.3, b=0.5, d=3, k=1000)

# delayed APQ average wait
# each column is a value of vectorized coord
test <- Avg_MM1.lam1(mu=1, lam1=seq(0,0.65,0.01), lam2=0.3, b=0.5, d=2, k=1000)
test2 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=seq(0,1,0.01), d=2, k=350)

plot(seq(0,0.65,0.01), test['Non.Preemptive.Class.1',], type='l', col='red')
lines(seq(0,0.65,0.01), test['Non.Preemptive.Class.2',], type='l', col='blue')
lines(seq(0,0.65,0.01), test['APQ.Class.1',], type='l', col='black')
lines(seq(0,0.65,0.01), test['APQ.Class.2',], type='l', col='orange')
lines(seq(0,0.65,0.01), test['Delayed.APQ.Class.1',], type='l', col='green')
lines(seq(0,0.65,0.01), test['Delayed.APQ.Class.2',], type='l', col='purple')

plot(seq(0,1,0.01), test2[1,], type='l', col='red', ylim=c(0,10))
lines(seq(0,1,0.01), test2[2,], col='blue')

#################
# Rho Values
#################

# rho 1
rho_seq <- seq(0,0.65,0.01)
adj_rho <- Avg_MM1.lam1(mu=1, lam1=rho_seq, lam2=0.3, b=0.5, d=2, k=1000)

plot(rho_seq, adj_rho['Non.Preemptive.Class.1',], type='l', col='black',
     main="Class 1, b=0.5 and a=1", xlab=expression(rho), ylab='Expected Waiting Time',
     xlim=c(0.05, 0.6), ylim=c(0,2.5), xaxt='n')
lines(rho_seq, adj_rho['Delayed.APQ.Class.1',], type='l', col='blue')
lines(rho_seq, adj_rho['APQ.Class.1',], type='l', col='red')
axis(1, at=seq(0,0.6,0.1), labels=(seq(0,0.6,0.1)+0.3))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

plot(rho_seq, adj_rho['Non.Preemptive.Class.2',], type='l', col='black',
     main="Class 2, b=0.5 and a=1", xlab=expression(rho), ylab='Expected Waiting Time',
     xlim=c(0.05, 0.6), ylim=c(0,2.5), xaxt='n')
lines(rho_seq, adj_rho['Delayed.APQ.Class.2',], type='l', col='blue')
lines(rho_seq, adj_rho['APQ.Class.2',], type='l', col='red')
axis(1, at=seq(0,0.6,0.1), labels=(seq(0,0.6,0.1)+0.3))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

# rho 2
rho_seq <- seq(0,0.65,0.01)
adj_rho <- Avg_MM1.lam1(mu=1, lam1=rho_seq, lam2=0.3, b=0.5, d=1.4, k=1000)

plot(rho_seq, adj_rho['Non.Preemptive.Class.1',], type='l', col='black',
     main="Class 1, b=0.5 and d=1.4", xlab=expression(rho), ylab='Expected Waiting Time',
     xlim=c(0.05, 0.6), ylim=c(0,2.5), xaxt='n')
lines(rho_seq, adj_rho['Delayed.APQ.Class.1',], type='l', col='blue')
lines(rho_seq, adj_rho['APQ.Class.1',], type='l', col='red')
axis(1, at=seq(0,0.6,0.1), labels=(seq(0,0.6,0.1)+0.3))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

plot(rho_seq, adj_rho['Non.Preemptive.Class.2',], type='l', col='black',
     main="Class 2, b=0.5 and a=0.7", xlab=expression(rho), ylab='Expected Waiting Time',
     xlim=c(0.05, 0.6), ylim=c(0,2.5), xaxt='n')
lines(rho_seq, adj_rho['Delayed.APQ.Class.2',], type='l', col='blue')
lines(rho_seq, adj_rho['APQ.Class.2',], type='l', col='red')
axis(1, at=seq(0,0.6,0.1), labels=(seq(0,0.6,0.1)+0.3))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

#################
# d Values
#################

# d 1
d_seq <- seq(0,20,0.2)
adj_d <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.5, d=d_seq, k=1000)

plot(d_seq, adj_d['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("Class 1, ", lambda_1, "=0.5", lambda_2, "=0.3", mu, "=1")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')
lines(d_seq, adj_d['Delayed.APQ.Class.1',], type='l', col='blue')
lines(d_seq, adj_d['APQ.Class.1',], type='l', col='red')
axis(1, at=seq(0,20,2))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

plot(d_seq, adj_d['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("Class 2, ", rho, "=0.8")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(5,10), xaxt='n')
lines(d_seq, adj_d['Delayed.APQ.Class.2',], type='l', col='blue')
lines(d_seq, adj_d['APQ.Class.2',], type='l', col='red')
axis(1, at=seq(0,20,2))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

# d 2
d_seq <- seq(0,20,0.25)
adj_d <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.8, d=d_seq, k=1000)

plot(d_seq, adj_d['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("Class 1, ", rho, "=0.8")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')
lines(d_seq, adj_d['Delayed.APQ.Class.1',], type='l', col='blue')
lines(d_seq, adj_d['APQ.Class.1',], type='l', col='red')
axis(1, at=seq(0,20,2))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

plot(d_seq, adj_d['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("Class 2, ", rho, "=0.8")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,10), xaxt='n')
lines(d_seq, adj_d['Delayed.APQ.Class.2',], type='l', col='blue')
lines(d_seq, adj_d['APQ.Class.2',], type='l', col='red')
axis(1, at=seq(0,20,2))
legend('topleft', 
       c('Non-Preemptive', 'Delayed APQ', 'APQ'), 
       lty=rep(1,3), 
       col=c('black', 'blue', 'red'))

# d 3
d_seq <- seq(0,15,0.3)
adj_d1 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.2, d=d_seq, k=1000)
adj_d2 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.4, d=d_seq, k=1000)
adj_d3 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.6, d=d_seq, k=1000)
adj_d4 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.8, d=d_seq, k=1000)
adj_d5 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=1, d=d_seq, k=1000)

plot(d_seq, adj_d1['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("M/M/1 Class 1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')
lines(d_seq, adj_d1['Delayed.APQ.Class.1',], type='l', col='darkblue')
points(0, adj_d1['APQ.Class.1',1], pch=20, col='darkblue')
lines(d_seq, adj_d2['Delayed.APQ.Class.1',], type='l', col='blue')
points(0, adj_d2['APQ.Class.1',1], pch=20, col='blue')
lines(d_seq, adj_d3['Delayed.APQ.Class.1',], type='l', col='darkmagenta')
points(0, adj_d3['APQ.Class.1',1], pch=20, col='darkmagenta')
lines(d_seq, adj_d4['Delayed.APQ.Class.1',], type='l', col='darkgreen')
points(0, adj_d4['APQ.Class.1',1], pch=20, col='darkgreen')
lines(d_seq, adj_d5['Delayed.APQ.Class.1',], type='l', col='darkorange')
points(0, adj_d5['APQ.Class.1',1], pch=20, col='darkorange')
lines(d_seq, adj_d1['First.Come.First.Served',], type='l', col='red')
axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=1.0', 'Delayed APQ, b=0.8', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.2', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'darkmagenta', 'blue', 'darkblue', 'black'))

plot(d_seq, adj_d1['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("M/M/1 Class 2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='d', ylab='Expected Waiting Time',
     ylim=c(3,9), xaxt='n')
lines(d_seq, adj_d1['Delayed.APQ.Class.2',], type='l', col='darkblue')
points(0, adj_d1['Delayed.APQ.Class.2',1], pch=20, col='darkblue')
lines(d_seq, adj_d2['Delayed.APQ.Class.2',], type='l', col='blue')
points(0, adj_d2['Delayed.APQ.Class.2',1], pch=20, col='blue')
lines(d_seq, adj_d3['Delayed.APQ.Class.2',], type='l', col='darkmagenta')
points(0, adj_d3['Delayed.APQ.Class.2',1], pch=20, col='darkmagenta')
lines(d_seq, adj_d4['Delayed.APQ.Class.2',], type='l', col='darkgreen')
points(0, adj_d4['Delayed.APQ.Class.2',1], pch=20, col='darkgreen')
lines(d_seq, adj_d5['Delayed.APQ.Class.2',], type='l', col='darkorange')
points(0, adj_d5['Delayed.APQ.Class.2',1], pch=20, col='darkorange')
lines(d_seq, adj_d1['First.Come.First.Served',], type='l', col='red')
axis(1, at=seq(0,15,1))
legend('bottomright', 
       c( 'First Come First Served', 'Delayed APQ, b=1.0', 'Delayed APQ, b=0.8', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.2', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'darkmagenta', 'blue', 'darkblue', 'black'))

#################
# b Values
#################

# b 1
b_seq <- seq(0,1,0.01)
adj_b1 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=2, k=500)
adj_b2 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=4, k=500)
adj_b3 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=6, k=500)
adj_b4 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=8, k=500)

plot(b_seq, adj_b1['Non.Preemptive.Class.1',], type='l', col='black',
     main=expression(paste("M/M/1 Class 1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='b', ylab='Expected Waiting Time',
     ylim=c(0,5))
lines(b_seq, adj_b1['Delayed.APQ.Class.1',], type='l', col='darkorange')
lines(b_seq, adj_b2['Delayed.APQ.Class.1',], type='l', col='darkgreen')
lines(b_seq, adj_b3['Delayed.APQ.Class.1',], type='l', col='blue')
lines(b_seq, adj_b4['Delayed.APQ.Class.1',], type='l', col='darkblue')
lines(b_seq, adj_b1['APQ.Class.1',], type='l', col='red')
legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=1', 'Delayed APQ, d=2', 'Delayed APQ, d=3', 'Delayed APQ, d=4', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'blue', 'darkblue', 'black'))

plot(b_seq, adj_b1['Non.Preemptive.Class.2',], type='l', col='black',
     main=expression(paste("M/M/1 Class 2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='b', ylab='Expected Waiting Time',
     ylim=c(3,9))
lines(b_seq, adj_b1['Delayed.APQ.Class.2',], type='l', col='darkorange')
lines(b_seq, adj_b2['Delayed.APQ.Class.2',], type='l', col='darkgreen')
lines(b_seq, adj_b3['Delayed.APQ.Class.2',], type='l', col='blue')
lines(b_seq, adj_b4['Delayed.APQ.Class.2',], type='l', col='darkblue')
lines(b_seq, adj_b1['APQ.Class.2',], type='l', col='red')
legend('bottomleft', 
       c( 'APQ', 'Delayed APQ, d=1', 'Delayed APQ, d=2', 'Delayed APQ, d=3', 'Delayed APQ, d=4', 'Non-Preemptive'), 
       lty=rep(1,7), 
       col=c('red', 'darkorange', 'darkgreen', 'blue', 'darkblue', 'black'))
# legend('topleft', 
#        c('Non-Preemptive', 'Delayed APQ, d=4', 'Delayed APQ, d=3', 'Delayed APQ, d=2', 'Delayed APQ, d=1', 'APQ'), 
#        lty=rep(1,6), 
#        col=c('black', 'darkblue', 'blue', 'darkgreen', 'darkorange', 'red'))

