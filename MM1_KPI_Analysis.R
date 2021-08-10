###############################################################################
### This provides parametrization optimization for various KPI levels.
###############################################################################

wdir <- '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/'

code.path <- paste0(wdir, 'code/')
source(paste0(code.path,"MM1_AvgWait.R"))
source(paste0(code.path,"MM1_Wait.R"))
source(paste0(code.path,"MM1_Class1_Approx.R"))

plot.path <- paste0(wdir, 'manuscript/plots/')

##################
## NPQ Boundary
##################

### for fixed lam1 find maximum lam2 such that Class-2 KPI can be achieved with NPQ
npq.lam1 <- function(kpi_level, kpi_prob, lam1, N=12, tol=1e-6){
  npq_diff <- function(lam2) {NP_CDF_class2(mu=1, lam1, lam2, N)(kpi_level) - kpi_prob}
  npq_diff <- Vectorize(npq_diff)
  
  err_try <- function(){
    ans <- uniroot(npq_diff, c(tol/10,1-lam1-tol/10), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(npq_diff(tol/10) < 0){
      return(0)
    } else {
      return(1-lam1)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
  
}
npq.lam1 <- Vectorize(npq.lam1, vectorize.args = 'lam1')

### for fixed lam1 find maximum lam2 such that Class-1 KPI can be achieved with NPQ
npq.lam1.class1 <- function(kpi_level, kpi_prob, lam1, tol=1e-6){
  npq_diff <- function(lam2) {NP_CDF_class1(mu=1, lam1, lam2)(kpi_level) - kpi_prob}
  npq_diff <- Vectorize(npq_diff)
  
  err_try <- function(){
    ans <- uniroot(npq_diff, c(tol/10,1-lam1-tol/10), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance -- returning anyways')
      return(ans$root) # can be dangerous behaviour, but estim.prec can be really big if the function is "too easy"
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(npq_diff(tol/10) < 0){
      return(0)
    } else {
      return(1-lam1)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
  
}
npq.lam1.class1 <- Vectorize(npq.lam1.class1, vectorize.args = 'lam1')

### for fixed lam2 find maximum lam1 such that KPI can be achieved with NPQ
npq.lam2 <- function(kpi_level, kpi_prob, lam2, N=12, tol=1e-6){
  npq_diff <- function(lam1) {NP_CDF_class2(mu=1, lam1, lam2, N)(kpi_level) - kpi_prob}
  npq_diff <- Vectorize(npq_diff)
  
  err_try <- function(){
    ans <- uniroot(npq_diff, c(tol/10,1-lam2-tol/10), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(npq_diff(tol/10) < 0){
      return(0)
    } else {
      return(1-lam2)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
  
}
npq.lam2 <- Vectorize(npq.lam2, vectorize.args = 'lam2')

##################
## FCFS Boundary
##################
# https://www.win.tue.nl/~iadan/que/h4.pdf
# Eq(14), page 6
# P(S+W > t) = exp(-(mu-lam)*t) ... Exp(mu-lam)
# P(W > t) = rho*exp(-mu*(1-rho)*t)

### for fixed lam1 find minimum lam2 such that KPI can't be achieved with FCFS
fcfs.lam1 <- function(kpi_level, kpi_prob, lam1, tol=1e-6){
  fcfs_diff <- function(lam2) {1 - (lam1 + lam2) * exp(-(1-lam1-lam2)*kpi_level) - kpi_prob}
  fcfs_diff <- Vectorize(fcfs_diff)
  
  err_try <- function(){
    ans <- uniroot(fcfs_diff, c(tol/10,1-lam1-tol/10), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(fcfs_diff(tol/10) < 0){
      return(0)
    } else {
      return(1-lam1)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
  
}
fcfs.lam1 <- Vectorize(fcfs.lam1, vectorize.args = 'lam1')

### for fixed lam2 find maximum lam1 such that KPI can't be achieved with FCFS
fcfs.lam2 <- function(kpi_level, kpi_prob, lam2, tol=1e-6){
  fcfs_diff <- function(lam1) {1 - (lam1 + lam2) * exp(-(1-lam1-lam2)*kpi_level) - kpi_prob}
  fcfs_diff <- Vectorize(fcfs_diff)
  
  err_try <- function(){
    ans <- uniroot(fcfs_diff, c(tol/10,1-lam2-tol/10), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(fcfs_diff(tol/10) < 0){
      return(0)
    } else {
      return(1-lam2)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
  
}
fcfs.lam2 <- Vectorize(fcfs.lam2, vectorize.args = 'lam2')

fcfs.cdf <- function(t, lam1, lam2){
  1 - (lam1 + lam2) * exp(-(1-lam1-lam2)*t)
} 

##################
## feasible region plotted together
##################

tol <- 1e-5
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

####
# Class-2

####
# Plot of NPQ lower and FCFS upper boundary

##
# 1 hour

fcfs1_0.95 <- fcfs.lam1(kpi_level=4, kpi_prob=0.95, lam1)
fcfs1_0.85 <- fcfs.lam1(kpi_level=4, kpi_prob=0.85, lam1)
fcfs1_0.6 <- fcfs.lam1(kpi_level=4, kpi_prob=0.6, lam1)

npq1_0.95 <- npq.lam1(kpi_level=4, kpi_prob=0.95, lam1)
npq1_0.85 <- npq.lam1(kpi_level=4, kpi_prob=0.85, lam1)
npq1_0.6 <- npq.lam1(kpi_level=4, kpi_prob=0.6, lam1)

pdf(paste0(plot.path,'MM1_feasible_bound_4.pdf'))

plot(lam1, fcfs1_0.95, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-2)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.95, type='l')
segments(0,fcfs1_0.95[1],0,1)

lines(lam1, npq1_0.85, type='l', lty=2)
lines(lam1, fcfs1_0.85, type='l', lty=2)

lines(lam1, npq1_0.6, type='l', lty=3, cex=0.5)
lines(lam1, fcfs1_0.6, type='l', lty=3, cex=0.5)

legend('topright', 
       c(expression(paste('P(W'[2],' < 4',mu,')=0.6')), expression(paste('P(W'[2],' < 4',mu,')=0.85')), expression(paste('P(W'[2],' < 4',mu,')=0.95'))), 
       lty=c(3,2,1),
       cex=1.5)

dev.off()

##
# 1/2 hour

fcfs1_0.95 <- fcfs.lam1(kpi_level=2, kpi_prob=0.95, lam1)
fcfs1_0.85 <- fcfs.lam1(kpi_level=2, kpi_prob=0.85, lam1)
fcfs1_0.6 <- fcfs.lam1(kpi_level=2, kpi_prob=0.6, lam1)

npq1_0.95 <- npq.lam1(kpi_level=2, kpi_prob=0.95, lam1)
npq1_0.85 <- npq.lam1(kpi_level=2, kpi_prob=0.85, lam1)
npq1_0.6 <- npq.lam1(kpi_level=2, kpi_prob=0.6, lam1)

pdf(paste0(plot.path,'MM1_feasible_bound_2.pdf'))

plot(lam1, fcfs1_0.95, type='l',
     main="M/M/1 FCFS 1/2 Hour Feasible Barrier (Class-2)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.95, type='l')
segments(0,fcfs1_0.95[1],0,1)

lines(lam1, npq1_0.85, type='l', lty=2)
lines(lam1, fcfs1_0.85, type='l', lty=2)

lines(lam1, npq1_0.6, type='l', lty=3, cex=0.5)
lines(lam1, fcfs1_0.6, type='l', lty=3, cex=0.5)

legend('topright', 
       c(expression(paste('P(W'[2],' < 2',mu,')=0.6')), expression(paste('P(W'[2],' < 2',mu,')=0.85')), expression(paste('P(W'[2],' < 2',mu,')=0.95'))), 
       lty=c(3,2,1),
       cex=1.5)

dev.off()

####
# Class-1

####
# Plot of NPQ upper and FCFS lower boundary

##
# 1 hour

fcfs1_0.90 <- fcfs.lam1(kpi_level=4, kpi_prob=0.90, lam1)
fcfs1_0.75 <- fcfs.lam1(kpi_level=4, kpi_prob=0.75, lam1)
fcfs1_0.50 <- fcfs.lam1(kpi_level=4, kpi_prob=0.50, lam1)

npq1_0.90 <- npq.lam1.class1(kpi_level=4, kpi_prob=0.90, lam1)
npq1_0.75 <- npq.lam1.class1(kpi_level=4, kpi_prob=0.75, lam1)
npq1_0.50 <- npq.lam1.class1(kpi_level=4, kpi_prob=0.50, lam1)

# these regions overlap, so for the paper we plot them separately

pdf(paste0(plot.path,'MM1_feasible_bound_4_class1_0.90.pdf'))

plot(lam1, fcfs1_0.90, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.90, type='l')
segments(0,fcfs1_0.90[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 4',mu,')=0.90')),
       lty=1,
       cex=1.5)

dev.off()

pdf(paste0(plot.path,'MM1_feasible_bound_4_class1_0.75.pdf'))

plot(lam1, fcfs1_0.75, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.75, type='l')
segments(0,fcfs1_0.75[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 4',mu,')=0.75')),
       lty=1,
       cex=1.5)

dev.off()

pdf(paste0(plot.path,'MM1_feasible_bound_4_class1_0.50.pdf'))

plot(lam1, fcfs1_0.50, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.50, type='l')
segments(0,fcfs1_0.50[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 4',mu,')=0.50')),
       lty=1,
       cex=1.5)

dev.off()

##
# 1/2 hour

fcfs1_0.90 <- fcfs.lam1(kpi_level=2, kpi_prob=0.90, lam1)
fcfs1_0.75 <- fcfs.lam1(kpi_level=2, kpi_prob=0.75, lam1)
fcfs1_0.50 <- fcfs.lam1(kpi_level=2, kpi_prob=0.50, lam1)

npq1_0.90 <- npq.lam1.class1(kpi_level=2, kpi_prob=0.90, lam1)
npq1_0.75 <- npq.lam1.class1(kpi_level=2, kpi_prob=0.75, lam1)
npq1_0.50 <- npq.lam1.class1(kpi_level=2, kpi_prob=0.50, lam1)

# these regions overlap, so for the paper we plot them separately

pdf(paste0(plot.path,'MM1_feasible_bound_2_class1_0.90.pdf'))

plot(lam1, fcfs1_0.90, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.90, type='l')
segments(0,fcfs1_0.90[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 2',mu,')=0.90')),
       lty=1,
       cex=1.5)

dev.off()

pdf(paste0(plot.path,'MM1_feasible_bound_2_class1_0.75.pdf'))

plot(lam1, fcfs1_0.75, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.75, type='l')
segments(0,fcfs1_0.75[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 2',mu,')=0.75')),
       lty=1,
       cex=1.5)

dev.off()

pdf(paste0(plot.path,'MM1_feasible_bound_2_class1_0.50.pdf'))

plot(lam1, fcfs1_0.50, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.50, type='l')
segments(0,fcfs1_0.50[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 2',mu,')=0.50')),
       lty=1,
       cex=1.5)

dev.off()

####################
## optimizing b and d in the case where FCFS possible and NPQ too penalizing
####################

### b and F2(x) increase together
### b and class-1 expectation increase together
### d goes down means class-1 expectation increases
### find b where F2(x) = p (smallest b that KPI is met)
optim.b <- function(kpi_level, kpi_prob, lam1, lam2, d, K=500, N=12, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(b) {DP_CDF_class2(mu, lam1, lam2, b, d, K, N)(kpi_level) - kpi_prob}
  cdf_diff <- Vectorize(cdf_diff)
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol,1-tol), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(cdf_diff(0.5) < 0){
      return(1)
    } else {
      return(0)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
}
optim.b <- Vectorize(optim.b, vectorize.args='d')

### b and F2(x) increase together
### b and class-1 waiting time decrease together
### find b where F1(x) = p (largest b that KPI is met)
optim.b.class1 <- function(kpi_level, kpi_prob, lam1, lam2, d, K=500, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(b) {DP_CDF_approx_class1(mu, lam1, lam2, b, d, K)(kpi_level) - kpi_prob}
  cdf_diff <- Vectorize(cdf_diff)
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol,1-tol), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(cdf_diff(0.5) < 0){
      return(0)
    } else {
      return(1)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
}
optim.b.class1 <- Vectorize(optim.b.class1, vectorize.args='d')

### find d where F2(x) = p (largest d that KPI is met)
optim.d <- function(kpi_level, kpi_prob, lam1, lam2, b, max_d=25, K=500, N=12, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(d) {DP_CDF_class2(mu, lam1, lam2, b, d, K, N)(kpi_level) - kpi_prob}
  cdf_diff <- Vectorize(cdf_diff)
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol,max_d-tol), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(cdf_diff(1) < 0){
      return(-1)
    } else {
      return(0)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
}
optim.d <- Vectorize(optim.d, vectorize.args='b')

### find d where F1(x) = p (smallest d that KPI is met)
optim.d.class1 <- function(kpi_level, kpi_prob, lam1, lam2, b, max_d=25, K=500, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(d) {DP_CDF_approx_class1(mu, lam1, lam2, b, d, K)(kpi_level) - kpi_prob}
  cdf_diff <- Vectorize(cdf_diff)
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol,max_d-tol), tol=tol)
    if (ans$estim.prec > tol){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(cdf_diff(1) < 0){
      return(-1)
    } else {
      return(0)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
}
optim.d.class1 <- Vectorize(optim.d.class1, vectorize.args='b')

### find d where there exists b such that F2(x) = p (largest d such that KPI can be met)
optim.db <- function(kpi_level, kpi_prob, lam1, lam2, max_d, K=500, N=12, tol=1e-4){
  
  optim_diff <- function(d) {optim.b(kpi_level, kpi_prob, lam1, lam2, d, K, N, tol) - (1-tol/10)}
  optim_diff <- Vectorize(optim_diff)
  
  err_try <- function(){
    ans <- uniroot(optim_diff, c(0,max_d), tol=tol)
    if (ans$estim.prec > 0.01){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(optim_diff(1) < 0){
      return(max_d)
    } else {
      return(0)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
} 

### find d where there exists b<1 such that F1(x) > p (biggest d such that KPI can't be met with b=1)
optim.db.class1 <- function(kpi_level, kpi_prob, lam1, lam2, max_d=25, K=500, tol=1e-4){
  
  optim_diff <- function(d) {optim.b.class1(kpi_level, kpi_prob, lam1, lam2, d, K, tol) - (1-tol/10)}
  optim_diff <- Vectorize(optim_diff)
  
  err_try <- function(){
    ans <- uniroot(optim_diff, c(0,max_d), tol=tol)
    if (ans$estim.prec > 0.01){
      print('not enough tolerance')
    } else {
      return(ans$root)
    }
  }
  
  err_catch <- function(){
    if(optim_diff(1) < 0){
      return(max_d)
    } else {
      return(0)
    }  
  }
  
  tryCatch(err_try(), error=function(cond){return(err_catch())})
} 

#####################
# Plot pairs of d and Class-1 expected waiting time at b*(d)
# Also plot d vs b*(d)
#####################

# helper to identify lam1,lam2 pairs, can also be done visually
tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

npq1_0.85 <- npq.lam1(kpi_level=4, kpi_prob=0.85, lam1)
fcfs1_0.85 <- fcfs.lam1(kpi_level=4, kpi_prob=0.85, lam1)

# lam1 = 0.55, lam2 = 0.05
d_max <- optim.db(4,0.85,0.55,0.05,1,K=1000,N=12)

d_vals <- c(seq(1e-6,d_max,0.025))
b_vals <- optim.b(4,0.85,0.55,0.05,d_vals,K=1000)

pdf(paste0(plot.path,'MM1_kpi_example_A_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-2 Constraint,"), paste("where ", lambda[1], "=0.55, ", lambda[2], "=0.05, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.55, 0.05, b_vals[i], d_vals[i], K=1000)[,'Delayed.APQ.Class.1']
}

pdf(paste0(plot.path,'MM1_kpi_example_A.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.55, ", lambda[2], "=0.05, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(1.476,1.48))

dev.off()

# lam1 = 0.3, lam2 = 0.3
d_max <- optim.db(4,0.85,0.3,0.3,5,K=1000,N=12)

d_vals <- c(seq(1e-6,d_max,0.1))
b_vals <- optim.b(4,0.85,0.3,0.3,d_vals,K=1000)

pdf(paste0(plot.path,'MM1_kpi_example_B_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-2 Constraint,"), paste("where ", lambda[1], "=0.3, ", lambda[2], "=0.3, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.3, 0.3, b_vals[i], d_vals[i], K=500)[,'Delayed.APQ.Class.1']
}

pdf(paste0(plot.path,'MM1_kpi_example_B.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.3, ", lambda[2], "=0.3, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(1.24,1.3))

dev.off()

# lam1 = 0.15, lam2 = 0.45
d_max <- optim.db(4,0.85,0.15,0.45,5,K=500,N=12)

d_vals <- c(seq(1e-6,d_max,0.25))
b_vals <- optim.b(4,0.85,0.15,0.45,d_vals,K=500)

pdf(paste0(plot.path,'MM1_kpi_example_C_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-2 Constraint,"), paste("where ", lambda[1], "=0.15, ", lambda[2], "=0.45, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.15, 0.45, b_vals[i], d_vals[i], K=500)[,'Delayed.APQ.Class.1']
}

pdf(paste0(plot.path,'MM1_kpi_example_C.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.15, ", lambda[2], "=0.45, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(0.74,1))

dev.off()

#####################
# Plot pairs of d and Class-2 expected waiting time at b*(d)
# Also plot d vs b*(d)
#####################

# helper to identify lam1,lam2 pairs, can also be done visually
tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

npq1_0.90 <- npq.lam1.class1(kpi_level=2, kpi_prob=0.90, lam1)
fcfs1_0.90 <- fcfs.lam1(kpi_level=2, kpi_prob=0.90, lam1)

npq1_0.90[11]
fcfs1_0.90[11]

plot(lam1, npq1_0.90)
plot(lam1, fcfs1_0.90)

# lam1 = 0.3, lam2 = 0.1
d_max <- optim.db.class1(2,0.90,0.3,0.1,10,K=500)

d_vals <- c(seq(1e-6,d_max,0.35))
b_vals <- optim.b.class1(2,0.90,0.3,0.1,d_vals,K=500)

pdf(paste0(plot.path,'MM1_kpi_class1_example_A_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-1 Constraint,"), paste("where ", lambda[1], "=0.3, ", lambda[2], "=0.1, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.3, 0.1, b_vals[i], d_vals[i], K=1000)[,'Delayed.APQ.Class.2']
}

pdf(paste0(plot.path,'MM1_kpi_class1_example_A.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-2 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.3, ", lambda[2], "=0.1, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(0.935,0.936))

dev.off()

# lam1 = 0.2, lam2 = 0.2
d_max <- optim.db.class1(2,0.90,0.2,0.2,5,K=500)

d_vals <- c(seq(1e-6,d_max,0.1))
b_vals <- optim.b.class1(2,0.90,0.2,0.2,d_vals,K=500)

pdf(paste0(plot.path,'MM1_kpi_class1_example_B_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-1 Constraint,"), paste("where ", lambda[1], "=0.2, ", lambda[2], "=0.2, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.2, 0.2, b_vals[i], d_vals[i], K=1000)[,'Delayed.APQ.Class.2']
}

pdf(paste0(plot.path,'MM1_kpi_class1_example_B.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-2 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.2, ", lambda[2], "=0.2, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(0.756,0.757))

dev.off()

# lam1 = 0.1, lam2 = 0.3
d_max <- optim.db.class1(2,0.90,0.1,0.3,10,K=500)

d_vals <- c(seq(1e-6,d_max,0.05))
b_vals <- optim.b.class1(2,0.90,0.1,0.3,d_vals,K=500)

pdf(paste0(plot.path,'MM1_kpi_class1_example_C_bvals.pdf'))

plot(d_vals-1e-6, b_vals,
     main=expression(atop(paste("M/M/1 Optimal b Values for Class-1 Constraint,"), paste("where ", lambda[1], "=0.1, ", lambda[2], "=0.3, ", mu, "=1"))),
     xlab='d Values', ylab='b*(d)', ylim=c(0,1))

dev.off()

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.1, 0.3, b_vals[i], d_vals[i], K=1000)[,'Delayed.APQ.Class.2']
}

pdf(paste0(plot.path,'MM1_kpi_class1_example_C.pdf'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-2 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.1, ", lambda[2], "=0.3, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)', ylim=c(0.696,0.697))

dev.off()

#####################
# Plot Class-2 CDF for various (d, b*(d)) pairs
#####################

# lam1 = 0.3, lam2 = 0.1
d_max <- optim.db.class1(2,0.90,0.3,0.1,10,K=500)

#d_vals <- c(seq(1e-6,d_max,0.35))
d_vals <- c(0.2,0.4,0.6)
b_vals <- optim.b.class1(2,0.90,0.3,0.1,d_vals,K=500)

tvals <- seq(0.01,15,0.1)
class2_cdf <- cbind(tvals, 'd1'=rep(0, length(tvals)), 'd2'=rep(0, length(tvals)), 'd3'=rep(0, length(tvals)))
for (i in c(1,2,3)){
  d <- d_vals[i]
  b <- b_vals[i]
  
  cdf_func <- DP_CDF_class2(mu=1, lam1=0.3, lam2=0.1, b, d)
  class2_cdf[,i+1] <- cdf_func(tvals)
}

pdf(paste0(plot.path,'MM1_kpi_class1_example_A_cdf.pdf'))

plot(tvals, class2_cdf[,2],
     main=expression(atop(paste("M/M/1 Class-2 Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.3, ", lambda[2], "=0.1, ", mu, "=1"))),
     xlab='d Values', ylab='Class-2 CDF at d and b*(d)', ylim=c(0.6,1), type='l')
lines(tvals, class2_cdf[,3], lty=2)
lines(tvals, class2_cdf[,4], lty=3)


dev.off()

