###############################################################################
### This provides parametrization optimization for various KPI levels.
### Prepared by Blair Bilodeau on December 15, 2019
###############################################################################

code.path <- '/Users/blairbilodeau/Documents/Research/NSERC_USRA_2018/Avg_Wait_Paper/Code/'
source(paste0(code.path,"MM1_AvgWait_FINAL.R"))
source(paste0(code.path,"MM1_Wait_FINAL.R"))

plot.path <- '/Users/blairbilodeau/Documents/Research/NSERC_USRA_2018/Avg_Wait_Paper/Manuscript/2019_12_27/Plots/'

##################
## NPQ Boundary
##################

### for fixed lam1 find maximum lam2 such that KPI can be achieved with NPQ
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


####
# Plot of NPQ boundary for 1-hour KPI
# to the left of the lines even NPQ can achieve KPI

png(paste0(plot.path,'MM1_npq_bound_4.png'))

tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

npq1 <- npq.lam1(kpi_level=4, kpi_prob=0.95, lam1)
npq2 <- npq.lam2(kpi_level=4, kpi_prob=0.95, lam2)
lam1.tot <- c(lam1, npq2)
lam2.tot <- c(npq1, lam2)

point_spacing <- 8
plot(lam1.tot, lam2.tot, type='l',
     main="M/M/1 NPQ 1 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1[seq(0, length(lam1), point_spacing)], npq1[seq(0, length(lam1), point_spacing)], type='p', pch=4)
lines(npq2[length(lam2)], lam2[length(lam2)], type='p', pch=4)

kpi_probs <- c(0.5, 0.8)
plot_symb <- c(20,15)
for (i in 1:length(kpi_probs)){
  npq1 <- npq.lam1(kpi_level=4, kpi_probs[i], lam1)
  npq2 <- npq.lam2(kpi_level=4, kpi_probs[i], lam2)
  lam1.tot <- c(lam1, npq2)
  lam2.tot <- c(npq1, lam2)
  
  point_spacing <- 8
  lines(lam1.tot, lam2.tot, type='l')
  lines(lam1[seq(0, length(lam1), point_spacing)], npq1[seq(0, length(lam1), point_spacing)], type='p', pch=plot_symb[i])
}

legend('topright', 
       c('P(W<4)=0.5', 'P(W<4)=0.8', 'P(W<4)=0.95'), 
       lty=rep(1,3), 
       pch=c(20,15,4))

dev.off()

####
# Plot of NPQ boundary for 1/2-hour KPI
# to the left of the lines even NPQ can achieve KPI

png(paste0(plot.path,'MM1_npq_bound_2.png'))

tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

npq1 <- npq.lam1(kpi_level=2, kpi_prob=0.95, lam1)
npq2 <- npq.lam2(kpi_level=2, kpi_prob=0.95, lam2)
lam1.tot <- c(lam1, npq2)
lam2.tot <- c(npq1, lam2)

point_spacing <- 8
plot(lam1.tot, lam2.tot, type='l',
     main="M/M/1 NPQ 1/2 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1[seq(0, length(lam1), point_spacing)], npq1[seq(0, length(lam1), point_spacing)], type='p', pch=4)
lines(npq2[length(lam2)], lam2[length(lam2)], type='p', pch=4)

kpi_probs <- c(0.5, 0.8)
plot_symb <- c(20,15)
for (i in 1:length(kpi_probs)){
  npq1 <- npq.lam1(kpi_level=2, kpi_probs[i], lam1)
  npq2 <- npq.lam2(kpi_level=2, kpi_probs[i], lam2)
  lam1.tot <- c(lam1, npq2)
  lam2.tot <- c(npq1, lam2)
  
  point_spacing <- 8
  lines(lam1.tot, lam2.tot, type='l')
  lines(lam1[seq(0, length(lam1), point_spacing)], npq1[seq(0, length(lam1), point_spacing)], type='p', pch=plot_symb[i])
}

legend('topright', 
       c('P(W<2)=0.5', 'P(W<2)=0.8', 'P(W<2)=0.95'), 
       lty=rep(1,3), 
       pch=c(20,15,4))

dev.off()


##################
## FCFS Boundary
##################
# https://www.win.tue.nl/~iadan/que/h4.pdf
# Eq(14), page 6

### for fixed lam1 find minimum lam2 such that KPI can't be achieved with FCFS
fcfs.lam1 <- function(kpi_level, kpi_prob, lam1){
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
fcfs.lam2 <- function(kpi_level, kpi_prob, lam2){
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

####
# Plot of FCFS boundary for 1-hour KPI
# to the right of the lines even FCFS can't achieve KPI

png(paste0(plot.path,'MM1_fcfs_bound_4.png'))

tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

fcfs1 <- fcfs.lam1(kpi_level=4, kpi_prob=0.95, lam1)
fcfs2 <- fcfs.lam2(kpi_level=4, kpi_prob=0.95, lam2)
lam1.tot <- c(lam1, fcfs2)
lam2.tot <- c(fcfs1, lam2)

point_spacing <- 8
plot(lam1.tot, lam2.tot, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1[seq(0, length(lam1), point_spacing)], fcfs1[seq(0, length(lam1), point_spacing)], type='p', pch=4)
lines(fcfs2[length(lam2)], lam2[length(lam2)], type='p', pch=4)

kpi_probs <- c(0.5, 0.8)
plot_symb <- c(20,15)
for (i in 1:length(kpi_probs)){
  fcfs1 <- fcfs.lam1(kpi_level=4, kpi_probs[i], lam1)
  fcfs2 <- fcfs.lam2(kpi_level=4, kpi_probs[i], lam2)
  lam1.tot <- c(lam1, fcfs2)
  lam2.tot <- c(fcfs1, lam2)
  
  point_spacing <- 8
  lines(lam1.tot, lam2.tot, type='l')
  lines(lam1[seq(0, length(lam1), point_spacing)], fcfs1[seq(0, length(lam1), point_spacing)], type='p', pch=plot_symb[i])
}

legend('topright', 
       c('P(W<4)=0.5', 'P(W<4)=0.8', 'P(W<4)=0.95'), 
       lty=rep(1,3), 
       pch=c(20,15,4))

dev.off()

####
# Plot of FCFS boundary for 1/2-hour KPI
# to the right of the lines even FCFS can't achieve KPI

png(paste0(plot.path,'MM1_fcfs_bound_2.png'))

tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

fcfs1 <- fcfs.lam1(kpi_level=2, kpi_prob=0.95, lam1)
fcfs2 <- fcfs.lam2(kpi_level=2, kpi_prob=0.95, lam2)
lam1.tot <- c(lam1, fcfs2)
lam2.tot <- c(fcfs1, lam2)

point_spacing <- 8
plot(lam1.tot, lam2.tot, type='l',
     main="M/M/1 FCFS 1/2 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1[seq(0, length(lam1), point_spacing)], fcfs1[seq(0, length(lam1), point_spacing)], type='p', pch=4)
lines(fcfs2[length(lam2)], lam2[length(lam2)], type='p', pch=4)

kpi_probs <- c(0.5, 0.8)
plot_symb <- c(20,15)
for (i in 1:length(kpi_probs)){
  fcfs1 <- fcfs.lam1(kpi_level=2, kpi_probs[i], lam1)
  fcfs2 <- fcfs.lam2(kpi_level=2, kpi_probs[i], lam2)
  lam1.tot <- c(lam1, fcfs2)
  lam2.tot <- c(fcfs1, lam2)
  
  point_spacing <- 8
  lines(lam1.tot, lam2.tot, type='l')
  lines(lam1[seq(0, length(lam1), point_spacing)], fcfs1[seq(0, length(lam1), point_spacing)], type='p', pch=plot_symb[i])
}

legend('topright', 
       c('P(W<2)=0.5', 'P(W<2)=0.8', 'P(W<2)=0.95'), 
       lty=rep(1,3), 
       pch=c(20,15,4))

dev.off()

##################
## feasible region plotted together
##################

tol <- 1e-6
lam1 <- seq(tol,1-tol,0.01)
lam2 <- seq(tol,1-tol,0.01)

##
# 1 hour

####
# Plot of NPQ lower and FCFS upper boundary

fcfs1_0.95 <- fcfs.lam1(kpi_level=4, kpi_prob=0.95, lam1)
fcfs1_0.8 <- fcfs.lam1(kpi_level=4, kpi_prob=0.8, lam1)
fcfs1_0.5 <- fcfs.lam1(kpi_level=4, kpi_prob=0.5, lam1)

npq1_0.95 <- npq.lam1(kpi_level=4, kpi_prob=0.95, lam1)
npq1_0.8 <- npq.lam1(kpi_level=4, kpi_prob=0.8, lam1)
npq1_0.5 <- npq.lam1(kpi_level=4, kpi_prob=0.5, lam1)

png(paste0(plot.path,'MM1_feasible_bound_4.png'))

plot(lam1, fcfs1_0.95, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.95, type='l')
segments(0,fcfs1_0.95[1],0,1)

lines(lam1, npq1_0.8, type='l', lty=2)
lines(lam1, fcfs1_0.8, type='l', lty=2)

lines(lam1, npq1_0.5, type='l', lty=3)
lines(lam1, fcfs1_0.5, type='l', lty=3)

legend('topright', 
       c('P(W<4)=0.5', 'P(W<4)=0.8', 'P(W<4)=0.95'), 
       lty=c(3,2,1))

dev.off()

##
# 1/2 hour

fcfs1_0.95 <- fcfs.lam1(kpi_level=2, kpi_prob=0.95, lam1)
fcfs1_0.8 <- fcfs.lam1(kpi_level=2, kpi_prob=0.8, lam1)
fcfs1_0.5 <- fcfs.lam1(kpi_level=2, kpi_prob=0.5, lam1)

npq1_0.95 <- npq.lam1(kpi_level=2, kpi_prob=0.95, lam1)
npq1_0.8 <- npq.lam1(kpi_level=2, kpi_prob=0.8, lam1)
npq1_0.5 <- npq.lam1(kpi_level=2, kpi_prob=0.5, lam1)

png(paste0(plot.path,'MM1_feasible_bound_2.png'))

plot(lam1, fcfs1_0.95, type='l',
     main="M/M/1 FCFS 1/2 Hour Feasible Barrier", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1))
lines(lam1, npq1_0.95, type='l')
segments(0,fcfs1_0.95[1],0,1)

lines(lam1, npq1_0.8, type='l', lty=2)
lines(lam1, fcfs1_0.8, type='l', lty=2)

lines(lam1, npq1_0.5, type='l', lty=3)
lines(lam1, fcfs1_0.5, type='l', lty=3)

legend('topright', 
       c('P(W<2)=0.5', 'P(W<2)=0.8', 'P(W<2)=0.95'), 
       lty=c(3,2,1))

dev.off()

####################
## optimizing b and d in the case where FCFS possible and NPQ too penalizing
####################

### b and F2(x) increase together
### b and class-1 expectation increase together
### d goes down means class-1 expectation increases
### find b where F(x) = p (smallest b that KPI is met)
optim.b <- function(kpi_level, kpi_prob, lam1, lam2, d, K=500, N=12, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(b) {DP_CDF_class2(mu, lam1, lam2, b, d, K, N)(kpi_level) - kpi_prob}
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol/10,1-tol/10), tol=tol)
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

### find d where F(x) = p (smallest d that KPI is met)
optim.d <- function(kpi_level, kpi_prob, lam1, lam2, b, max_d=25, K=500, N=12, tol=1e-6){
  mu <- 1
  
  cdf_diff <- function(d) {DP_CDF_class2(mu, lam1, lam2, b, d, K, N)(kpi_level) - kpi_prob}
  
  err_try <- function(){
    ans <- uniroot(cdf_diff, c(tol/10,max_d-tol/10), tol=tol)
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

### find d where there exists b such that F(x) = p (largest d such that KPI can be met)
optim.db <- function(kpi_level, kpi_prob, lam1, lam2, max_d, K=500, N=12, tol=1e-4){
  
  optim_diff <- function(d) {optim.b(kpi_level, kpi_prob, lam1, lam2, d, K, N, tol) - (1-tol/10)}
  
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

npq1_0.8 <- npq.lam1(kpi_level=4, kpi_prob=0.8, lam1)
fcfs1_0.8 <- fcfs.lam1(kpi_level=4, kpi_prob=0.8, lam1)

#####################
# Plot pairs of d and Class-1 expected waiting time at d*(b)
#####################

##### P(W<4)>0.8 example A
match(0.65+tol,lam1)
npq1_0.8[match(0.65+tol,lam1)]
fcfs1_0.8[match(0.65+tol,lam1)]

# lam1 = 0.65, lam2 = 0.025
d_max <- optim.db(4,0.8,0.65,0.025,1,K=1000,N=12)

d_vals <- c(seq(1e-6,d_max,0.025))
b_vals <- optim.b(4,0.8,0.65,0.025,d_vals,K=1000)

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.65, 0.025, b_vals[i], d_vals[i], K=1000)[,'Delayed.APQ.Class.1']
}

png(paste0(plot.path,'MM1_kpi_example_A.png'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.65, ", lambda[2], "=0.025, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)')

dev.off()

##### P(W<4)>0.8 example B
match(0.32+tol,lam1)
npq1_0.8[match(0.32+tol,lam1)]
fcfs1_0.8[match(0.32+tol,lam1)]

# lam1 = 0.32, lam2 = 0.32
d_max <- optim.db(4,0.8,0.32,0.32,5)

d_vals <- c(seq(1e-6,d_max,0.1))
b_vals <- optim.b(4,0.8,0.32,0.32,d_vals)

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.32, 0.32, b_vals[i], d_vals[i], K=500)[,'Delayed.APQ.Class.1']
}

png(paste0(plot.path,'MM1_kpi_example_B.png'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.32, ", lambda[2], "=0.32, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)')

dev.off()

##### P(W<4)>0.8 example C
match(0.2+tol,lam1)
npq1_0.8[match(0.2+tol,lam1)]
fcfs1_0.8[match(0.2+tol,lam1)]

# lam1 = 0.2, lam2 = 0.45
d_max <- optim.db(4,0.8,0.2,0.45,5)

d_vals <- c(seq(1e-6,d_max,0.1))
b_vals <- optim.b(4,0.8,0.2,0.45,d_vals)

avg_mm1 <- numeric(length(d_vals))
for (i in 1:length(d_vals)){
  avg_mm1[i] <- Avg_MM1(1, 0.2, 0.45, b_vals[i], d_vals[i], K=500)[,'Delayed.APQ.Class.1']
}

png(paste0(plot.path,'MM1_kpi_example_C.png'))

plot(d_vals-1e-6, avg_mm1,
     main=expression(atop(paste("M/M/1 Class-1 Expected Waiting Time for Optimal d and b Values,"), paste("where ", lambda[1], "=0.2, ", lambda[2], "=0.45, ", mu, "=1"))),
     xlab='d Values', ylab='Expected Waiting Time at d and b*(d)')

dev.off()


