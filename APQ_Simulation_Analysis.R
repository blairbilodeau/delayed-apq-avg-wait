
#########################
# params
#########################

wdir <- '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/'

code.path <- paste0(wdir, 'code/')
plot.path <- paste0(wdir, 'manuscript/plots/')
source(paste0(code.path,"MM1_Wait.R"))

path <- paste0(code.path, 'APQ_Simulation_Results/')

num_customers <- 6000
burn_in <- 1500 
burn_out <- 500
num_sims <- 50

#########################
# extraction functions
#########################

get.file <- function(num_customers, lam1, lam2, b, d, sim){
  filename <- paste0('APQsim', '_ncust_', toString(num_customers), '__', 'lam1_', toString(lam1), '_lam2_', toString(lam2), '__', 'b_', toString(b), '_d_', toString(d), '__', 'r_', toString(sim), '.csv')
  cust.df <- read.csv(paste0(path, filename))
  cust.df[burn_in:(nrow(cust.df)-burn_out),]
}

get.wait1cdf <- function(cust.df){
  wait <- cust.df[cust.df$level==1,'wait']
  m = length(wait)
  prob <- function(x){
    sum(wait <= x) / m
  }
  Vectorize(prob)
}

get.wait2cdf <- function(cust.df){
  wait <- cust.df[cust.df$level==2,'wait']
  m = length(wait)
  prob <- function(x){
    sum(wait <= x) / m
  }
  Vectorize(prob)
}

get.meancdf <- function(cdflist){
  prob <- function(x){
    mean(sapply(cdflist, mapply, x))
  }
  Vectorize(prob)
}

get.quantilecdf <- function(cdflist, quant){
  prob <- function(x){
    quantile(sapply(cdflist, mapply, x), quant)
  }
  Vectorize(prob)
}

#########################
# inital sanity checks
#########################

# mu <- 1
# lam1 <- 0.2
# lam2 <- 0.7
# 
# fcfs.waitcdf <- FCFS_CDF(mu,lam1+lam2)
# fcfs.waitcdf <- DP_CDF_class2(mu,lam1,lam2,0.6,0)
# fcfs.waitcdf <- NP_CDF_class2(mu,lam1,lam2)
# 
# fcfs.wait1cdfs <- c()
# 
# for (sim in seq(0,num_sims-1)){
#   fcfs.cust.df <- get.file(num_customers, lam1, lam2, b=0, d=0, sim)
#   fcfs.wait1cdfs <- c(fcfs.wait1cdfs, get.wait2cdf(fcfs.cust.df))
# }
# 
# fcfs.meanwait1cdf <- get.meancdf(fcfs.wait1cdfs)
# fcfs.minwait1cdf <- get.quantilecdf(fcfs.wait1cdfs,0.05)
# fcfs.maxwait1cdf <- get.quantilecdf(fcfs.wait1cdfs,0.95)
# 
# x <- seq(0.1,20,0.1)
# plot(x, fcfs.waitcdf(x), type='l', ylim=c(0,1),
#      main=bquote(atop("M/M/1 APQ Class-1 Exact vs. Simulated CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu))),
#      xlab='Waiting Time (t)', ylab='P(W < t)')
# lines(x, fcfs.meanwait1cdf(x), lty=3)
# lines(x, fcfs.minwait1cdf(x), lty=3)
# lines(x, fcfs.maxwait1cdf(x), lty=3)
# 
# 
# df1 <- fcfs.cust.df[fcfs.cust.df$level==1,]
# 
# cummean <- cumsum(df1$wait) / seq_along(1:nrow(df1))
# plot(1:nrow(df1), cummean)
# 
# df2 <- fcfs.cust.df[fcfs.cust.df$level==2,]
# 
# cummean <- cumsum(df2$wait) / seq_along(1:nrow(df2))
# plot(1:nrow(df2), cummean)

#########################
# comparisons as sanity check
#########################

mu <- 1
lam.pairs <- list(c(0.5,0.3), c(0.2,0.7))

for (lam.pair in lam.pairs){
  lam1 <- lam.pair[1]
  lam2 <- lam.pair[2]
  
  lam <- lam1+lam2
  
  fcfs.waitcdf <- FCFS_CDF(mu,lam)
  
  # simulated fcfs cdfs
  fcfs.wait1cdfs <- c()
  fcfs.wait2cdfs <- c()
  
  for (sim in seq(0,num_sims-1)){
    fcfs.cust.df <- get.file(num_customers, lam1, lam2, b=1, d=0, sim)
    fcfs.wait1cdfs <- c(fcfs.wait1cdfs, get.wait1cdf(fcfs.cust.df))  
    fcfs.wait2cdfs <- c(fcfs.wait2cdfs, get.wait2cdf(fcfs.cust.df))  
  }
  
  fcfs.meanwait1cdf <- get.meancdf(fcfs.wait1cdfs)
  fcfs.meanwait2cdf <- get.meancdf(fcfs.wait2cdfs)
  
  # simulated npq cdfs
  npq.wait1cdfs <- c()
  npq.wait2cdfs <- c()
  
  for (sim in seq(0,num_sims-1)){
    npq.cust.df <- get.file(num_customers, lam1, lam2, b=0, d=0, sim)
    npq.wait1cdfs <- c(npq.wait1cdfs, get.wait1cdf(npq.cust.df))  
    npq.wait2cdfs <- c(npq.wait2cdfs, get.wait2cdf(npq.cust.df))  
  }
  
  npq.meanwait1cdf <- get.meancdf(npq.wait1cdfs)
  npq.meanwait2cdf <- get.meancdf(npq.wait2cdfs)
  
  npq.wait1cdf <- NP_CDF_class1(mu,lam1,lam2)
  npq.wait2cdf <- NP_CDF_class2(mu,lam1,lam2)
  
  # class 1 plot
  
  plot.name <- paste0('simulated_class1_lam1_',toString(lam1),'_lam2_',toString(lam2),'.pdf')
  pdf(paste0(plot.path,plot.name))
  
  x <- seq(0,20,0.1)
  plot(x, fcfs.waitcdf(x), type='l', ylim=c(0,1),
       main=bquote(atop("M/M/1 APQ Class-1 Exact vs. Simulated CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu))), 
       xlab='Waiting Time (t)', ylab='P(W < t)')
  lines(x, fcfs.meanwait1cdf(x), lty=3)
  
  lines(x, npq.wait1cdf(x))
  lines(x, npq.meanwait1cdf(x), lty=3)
  
  # simulated apq cdfs
  bvals <- c(0.2, 0.4, 0.6, 0.8)
  
  apq.meanwait1cdf.list <- c()
  apq.wait1cdf.list <- c()
  
  for (b in bvals){
    apq.wait1cdfs <- c()
    apq.wait1cdf <- AP_CDF_class1(lam1, lam2, b)
    
    for (sim in seq(0,num_sims-1)){
      apq.cust.df <- get.file(num_customers, lam1, lam2, b, 0, sim)
      apq.wait1cdfs <- c(apq.wait1cdfs, get.wait1cdf(apq.cust.df))  
    }
    
    apq.meanwait1cdf <- get.meancdf(apq.wait1cdfs)
    
    lines(x, apq.wait1cdf(x))
    lines(x, apq.meanwait1cdf(x), lty=3)
  }
  
  text(15,0.2,labels=paste0("b = (0, ", toString(bvals), ", 1)"), cex=0.75)
  
  dev.off()
  
  # class 2 plot
  
  plot.name <- paste0('simulated_class2_lam1_',toString(lam1),'_lam2_',toString(lam2),'.pdf')
  pdf(paste0(plot.path,plot.name))
  
  x <- seq(0.1,20,0.1)
  plot(x, fcfs.waitcdf(x), type='l', ylim=c(0,1),
       main=bquote(atop("M/M/1 APQ Class-2 Exact vs. Simulated CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu))), 
       xlab='Waiting Time (t)', ylab='P(W < t)')
  lines(x, fcfs.meanwait2cdf(x), lty=3)
  
  lines(x, npq.wait2cdf(x))
  lines(x, npq.meanwait2cdf(x), lty=3)
  
  # simulated apq cdfs
  bvals <- c(0.2, 0.4, 0.6, 0.8)
  
  apq.meanwait2cdf.list <- c()
  apq.wait2cdf.list <- c()
  
  for (b in bvals){
    apq.wait2cdfs <- c()
    apq.wait2cdf <- DP_CDF_class2(mu, lam1, lam2, b, 0, K=500)
    
    for (sim in seq(0,num_sims-1)){
      apq.cust.df <- get.file(num_customers, lam1, lam2, b, 0, sim)
      apq.wait2cdfs <- c(apq.wait2cdfs, get.wait2cdf(apq.cust.df))  
    }
    
    apq.meanwait2cdf <- get.meancdf(apq.wait2cdfs)
    
    lines(x, apq.wait2cdf(x))
    lines(x, apq.meanwait2cdf(x), lty=3)
  }
  
  text(15,0.2,labels=paste0("b = (0, ", toString(bvals), ", 1)"), cex=0.75)
  
  dev.off()
}

#########################
# d-APQ Class-1 simulation comparison
#########################

mu <- 1
lam.pairs <- list(c(0.5,0.3), c(0.2,0.7))
dvals <- c(2,6)

bvals <- c(0.2,0.4,0.6,0.8)
t <- c(0.001, seq(0.01, 30, 0.1))

for (lam.pair in lam.pairs){
  lam1 <- lam.pair[1]
  lam2 <- lam.pair[2]
  
  for (d in dvals){
    plot.name <- paste0('diff_simulated_d_',toString(d),'_lam1_',toString(lam1),'_lam2_',toString(lam2),'.pdf')
    pdf(paste0(plot.path,plot.name))
    
    # use FCFS
    plot(t, rep(0,length(t)), type='l',ylim=c(0,0.1),
         main=bquote(atop("M/M/1 Delayed APQ Class-1 Simulated vs. Approximate CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu) ~ ", " ~ d == .(d))),
         ylab = "Absolute Difference",
         xlab = "Time (t)")
    
    # inbetween FCFS and NPQ
    for (i in 1:length(bvals)){
      b <- bvals[i]
      
      dapq.wait1approxcdf <- DP_CDF_approx_class1(mu, lam1, lam2, b, d, K=500)
      dapq.wait1cdfs <- c()
      
      for (sim in seq(0,num_sims-1)){
        dapq.cust.df <- get.file(num_customers, lam1, lam2, b, d, sim)
        dapq.wait1cdfs <- c(dapq.wait1cdfs, get.wait1cdf(dapq.cust.df))  
      }
      
      dapq.meanwait1cdf <- get.meancdf(dapq.wait1cdfs)
      
      lines(t, abs(dapq.wait1approxcdf(t) - dapq.meanwait1cdf(t)), lty=i)
      
      #text(xtextvals[i],ytextvals[i],labels=bquote(~ b == .(round(b,1))), cex=0.75)
    }
    
    #! check bvals !#
    legend('topright', 
           c('b=0.2', 'b=0.4', 'b=0.6', 'b=0.8'), 
           lty=c(1,2,3,4))
    
    dev.off()
    
  }

}






