###############################################################################
### Plots comparisons for an exponential approximation of the waiting time for Class-1 customers in 
#   an M/M/1 Delayed APQ.
###############################################################################

wdir <- '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/'

code.path <- paste0(wdir, 'code/')
source(paste0(code.path,"MM1_Wait.R"))
source(paste0(code.path,"MM1_AvgWait.R"))

plot.path <- paste0(wdir, 'manuscript/plots/')


#########################
# Plotting
#########################

####
# Exact vs. approx
####

bvals <- c(0.2,0.4,0.6,0.8)
#xtextvals <- c(5, 10, 15, 20, 25, 30)
#ytextvals <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1)

mu <- 1
lam1 <- 0.7
lam2 <- 0.2
rho <- lam1 + lam2
rho2 <- lam2
lam <- lam1 + lam2

t <- c(0.001, seq(0.01, 30, 0.1))

plot.name <- paste0('approx_lam1_',toString(lam1),'_lam2_',toString(lam2),'.pdf')
pdf(paste0(plot.path,plot.name))

# use FCFS
fcfs.waitcdf <- function(x) 1 - rho*exp(-mu*(1-rho)*x)
W1tvec <- fcfs.waitcdf(t)
plot(t, W1tvec, type='l',ylim=c(0,1), lty=3,
     main=bquote(atop("M/M/1 APQ Class-1 Exact vs. Approximate CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu))),
     ylab = expression(paste("P(", W[1], "< t)")),
     xlab = "Time (t)")
text(25,0.2,labels=paste0("b = (0, ", toString(bvals), ", 1)"), cex=0.75)

# inbetween FCFS and NPQ
for (i in 1:length(bvals)){
  b <- bvals[i]

  W1tvec <- AP_CDF_class1(lam1, lam2, b, N=12)(t)
  lines(t, W1tvec)
  
  rho1A <- lam1*(1-b)
  EW_APQ2 <- lam / (mu^2 * (1-rho) * (1-rho1A))
  EW_APQ1 <- lam / mu^2 / (1-rho) - rho2 * EW_APQ2 * (1-b)
  W1t.approx <- function(t) {1 - (lam)*exp(-((lam)/EW_APQ1)*t)}
  lines(t, W1t.approx(t), lty=2)
  
  #text(xtextvals[i],ytextvals[i],labels=bquote(~ b == .(round(b,1))), cex=0.75)
}

# use NPQ
W1tvec <- NP_CDF_class1(mu,lam1,lam2)(t)
lines(t, W1tvec,lty=3)

dev.off()

####
# Absolute differences
####

bvals <- c(0.2,0.4,0.6,0.8)
#xtextvals <- c(5, 10, 15, 20, 25, 30)
#ytextvals <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1)

mu <- 1
lam1 <- 0.7
lam2 <- 0.2
rho <- lam1 + lam2
rho2 <- lam2
lam <- lam1 + lam2

t <- c(0.001, seq(0.01, 30, 0.1))

plot.name <- paste0('diff_lam1_',toString(lam1),'_lam2_',toString(lam2),'.pdf')
pdf(paste0(plot.path,plot.name))

# use FCFS
fcfs.waitcdf <- function(x) 1 - rho*exp(-mu*(1-rho)*x)
W1tvec <- fcfs.waitcdf(t)
plot(t, rep(0,length(t)), type='l',ylim=c(0,0.1),
     main=bquote(atop("M/M/1 APQ Class-1 Exact vs. Approximate CDF, ", "where " ~ lambda[1] == .(lam1) ~ ", " ~ lambda[2] == .(lam2) ~ ", " ~ mu == .(mu))),
     ylab = "Absolute Difference",
     xlab = "Time (t)")

# inbetween FCFS and NPQ
for (i in 1:length(bvals)){
  b <- bvals[i]
  
  W1tvec <- AP_CDF_class1(lam1, lam2, b, N=12)(t)
  
  rho1A <- lam1*(1-b)
  EW_APQ2 <- lam / (mu^2 * (1-rho) * (1-rho1A))
  EW_APQ1 <- lam / mu^2 / (1-rho) - rho2 * EW_APQ2 * (1-b)
  W1t.approx <- function(t) {1 - (lam)*exp(-((lam)/EW_APQ1)*t)}
  lines(t, abs(W1tvec - W1t.approx(t)), lty=i)
  
  #text(xtextvals[i],ytextvals[i],labels=bquote(~ b == .(round(b,1))), cex=0.75)
}

#! check bvals !#
legend('topright', 
       c('b=0.2', 'b=0.4', 'b=0.6', 'b=0.8'), 
       lty=c(1,2,3,4))

dev.off()
