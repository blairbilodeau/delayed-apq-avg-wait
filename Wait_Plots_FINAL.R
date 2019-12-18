###############################################################################
### This creates the CDF plots for the associated manuscript.
### Prepared by Blair Bilodeau on December 15, 2019
###############################################################################

code.path <- '/Users/blairbilodeau/Documents/Research/NSERC_USRA_2018/Avg_Wait_Paper/Code/'
source(paste0(code.path,"MM1_Wait_FINAL.R"))

plot.path <- '/Users/blairbilodeau/Documents/Research/NSERC_USRA_2018/Avg_Wait_Paper/Manuscript/2019_12_15/Plots/'

###############################################################################
### M/M/1
###############################################################################

point_spacing <- 4

# params
mu <- 1
lam1 <- 0.5
lam2 <- 0.3
b <- 0.8
K <- 1000
N <- 12

# functions
np.F <- NP_CDF_class2(mu, lam1, lam2, N)
dp.F <- DP_CDF_class2(mu, lam1, lam2, b, 0, K, N)
d6.F <- DP_CDF_class2(mu, lam1, lam2, b, 6, K, N)
d10.F <- DP_CDF_class2(mu, lam1, lam2, b, 10, K, N)

# data
x <- seq(0.01,50,0.5)
np <- np.F(x)
dp <- dp.F(x)
d6 <- d6.F(x)
d10 <- d10.F(x)

# d=6, d=10
png(paste0(path,'MM1_class2.png'))

plot(x, np, type='l', lty=1,  
     main=expression(paste("M/M/1 Class-2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1, ", b, "=0.8")), 
     xlab='Waiting Time (t)', ylab='P(W < t)', ylim=c(0,1))
points(x[seq(point_spacing, length(x), point_spacing)], np[seq(point_spacing, length(x), point_spacing)], pch=4, cex=0.8)

lines(x, d6, type='l', lty=2)
points(x[seq(point_spacing, length(x), point_spacing)], d6[seq(point_spacing, length(x), point_spacing)], pch=20)

lines(x, d10, type='l', lty=2)
points(x[seq(point_spacing, length(x), point_spacing)], d10[seq(point_spacing, length(x), point_spacing)], pch=18)

lines(x, dp, type='l', lty=1)

abline(v=6, lty=3)
abline(v=10, lty=3)

axis(1, at=c(6, seq(0,50,10)))
legend('bottomright', 
       c('APQ', 'Delayed APQ, d=6', 'Delayed APQ, d=10', 'Non-preemptive'), 
       lty = c(1,2,2,1), 
       pch=c(NA,20,18,4))

dev.off()