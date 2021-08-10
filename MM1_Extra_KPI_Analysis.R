###############################################################################
### This provides parametrization optimization for various KPI levels.
# This contains additional plots not included in the manuscript
###############################################################################

wdir <- '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/'

code.path <- paste0(wdir, 'code/')
source(paste0(code.path,"MM1_AvgWait.R"))
source(paste0(code.path,"MM1_Wait.R"))
source(paste0(code.path,"MM1_Class1.R"))
source(paste0(code.path,"MM1_KPI_Analysis.R"))

plot.path <- paste0(wdir, 'manuscript/plots/')


####
# Plot of NPQ boundary for 1-hour KPI
# to the left of the lines even NPQ can achieve KPI

pdf(paste0(plot.path,'MM1_npq_bound_4.pdf'))

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
       pch=c(20,15,4),
       cex=1.5)

dev.off()

####
# Plot of NPQ boundary for 1/2-hour KPI
# to the left of the lines even NPQ can achieve KPI

pdf(paste0(plot.path,'MM1_npq_bound_2.pdf'))

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
       pch=c(20,15,4),
       cex=1.5)

dev.off()

####
# Plot of FCFS boundary for 1-hour KPI
# to the right of the lines even FCFS can't achieve KPI

pdf(paste0(plot.path,'MM1_fcfs_bound_4.pdf'))

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
       pch=c(20,15,4),
       cex=1.5)

dev.off()

####
# Plot of FCFS boundary for 1/2-hour KPI
# to the right of the lines even FCFS can't achieve KPI

pdf(paste0(plot.path,'MM1_fcfs_bound_2.pdf'))

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
       pch=c(20,15,4),
       cex=1.5)

dev.off()

##
# combine class-1 and class-2
# note that the overlap is super tiny, so we chose not to include it
fcfs1_class2 <- fcfs.lam1(kpi_level=4, kpi_prob=0.85, lam1)
fcfs1_class1 <- fcfs.lam1(kpi_level=2, kpi_prob=0.9, lam1)
npq1_class2 <- npq.lam1(kpi_level=4, kpi_prob=0.85, lam1)
npq1_class1 <- npq.lam1.class1(kpi_level=2, kpi_prob=0.90, lam1)

plot(lam1, fcfs1_class1, type='l',
     main="M/M/1 FCFS 1 Hour Feasible Barrier (Class-1)", 
     xlab=expression(paste(lambda[1])), ylab=expression(paste(lambda[2])),
     xlim=c(0,1), ylim=c(0,1),
     lty=1)
lines(lam1, fcfs1_class2, type='l', lty=1, col='blue')
lines(lam1, npq1_class1, type='l', lty=3, col='black')
lines(lam1, npq1_class2, type='l', lty=3, col='blue')
segments(0,fcfs1_0.50[1],0,1)
legend('topright', 
       expression(paste('P(W'[1],' < 2',mu,')=0.50')),
       lty=1,
       cex=1.5)