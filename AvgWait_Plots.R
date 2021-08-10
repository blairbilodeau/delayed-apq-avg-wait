###############################################################################
### This creates the average waiting time plots (Section 3.2) for the associated manuscript.
###############################################################################

wdir <- '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/'

code.path <- paste0(wdir, 'code/')
source(paste0(code.path,"MM1_AvgWait.R"))
source(paste0(code.path,"MD1_AvgWait.R"))

plot.path <- paste0(wdir, 'manuscript/plots/')

###############################################################################
### M/M/1
###############################################################################

#################
# b Values
#################

point_spacing <- 4
b_seq <- seq(0,1,0.01)
adj_b1 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=2, K=500)
adj_b2 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=4, K=500)
adj_b3 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=6, K=500)
adj_b4 <- Avg_MM1.b(mu=1, lam1=0.5, lam2=0.3, b=b_seq, d=8, K=500)

# b Class-1
pdf(paste0(plot.path,'MM1_bvals_class1.pdf'))

plot(b_seq, adj_b1['Non.Preemptive.Class.1',], type='l', lty=1,
     main=expression(paste("M/M/1 Class-1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Accumulation Rate (b)', ylab='Expected Waiting Time',
     ylim=c(0,5))

lines(b_seq, adj_b1['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=20)

lines(b_seq, adj_b2['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b2['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=18)

lines(b_seq, adj_b3['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b3['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=15, cex=0.8)

lines(b_seq, adj_b4['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b4['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=17, cex=0.8)

lines(b_seq, adj_b1['APQ.Class.1',], type='l', lty=1)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=4, cex=0.8)

legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=2', 'Delayed APQ, d=4', 'Delayed APQ, d=6', 'Delayed APQ, d=8', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,1), 
       pch=c(4,20,18,15,17,NA),
       cex=1.2)

dev.off()
# use dev.new() to see plot

# b Class-2
pdf(paste0(plot.path,'MM1_bvals_class2.pdf'))

plot(b_seq, adj_b1['Non.Preemptive.Class.2',], type='l', lty=1,
     main=expression(paste("M/M/1 Class-2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Accumulation Rate (b)', ylab='Expected Waiting Time',
     ylim=c(2,9))

lines(b_seq, adj_b1['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=20)

lines(b_seq, adj_b2['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b2['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=18)

lines(b_seq, adj_b3['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b3['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=15, cex=0.8)

lines(b_seq, adj_b4['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b4['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=17, cex=0.8)

lines(b_seq, adj_b1['APQ.Class.2',], type='l', lty=1)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=4, cex=0.8)

legend('bottomleft', 
       c( 'APQ', 'Delayed APQ, d=2', 'Delayed APQ, d=4', 'Delayed APQ, d=6', 'Delayed APQ, d=8', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,1), 
       pch=c(4,20,18,15,17,NA),
       cex=1.2)

dev.off()

#################
# d Values
#################

point_spacing <- 4
d_seq <- seq(0,15,0.3)
adj_d1 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.2, d=d_seq, K=1000)
adj_d2 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.4, d=d_seq, K=1000)
adj_d3 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.6, d=d_seq, K=1000)
adj_d4 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=0.8, d=d_seq, K=1000)
adj_d5 <- Avg_MM1.d(mu=1, lam1=0.5, lam2=0.3, b=1, d=d_seq, K=1000)

# d Class-1
pdf(paste0(plot.path,'MM1_dvals_class1.pdf'))

plot(d_seq, adj_d1['Non.Preemptive.Class.1',], type='l', lty=1,
     main=expression(paste("M/M/1 Class-1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Delay Length (d)', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')

lines(d_seq, adj_d1['Delayed.APQ.Class.1',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d1['Delayed.APQ.Class.1',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=20)
points(0, adj_d1['APQ.Class.1',1], pch=20)

lines(d_seq, adj_d2['Delayed.APQ.Class.1',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d2['Delayed.APQ.Class.1',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=18)
points(0, adj_d2['APQ.Class.1',1], pch=18)

lines(d_seq, adj_d3['Delayed.APQ.Class.1',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d3['Delayed.APQ.Class.1',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=15, cex=0.8)
points(0, adj_d3['APQ.Class.1',1], pch=15, cex=0.8)

lines(d_seq, adj_d4['Delayed.APQ.Class.1',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d4['Delayed.APQ.Class.1',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=17, cex=0.8)
points(0, adj_d4['APQ.Class.1',1], pch=17, cex=0.8)

lines(d_seq, adj_d5['Delayed.APQ.Class.1',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d5['Delayed.APQ.Class.1',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=25, cex=0.6, bg='black')
points(0, adj_d5['APQ.Class.1',1], pch=4, cex=0.8)

lines(d_seq, adj_d1['First.Come.First.Served',], type='l', lty=1)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d1['First.Come.First.Served',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=4, cex=0.8)

axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=0.2', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.8', 'Delayed APQ, b=1.0', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,2,1), 
       pch=c(4,20,18,15,17,25,NA),
       pt.bg=rep('black',7),
       cex=1.2,
       bg='white')

dev.off()

# d Class-2
pdf(paste0(plot.path,'MM1_dvals_class2.pdf'))

plot(d_seq, adj_d1['Non.Preemptive.Class.2',], type='l', lty=1,
     main=expression(paste("M/M/1 Class-2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Delay Length (d)', ylab='Expected Waiting Time',
     ylim=c(3,9), xaxt='n')

lines(d_seq, adj_d1['Delayed.APQ.Class.2',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d1['Delayed.APQ.Class.2',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=20)
points(0, adj_d1['APQ.Class.2',1], pch=20)

lines(d_seq, adj_d2['Delayed.APQ.Class.2',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d2['Delayed.APQ.Class.2',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=18)
points(0, adj_d2['APQ.Class.2',1], pch=18)

lines(d_seq, adj_d3['Delayed.APQ.Class.2',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d3['Delayed.APQ.Class.2',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=15, cex=0.8)
points(0, adj_d3['APQ.Class.2',1], pch=15, cex=0.8)

lines(d_seq, adj_d4['Delayed.APQ.Class.2',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d4['Delayed.APQ.Class.2',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=17, cex=0.8)
points(0, adj_d4['APQ.Class.2',1], pch=17, cex=0.8)

lines(d_seq, adj_d5['Delayed.APQ.Class.2',], type='l', lty=2)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d5['Delayed.APQ.Class.2',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=25, cex=0.6, bg='black')
points(0, adj_d5['APQ.Class.2',1], pch=4, cex=0.8)

lines(d_seq, adj_d1['First.Come.First.Served',], type='l', lty=1)
lines(d_seq[seq(point_spacing, length(d_seq), point_spacing)], adj_d1['First.Come.First.Served',][seq(point_spacing, length(d_seq), point_spacing)], type='p', pch=4, cex=0.8)

axis(1, at=seq(0,15,1))
legend('bottomright', 
       c( 'First Come First Served', 'Delayed APQ, b=0.2', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.8', 'Delayed APQ, b=1.0', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,2,1), 
       pch=c(4,20,18,15,17,25,NA),
       pt.bg=rep('black',7),
       cex=1.2,
       bg='white')

dev.off()

###############################################################################
### M/D/1
###############################################################################

#################
# b Values
#################

point_spacing <- 4
b_seq <- seq(0,1,0.02)
adj_b0 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=0)
adj_b1 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=2)
adj_b2 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=4)
adj_b3 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=6)
adj_b4 <- Avg_MD1.b(lam1=0.5, lam2=0.3, b=b_seq, d=8)

# b Class-1
pdf(paste0(plot.path,'MD1_bvals_class1.pdf'))

plot(b_seq, adj_b1['Non.Preemptive.Class.1',], type='l', lty=1,
     main=expression(paste("M/D/1 Class-1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Accumulation Rate (b)', ylab='Expected Waiting Time',
     ylim=c(0,5))

lines(b_seq, adj_b1['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=20)

lines(b_seq, adj_b2['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b2['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=18)

lines(b_seq, adj_b3['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b3['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=15, cex=0.8)

lines(b_seq, adj_b4['Delayed.APQ.Class.1',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b4['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=17, cex=0.8)

lines(b_seq, adj_b0['Delayed.APQ.Class.1',], type='l', lty=1)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b0['Delayed.APQ.Class.1',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=4, cex=0.8)

legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=2', 'Delayed APQ, d=4', 'Delayed APQ, d=6', 'Delayed APQ, d=8', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,1), 
       pch=c(4,20,18,15,17,NA),
       cex=1.2)

dev.off()
# use dev.new() to see plot

# b Class-2
pdf(paste0(plot.path,'MD1_bvals_class2.pdf'))

plot(b_seq, adj_b1['Non.Preemptive.Class.2',], type='l', lty=1,
     main=expression(paste("M/D/1 Class-2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Accumulation Rate (b)', ylab='Expected Waiting Time',
     ylim=c(2,9))

lines(b_seq, adj_b1['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b1['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=20)

lines(b_seq, adj_b2['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b2['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=18)

lines(b_seq, adj_b3['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b3['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=15, cex=0.8)

lines(b_seq, adj_b4['Delayed.APQ.Class.2',], type='l', lty=2)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b4['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=17, cex=0.8)

lines(b_seq, adj_b0['Delayed.APQ.Class.2',], type='l', lty=1)
lines(b_seq[seq(point_spacing, length(b_seq), point_spacing)], adj_b0['Delayed.APQ.Class.2',][seq(point_spacing, length(b_seq), point_spacing)], type='p', pch=4, cex=0.8)

legend('topleft', 
       c( 'APQ', 'Delayed APQ, d=2', 'Delayed APQ, d=4', 'Delayed APQ, d=6', 'Delayed APQ, d=8', 'Non-Preemptive'), 
       lty=c(1,2,2,2,2,1), 
       pch=c(4,20,18,15,17,NA),
       cex=1.2)

dev.off()

#################
# d Values
#################

point_spacing <- 4
d_seq <- seq(0,15,1)
adj_d1 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.2, d=d_seq)
adj_d2 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.4, d=d_seq)
adj_d3 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.6, d=d_seq)
adj_d4 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=0.8, d=d_seq)
adj_d5 <- Avg_MD1.d(lam1=0.5, lam2=0.3, b=1, d=d_seq)

# d Class-1
pdf(paste0(plot.path,'MD1_dvals_class1.pdf'))

plot(d_seq, adj_d1['Non.Preemptive.Class.1',], type='l', lty=1,
     main=expression(paste("M/D/1 Class-1, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Delay Length (d)', ylab='Expected Waiting Time',
     ylim=c(0,5), xaxt='n')

points(d_seq, adj_d1['Delayed.APQ.Class.1',], pch=20)

points(d_seq, adj_d2['Delayed.APQ.Class.1',], pch=18) 

points(d_seq, adj_d3['Delayed.APQ.Class.1',], pch=15, cex=0.8)

points(d_seq, adj_d4['Delayed.APQ.Class.1',], pch=17, cex=0.8)

points(d_seq, adj_d5['Delayed.APQ.Class.1',], pch=25, cex=0.6, bg='black') 

lines(d_seq, adj_d1['First.Come.First.Served',], type='l', lty=1)
points(d_seq, adj_d1['First.Come.First.Served',], pch=4, cex=0.8)

axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=0.2', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.8', 'Delayed APQ, b=1.0', 'Non-Preemptive'), 
       lty = c(1,NA,NA,NA,NA,NA,1),
       pch = c(4,20,18,15,17,25,NA), 
       pt.bg=rep('black',7),
       cex=1.2,
       bg='white')

dev.off()

# d Class-2
pdf(paste0(plot.path,'MD1_dvals_class2.pdf'))

plot(d_seq, adj_d1['Non.Preemptive.Class.2',], type='l', lty=1,
     main=expression(paste("M/D/1 Class-2, ", lambda[1], "=0.5, ", lambda[2], "=0.3, ", mu, "=1")), 
     xlab='Delay Length (d)', ylab='Expected Waiting Time',
     ylim=c(2,9), xaxt='n')

points(d_seq, adj_d1['Delayed.APQ.Class.2',], pch=20)

points(d_seq, adj_d2['Delayed.APQ.Class.2',], pch=18) 

points(d_seq, adj_d3['Delayed.APQ.Class.2',], pch=15, cex=0.8)

points(d_seq, adj_d4['Delayed.APQ.Class.2',], pch=17, cex=0.8)

points(d_seq, adj_d5['Delayed.APQ.Class.2',], pch=25, cex=0.6, bg='black') 

lines(d_seq, adj_d1['First.Come.First.Served',], type='l', lty=1)
points(d_seq, adj_d1['First.Come.First.Served',], pch=4, cex=0.8)

axis(1, at=seq(0,15,1))
legend('topright', 
       c( 'First Come First Served', 'Delayed APQ, b=0.2', 'Delayed APQ, b=0.4', 'Delayed APQ, b=0.6', 'Delayed APQ, b=0.8', 'Delayed APQ, b=1.0', 'Non-Preemptive'), 
       lty = c(1,NA,NA,NA,NA,NA,1),
       pch = c(4,20,18,15,17,25,NA), 
       pt.bg=rep('black',7),
       cex=1.2,
       bg='white')

dev.off()