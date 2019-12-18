###############################################################################
### This computes the waiting time for Class-2 customers in 
#   an M/M/1 Delayed APQ.
### Prepared by Blair Bilodeau on December 15, 2019
###############################################################################

#########
# LST Inversion Algorithms
#########

# gaver-stehfest (returns PDF)
GS <- function(L, N_2){
  f <- function(t){  
    N <- 2*N_2
    fn_sum <- 0
    for (i in 1:N){
      k <- floor((i+1)/2):min(i,N_2)
      Vi <- (-1)^(N_2 + i) * sum( (k^(N_2) * factorial(2*k)) / (factorial(N_2 - k) * factorial(k) * factorial(k-1) * factorial(i-k) * factorial(2*k - i)) )
      fn_sum <- fn_sum + Vi * L(i * log(2) / t)
    }
    log(2) * fn_sum / t 
  }
  Vectorize(f)
}

# coefficients from GaverStehfest
# use this to get the CDF rather than PDF
GS_coef <- matrix(0L, nrow=18, ncol=9)
GS_coef[1:2,1] <- c(2,-2)
GS_coef[1:4,2] <- c(-2,26,-48,24)
GS_coef[1:6,3] <- c(1,-49,366,-858,810,-270)
GS_coef[1:8,4] <- c(-3.333333e-01,
                    4.833333e+01,
                    -9.060000e+02,
                    5.464666e+03,
                    -1.437666e+04,
                    1.873000e+04,
                    -1.194666e+04,
                    2.986666e+03)
GS_coef[1:10,5] <- c(8.333333e-02,
                     -3.208333e+01,
                     1.279000e+03,
                     -1.562366e+04,
                     8.424416e+04,
                     -2.369575e+05,
                     3.759116e+05,
                     -3.400716e+05,
                     1.640625e+05,
                     -3.281250e+04)
GS_coef[1:12,6] <- c(-1.6666666e-02,
                     1.6016666e+01,
                     -1.2470000e+03,
                     2.7554333e+04,
                     -2.6328083e+05,
                     1.3241387e+06,
                     -3.8917055e+06,
                     7.0532863e+06,
                     -8.0053365e+06,
                     5.5528305e+06,
                     -2.1555072e+06,
                     3.5925120e+05)
GS_coef[1:14,7] <- c( 2.77777777777e-03,
                      -6.40277777777e+00,
                      9.24050000000e+02,
                      -3.45979277777e+04,
                      5.40321111111e+05,
                      -4.39834636666e+06,
                      2.10875917777e+07,
                      -6.39449130444e+07,
                      1.27597579550e+08,
                      -1.70137188083e+08,
                      1.50327467033e+08,
                      -8.45921615000e+07,
                      2.74788847666e+07,
                      -3.92555496666e+06)
GS_coef[1:16,8] <- c(-3.968253968253968E-04,
                     2.133730158730159E+00,
                     -5.510166666666667E+02,
                     3.350016111111111E+04,
                     -8.126651111111111E+05,
                     1.007618376666667E+07,
                     -7.324138297777778E+07,
                     3.390596320730159E+08,
                     -1.052539536278571E+09,
                     2.259013328583333E+09,
                     -3.399701984433333E+09,
                     3.582450461700000E+09,
                     -2.591494081366667E+09,
                     1.227049828766667E+09,
                     -3.427345554285714E+08,
                     4.284181942857143E+07)
GS_coef[1:18,9] <- c( 4.960317460317460E-05,
                      -6.095734126984128E-01,
                      2.745940476190476E+02,
                      -2.630695674603174E+04,
                      9.572572013888889E+05,
                      -1.735869484583333E+07,
                      1.824212226472222E+08,
                      -1.218533288309127E+09,
                      5.491680025283035E+09,
                      -1.736213111520684E+10,
                      3.945509690352738E+10,
                      -6.526651698517500E+10,
                      7.873006832822083E+10,
                      -6.855644419612083E+10,
                      4.198434347505357E+10,
                      -1.716093471183929E+10,
                      4.204550039102679E+09,
                      -4.671722265669643E+08)

# gaver-stehfest (returns CDF)
GS_CDF <- function(L, N){ # need even N
  coef_col <- N / 2
  i <- 1:N
  cdf <- function(x) sum(GS_coef[i, coef_col] / i * L(log(2) * i/x))
  Vectorize(cdf)
}


#########
## Function to compute the waiting time CDF LST
## K is truncation level of infinite sum
MM1_class2 <- function(mu, lam1, lam2, b, d, K){
  
  # constants
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
  
  # laplace transforms
  eta1A <- function(s) (s + mu + lam1A - sqrt((s+mu+lam1A)^2 - 4*mu*lam1A)) / (2*mu*lam1A)
  eta1 <- function(s) (s + mu + lam1 - sqrt((s+mu+lam1)^2 - 4*mu*lam1)) / (2*mu*lam1)
  B <- function(s) mu / (mu+s)
  
  # numerical params
  x <- matrix(0L, nrow=K, ncol=K)
  
  # setup first 2 rows
  x[1,1] <- r - p
  x[2,1] <- q*r*p
  x[2,2] <- r^2 - p^2
  
  # compute lower triangle x[k,ell]
  for (row in 3:K){
    x[row,1] <- q * x[row-1,2]
    for (col in 2:(K-1)){
      x[row,col] <- p * x[row-1,col-1] + q * x[row-1,col+1]
    }
    x[row,row] <- r^row - p^row
  }
  
  # delayed APQ greater than d
  LST_DP_d <- function(s){
    
    # multiply with J
    eta1As <- eta1A(s)
    J <- matrix(eta1As^(1:K),nrow=K,ncol=1)
    xK <- x %*% J
    if (abs(xK[K]) > 1e-10) print('ERROR: x_l did not converge')
    
    # compute poisson probabilities
    pK <- dpois(1:K, nu*d)
    
    # combine terms
    T1 <- sum(pK * xK)
    T2 <- rho * exp(-nu*d) * eta1As / (1-rho*eta1As) * exp(r*nu*d*eta1As)
    exp(-s*d) * (1-rho) * (T1 + T2)
  }
  
  # non-preemptive greater than d
  LST_NP_d <- function(s){
    
    # multiply with J
    eta1s <- eta1(s)
    J <- matrix(eta1s^(1:K),nrow=K,ncol=1)
    xK <- x %*% J
    if (abs(xK[K]) > 1e-10) print('ERROR: x_l did not converge')
    
    # compute poisson probabilities
    pK <- dpois(1:K, nu*d)
    
    # combine terms
    T1 <- sum(pK * xK)
    T2 <- rho * exp(-nu*d) * eta1s / (1-rho*eta1s) * exp(r*nu*d*eta1s)
    exp(-s*d) * (1-rho) * (T1 + T2)
  }
  
  # non-preemptive full
  LST_NP <- function(s){
    eta1s <- eta1(s)
    Bs <- B(s + lam1 - lam1*eta1s)
    (1-rho) * (s+lam1-lam1*eta1s) / (lam2*Bs - lam2 + s)
  }
  
  c(Vectorize(LST_DP_d), Vectorize(LST_NP_d), Vectorize(LST_NP))
}

## Use above computation of CDF LST to get CDF explicitly for NPQ
NP_CDF_class2 <- function(mu, lam1, lam2, N){
  f <- function(x){
    GS_CDF(MM1_class2(mu, lam1, lam2, b=0, d=0, K=500)[[3]], N)(x)
  }
  Vectorize(f)
}

## Use above computation of CDF LST to get CDF explicitly for Delayed APQ
DP_CDF_class2 <- function(mu, lam1, lam2, b, d, K, N){
  DP_d <- function(x){
    GS_CDF(MM1_class2(mu, lam1, lam2, b, d, K)[[1]], N)(x)
  }
  DP_d <- Vectorize(DP_d)
  
  NP_d <- function(x){
    GS_CDF(MM1_class2(mu, lam1, lam2, b, d, K)[[2]], N)(x)
  }
  NP_d <- Vectorize(NP_d)
  
  f <- function(x){
    NP_CDF_class2(mu, lam1, lam2, N)(x) - NP_d(x) + DP_d(x)
  }
  
  Vectorize(f)
}

