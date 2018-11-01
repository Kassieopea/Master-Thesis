## estimate Gaussian DTSM using block-wise Metropolis-Hastings algorithm
## parameterization in terms of RISK PRICES
## theta = (kinfQ, lamQ, Sigma, lam0, lam1, sigw)
## lambda = vec(lam0, lam1)
setwd("C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/R_MEF")

rm(list=ls())

# load("Testdata_N3.RData")
source("R/rrp_functions_MEF.r")
source("R/jsz_fns_MEF.r")
source("R/estimation_fns_MEF.r")
set.seed(4493)
library(corpcor)
 
# ###load data, init parameters
init(N=3)
# browser()
## model specification
model <- "M3"
gamma <- numeric(N+(N*(2+N)))

# switch(model,
#        M0 = {gamma[] <- 1},
#        M1 = {gamma[c(4,13,17)] <- 1},
#        M2 = {gamma[c(10,17)] <- 1},
#        # M3 = {gamma[c(3,4,12,13,17)] <- 1})
#        M3 = {gamma[c(10,15,17)] <- 1})
switch(model,
       M0 = {gamma[] <- 1},
       M1 = {gamma[c(3,4,8,12,13,17)] <- 1},
       M2 = {gamma[c(3,4,8,12,13,14,17)] <- 1},
       # M3 = {gamma[c(3,4,12,13,17)] <- 1})
       M3 = {gamma[c(3,4,8,12,15,17)] <- 1})
cat("model to estimate: ", model, "\n", gamma, "\n")

#####get starting values
pars.mle <- estML(Y, W, mats, gamma)
pars <- pars.mle

## priors
# browser()
priors <- getPriors()

## diagonalized g-prior for parameters
# browser()
crossterm <- cbind(diag(N), matrix(rep(0, (2*N)), nrow =N))
#mu.Q <- rbind(pars$loads$K0Q.cP, matrix(rep(0,2), nrow = 2))
mu.Q <- pars$loads$K0Q.cP
phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
rvar.res <- getLambdaRVAR(rep(1, N+(N*(2+N))), cP=cbind(Y%*%t(W), MF), mu.Q, phi.Q, pars$OmegaInv, crossterm)
# rvar.res <- getLambdaRVAR(rep(1, N*(N+1)), cP=Y%*%t(W), pars$loads$K0Q.cP, pars$loads$K1Q.cP, OmegaInv)
priors$lambda.sd <- sqrt(priors$g)*sqrt(diag(rvar.res$cov.mat))

## show significance
# browser()
getCondPostLambda(rvar.res)

if (all(gamma==1)) {
  cat("Estimating unrestricted model...\n")
} else {
  cat("Estimating restricted model (", vec2str(gamma), ") ...\n")
}

results.file <- getResultsFileName(paste("THESIS_mcmc", model, sep="_"), N, gamma)

########################################
## Metropolis-Hastings algorithm

Scale.Omega.1 <<- 5000
Scale.Omega.2 <<-100000
Scale.Theta <<- 100
Scale.Phi <<- 25

# Scale.Omega.1 <<- 1
# Scale.Omega.2 <<-100000
# Scale.Theta <<- 25
# Scale.Phi <<- 50

## number of iterations
M <- 15000

nu <<- 5

#browser()
## draws
Omega.i <- matrix(NA, M, (N+2)^2)
lambda.i <- matrix(NA, M, N+(N*(2+N)))
lamQ.i <- matrix(NA, M, N)
kinfQ.i <- matrix(NA, M, 1)
sige2.i <- matrix(NA, M, 1)
Phi.i <- matrix(NA, M, (N+2)^2)
Mu.i <- matrix(NA, M, (N+2))

## acceptance probabilities
alpha.Omega <- numeric(M)
alpha.thetaQ <- numeric(M)
alpha.par <- numeric(M)

count.nonstat <- 0
# 
save.image("Testdata_N3.RData")
# browser()

cat("Running MCMC algorithm...\n")
tic()
# browser()
for (m in 1:M) {
  if( m %% 1000 == 0) {
    cat("*** iteration ",m,"\n")
    cat("acceptance probabilities:\n")
    cat("thetaQ => ", round(100*mean(alpha.thetaQ[1:(m-1)])), "\n")
    cat("Omega  => ", round(100*mean(alpha.Omega[1:(m-1)])), "\n")
    cat("PAR => ", round(100*mean(alpha.par[1:(m-1)])), "\n")
    cat("current amount of non-stat: ", count.nonstat, "\n")
    cat("lambda MEAN :", colMeans(lambda.i[1:(m-1),]), "\n\n")
    
  }
  pars.save <- pars
  # browser()
  # cat("HELLO 11111")
  pars <- drawLambda(pars)
  # browser()
  pars <- drawThetaQ(pars)
  # browser()
  pars <- drawPhi(pars)
  pars <- drawOmega(pars)
  # browser()
  pars <- drawSige2(pars)
  
  ## if non-stationary, retain previous iteration's values
  if (max(abs(eigen(pars$Phi[1:N,1:N], only.values=TRUE)$values))>=1) {
    pars <- pars.save
    count.nonstat <- count.nonstat + 1
    # cat("eigenvalue: " , max(abs(eigen(pars$Phi[1:N,1:N], only.values=TRUE)$values)), "\n")
    # cat("lambda par:", pars$lambda, "\n\n")
  }
  # browser()
  lambda.i[m,] <- pars$lambda
  Omega.i[m,] <- as.vector(pars$Omega)
  kinfQ.i[m] <- pars$kinfQ
  lamQ.i[m,] <- pars$lamQ
  sige2.i[m] <- pars$sige2
  Phi.i[m,] <- as.vector(pars$Phi)
  Mu.i[m,] <- pars$mu
}
toc()

alpha.kinfQ <- alpha.thetaQ
alpha.lamQ <- alpha.thetaQ

save(file=results.file, gamma, priors, M, N, mats, W, dates, Omega.i, lamQ.i, kinfQ.i, lambda.i, sige2.i, Phi.i, Mu.i, alpha.Omega, alpha.thetaQ, alpha.kinfQ, alpha.lamQ, alpha.par)
