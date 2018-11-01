## estimate Gaussian DTSM using block-wise Metropolis-Hastings algorithm
## RJMCMC - Reversible-Jump Markov Chain Monte Carlo

rm(list=ls())
# load("Testdata_RJMCMC.RData")
source("R/rrp_functions_MEF.r")
source("R/jsz_fns_MEF.r")
source("R/estimation_fns_MEF.r")
library(corpcor)
set.seed(4493)
init(N=3)

## starting values
pars.mle <- estML(Y, W, mats)
pars <- pars.mle

## priors
priors <- getPriors()

## diagonalized g-prior for parameters
crossterm <<- cbind(diag(N), matrix(rep(0, (N*2)), nrow =N))
phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (N*2)), nrow = N))
rvar.res <- getLambdaRVAR(rep(1, (N+(N+2)*N)), cP=cbind(Y%*%t(W),MF), pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
priors$lambda.sd <- sqrt(priors$g)*sqrt(diag(rvar.res$cov.mat))

results.file <- getResultsFileName("THESIS_rjmcmc", N)

M <- 55000 ## number of iterations

## draws
Omega.i <- matrix(NA, M, (N+2)^2)
lambda.i <- matrix(NA, M, (N+(N+2)*N))
gamma.i <- matrix(NA, M, (N+(N+2)*N))
lamQ.i <- matrix(NA, M, N)
kinfQ.i <- matrix(NA, M, 1)
sige2.i <- matrix(NA, M, 1)
Phi.i <- matrix(NA, M, (N+2)^2)
Mu.i <- matrix(NA, M, (N+2))

## acceptance probabilities
alpha.Omega <- numeric(M)
alpha.thetaQ <- numeric(M)
alpha.jump <- numeric(M)*NA
alpha.par <- numeric(M)
count.nonstat <- 0

save.image("Testdata_RJMCMC.RData")
Scale.Omega.1 <<- 5000
Scale.Omega.2 <<-100000
Scale.Theta <<- 100
Scale.Phi <<- 25


cat("Running MCMC algorithm...\n")
for (m in 1:M) {
  if( m %% 1000 == 0) {
    cat("*** iteration ",m,"\n")
    cat("acceptance probabilities:\n")
    cat("jump   => ", round(100*mean(alpha.jump[!is.na(alpha.jump)])), "\n")
    cat("thetaQ => ", round(100*mean(alpha.thetaQ[1:(m-1)])), "\n")
    cat("Omega  => ", round(100*mean(alpha.Omega[1:(m-1)])), "\n")
    cat("Phi  => ", round(100*mean(alpha.par[1:(m-1)])), "\n")
    cat("current model: ", pars$gamma, "\n")
    cat("current amount of non-stat: ", count.nonstat, "\n")
    }
  # browser()
  pars.save <- pars
  pars <- drawJump(pars)
  pars <- drawOmega(pars)
  pars <- drawThetaQ(pars)
  pars <- drawPhi(pars)
  pars <- drawSige2(pars)
  
  ## if non-stationary, retain previous iteration's values
  if (max(abs(eigen(pars$Phi[1:N,1:N], only.values=TRUE)$values))>=1) {
    # if (max(abs(diag(pars$Phi)))>=1) {
    pars <- pars.save
    count.nonstat <- count.nonstat + 1
  }
  
  gamma.i[m,] <- pars$gamma
  lambda.i[m,] <- pars$lambda
  Omega.i[m,] <- as.vector(pars$Omega)
  kinfQ.i[m] <- pars$kinfQ
  lamQ.i[m,] <- pars$lamQ
  sige2.i[m] <- pars$sige2
  Phi.i[m,] <- as.vector(pars$Phi)
  Mu.i[m,] <- pars$mu
}

save(file=results.file, priors, M, N, mats, W, dates, Omega.i, lamQ.i, kinfQ.i, lambda.i, gamma.i, sige2.i, Phi.i, Mu.i, alpha.Omega, alpha.thetaQ, alpha.jump, alpha.par)
