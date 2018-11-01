## MCMC (URP) for specific model on simulated data
# 
# install.packages('MCMCpack')
# install.packages('mvtnorm')
# install.packages('numDerim')
setwd("C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/R_MEF")
library(MCMCpack)
library(mvtnorm)
library(corpcor)
#library(numDerim)

# rm(list=ls())
# load("SIM_TEST.RData")
options(error=recover)
source("R/rrp_functions_MEF.r")
source("R/jsz_fns_MEF.r")
source("R/estimation_fns_MEF.r")
source("R/Bauer_functions.r")
source("R/jsz_fns_MEF_BAUER.r")
set.seed(4493)
# load("true_pars_old.RData")
# true model -- set DGP parameters
true.pars <- initSimulation(restr = TRUE)

## restr = TRUE runs the simulation using a DGP with restricted risk prices. To instead use a DGP with unrestricted risk prices, use restr = FALSE (and change change character string in call to `getResultsFileName` to include "urp" instead of "rrp"

T <- 300
M <- 10000 	## number of iterations
B <- 5000  	## number of burn-in iterations
n <- 100

## priors
priors <- getPriors()

results.file <- getResultsFileName("THESIS_sim_rrp_mcmc", N)

## draws

# Omega.i.part1 <- matrix(NA, n*M, 12)
# Omega.i.part2 <- matrix(NA, n*M, 13)
Omega.i <- matrix(NA, n*M, (N+2)^2)
lambda.i <- matrix(NA, n*M, N+(N*(2+N)))
lamQ.i <- matrix(NA, n*M, N)
kinfQ.i <- matrix(NA, n*M, 1)
sige2.i <- matrix(NA, n*M, 1)
Phi.i <- matrix(NA, n*M, (N+2)^2)
Mu.i <- matrix(NA, n*M, (N+2))
## pars.mle <- vector(mode='list', length(n))
Ysim <- vector(mode='list', length(n))


Scale.Omega.1 <<- 1
Scale.Omega.2 <<-1
Scale.Theta <<- 1
Scale.Phi <<- 1
# save.image("SIM_TEST.RData")

ind <- 0
count_error <- 0
for (i in 1:n) {
  save.image("SIM_TEST.RData")
  cat("************************************\n")
  cat(" simulation ", i, "\n")
  cat("************************************\n")
  
  ## acceptance probabilities 
  alpha.Omega <- numeric(B+M)
  alpha.thetaQ <- numeric(B+M)
  alpha.par <- numeric(B+M)
  count.nonstat <- 0
  # browser()
  SIM <- simulateYields(true.pars, T)
  # browser()
  # SIM_Y <- simulateYields.OLD(true.pars.OLD, T)
  Y <- SIM$Y
  MF <- SIM$MF
  M_sim1 <- MF[,1]%*% MF_load1
  M_sim2<- MF[,2]%*%MF_load2
  
  Ysim[[i]] <- Y
  
  ZP <<- cbind(Y%*%t(W), M_sim1 %*% t(MF_load1),M_sim2 %*% t(MF_load2) )
  ## starting values
  # browser()
  possibleError <-  tryCatch(
    estML(Y,W,mats),
    # estML(Y,W,c(1,2,3,4)),
    error = function(e) e
  )
  # browser()
  if(inherits(possibleError, "error")) next
    
    # pars.mle <- estML(Y, W, mats)
  pars.mle <- possibleError
  
  pars <- pars.mle
  # if (max(abs(eigen(pars$Phi, only.values=TRUE)$values))>=1) {
  #   pars$Phi[4,4] <- ar.ols(ZP[,4], order = 1)$ar
  # }
  
  ## ## save ML estimates
  ## pars.mle[[i]] <- pars[c('Omega', 'Sigma', 'lamQ', 'kinfQ', 'mu', 'Phi', 'sige2', 'loads', 'lambda', 'gamma')]
  
  ## prior for lambda
  crossterm <- cbind(diag(N), matrix(rep(0, (2*N)), nrow =N))
  phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
  rvar.res <- getLambdaRVAR(rep(1, N+N*(N+2)), cP=cbind(Y%*%t(W), MF), pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
  ## pars.mle[[i]]$lambda.se <- sqrt(diag(rvar.res$cov.mat))
  priors$lambda.sd <- sqrt(priors$g)*sqrt(diag(rvar.res$cov.mat))   ## diagonalized g-prior

  ## show posterior for lambda (given opt. values for other params)
  post.mom <- getCondPostLambda(rvar.res)

  for (m in 1:(B+M)) {
      if( m %% 1000 == 0) {
          cat("*** iteration ",m,"\n")
          cat("acceptance probabilities:\n")
          cat("thetaQ => ", round(100*mean(alpha.thetaQ[1:(m-1)])), "\n")
          cat("Omega  => ", round(100*mean(alpha.Omega[1:(m-1)])), "\n")
          cat("Phi  => ", round(100*mean(alpha.par[1:(m-1)])), "\n")
          cat("current amount of non-stat: ", count.nonstat, "\n")
          cat("current amount of error in drawing: ", count_error, "\n")
      }
      # browser()
      pars.save <- pars
      pars <- drawLambda(pars)
      possibleError2 <-  tryCatch(
        drawThetaQ(pars),
        # estML(Y,W,c(1,2,3,4)),
        error = function(e) e
      )
      # browser()
      if(inherits(possibleError2, "error")) {
        count_error <- count_error + 1
        next
      }
      pars <- possibleError2
      possibleError3 <-  tryCatch(
        drawOmega(pars),
        # estML(Y,W,c(1,2,3,4)),
        error = function(e) e
      )
      # browser()
      if(inherits(possibleError3, "error")) {
        count_error <- count_error + 1
        next
      }
      pars <- possibleError3
      possibleError4 <-  tryCatch(
        drawSige2(pars),
        # estML(Y,W,c(1,2,3,4)),
        error = function(e) e
      )
      # browser()
      if(inherits(possibleError4, "error")) {
        count_error <- count_error + 1
        next
      }
      pars <- possibleError4
      pars <- drawPhi(pars)
      # pars <- drawOmega(pars)
      
      # pars <- drawSige2(pars)
      # pars <- drawPhi(pars)

      ## if non-stationary, retain previous iteration's values
      if (max(abs(eigen(pars$Phi[1:N,1:N], only.values=TRUE)$values))>=1) {
      # if (max(abs(eigen(pars$Phi, only.values=TRUE)$values))>=1) {
          pars <- pars.save
          count.nonstat <- count.nonstat + 1
      }

      if (m>B) {
          ind <- ind+1
          lambda.i[ind,] <- pars$lambda
          Omega.i[ind,] <- as.vector(pars$Omega)
          kinfQ.i[ind] <- pars$kinfQ
          lamQ.i[ind,] <- pars$lamQ
          sige2.i[ind] <- pars$sige2
          Phi.i[ind,] <- as.vector(pars$Phi)
          Mu.i[ind,] <- pars$mu
      }
  }
}

save(file=results.file, true.pars, priors, Ysim, W, B, M, N, T, n, n.per, mats, Mu.i,lamQ.i, kinfQ.i, lambda.i, sige2.i, Omega.i, Phi.i)
