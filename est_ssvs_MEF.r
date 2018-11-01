## estimate Gaussian DTSM using block-wise Metropolis-Hastings algorithm
## SSVS -- Stochastic Search Variables Selection (George, McCulloch, 1993)

rm(list=ls())
# load("Testdata_SSVS.RData")
source("R/rrp_functions_MEF.r")
source("R/jsz_fns_MEF.r")
source("R/estimation_fns_MEF.r")
set.seed(4493)
library(corpcor)
init(N=3)
pars.mle <- estML(Y, W, mats)
pars <- pars.mle

## priors
priors <- getPriors()

## get cond. posterior SD for lambda at MLE estimates
crossterm <- cbind(diag(N), matrix(rep(0, (2*N)), nrow =N))
mu.Q <- pars$loads$K0Q.cP
phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
rvar.res <- getLambdaRVAR(rep(1, N+(N*(2+N))), cP=cbind(Y%*%t(W), MF), mu.Q, phi.Q, pars$OmegaInv, crossterm)


## SSVS prior
priors$tau1 <-priors$c1 * sqrt(diag(rvar.res$cov.mat))  # for included elements of lambda
priors$tau0 <- priors$c0 * sqrt(diag(rvar.res$cov.mat)) # for excluded elements of lambda
cat("c0 =", priors$c0, " c1 =", priors$c1, "\n")

## look at prior density implied by c0 and c1
cat("c1^2/c0^2 =", priors$c1^2/priors$c0^2, "\n")

results.file <- getResultsFileName("THESIS_ssvs", N)

Scale.Omega.1 <<- 5000
Scale.Omega.2 <<-100000
Scale.Theta <<- 100
Scale.Phi <<- 25

save.image("Testdata_SSVS.RData")
########################################
## Metropolis-Hastings algorithm

## number of iterations
M <- 55000

## draws
Omega.i <- matrix(NA, M, (N+2)^2)
lambda.i <- matrix(NA, M, N+N*(N+2))
gamma.i <- matrix(NA, M, N+N*(N+2))
lamQ.i <- matrix(NA, M, N)
kinfQ.i <- matrix(NA, M, 1)
sige2.i <- matrix(NA, M, 1)
Phi.i <- matrix(NA, M, (N+2)^2)
Mu.i <- matrix(NA, M, (N+2))

## acceptance probabilities
alpha.Omega <- numeric(M)
alpha.thetaQ <- numeric(M)
alpha.par <- numeric(M)

## probability of inclusion
prob.gamma <- matrix(NA, M, N+N*(N+2))

count.nonstat <- 0
# browser()
for (m in 1:M) {
    if( m %% 1000 == 0) {
        cat("*** iteration ",m,"\n")
        cat("acceptance probabilities:\n")
        cat("thetaQ  => ", round(100*mean(alpha.thetaQ[1:(m-1)])), "\n")
        cat("Omega   => ", round(100*mean(alpha.Omega[1:(m-1)])), "\n")
        cat("Phi   => ", round(100*mean(alpha.par[1:(m-1)])), "\n")
        cat("current model: ", pars$gamma, "\n")
        cat("prob. incl.: ", round(100*colMeans(prob.gamma[(m-19):(m-1),])), "\n")
        cat("current amount of non-stat: ", count.nonstat, "\n")
    }
    # browser()
    pars.save <- pars
    pars <- drawLambdaSSVS(pars)
    pars <- drawGammaSSVS(pars)
    pars <- drawOmega(pars)
    pars <- drawThetaQ(pars)
    pars <- drawPhi(pars)
    pars <- drawSige2(pars)

    ## if non-stationary, retain previous iteration's values
    if (max(abs(eigen(pars$Phi[1:N,1:N], only.values=TRUE)$values))>=1) {
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
    Mu.i[m,]<- pars$mu
}

save(file=results.file, priors, M, N, mats, W, dates, Omega.i, lamQ.i, kinfQ.i, lambda.i, gamma.i, sige2.i, Phi.i, Mu.i, alpha.Omega, alpha.thetaQ, prob.gamma)
