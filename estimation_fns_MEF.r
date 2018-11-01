getPriors <- function()
    list(g = 100,  ## Zellner's g-prior
         c0 = 0.01,
         c1 = 100,
         p = rep(.5, N+(N*(2+N))), ## probability of inclusion
         ## sig.kinfQ = 1/n.per/100,
         lamQ.min = -1,
         lamQ.max = 0)

checkKKT <- function(theta, obj, ...) {
    require(numDeriv)
    y <- obj(theta, ...)
    kkttol <- 10*.Machine$double.eps^(1/4)
    kkt2tol <- 100* (.Machine$double.eps^(1/4))
    ngatend <- grad(obj, theta, method="Richardson", side=NULL, method.args=list(), ...)
    cat("Gradient:")
    print(ngatend)
    kkt1 <- max(abs(ngatend)) <= kkttol*(1.0+abs(y))
    cat("kkt1 = ", kkt1, "\n")
    nhatend <- hessian(obj, theta,  method="Richardson", method.args=list(), ...)
    hev <- eigen(nhatend)$values
    cat("Eigenvalues:", hev, "\n")
    negeig <- (hev <= -kkttol*(1+abs(y)))
    cat("negeig = ", negeig, "\n")
    evratio <- tail(hev, 1)/hev[1]
    cat("evratio =", evratio, "\n")
    cat("evratio requirement >", kkt2tol,"\n")
    kkt2 <- (evratio > kkt2tol) && (!negeig)
    cat("kkt2 =", kkt2, "\n")
}

getOptim <- function(theta, obj, ...) {
    require(ucminf)
    obj <- match.fun(obj)
    cat('Starting optimization...\n')
    cat("Starting values:\n")
    # print(theta)
    cat("LLK at starting values:", -obj(theta, ...), "\n")

    cat("1) optimization with Nelder-Mead\n")
    # print(theta)
    # print("THETA BEFORE FIRST METHOD")
    # 
    # browser()
    rval <- optim(theta, obj, gr=NULL, ..., control=list(maxit=1000000)) ## , parscale=myparscale))
    # rval <- optim(theta, obj, gr=NULL, ..., method = "SANN", hessian = TRUE)
    theta <- rval$par
    if (rval$convergence>0)
        cat(rval$convergence, "\n")
        if(rval$convergence == 10){
          rval <- optim(theta, obj, gr=NULL, ..., method = "SANN", hessian = TRUE)
          cat("Used the SANN method instead ", rval$convergence, "\n")
        }
        # 
        #stop("optimization not converged")

    cat("2) optimization with gradient-based algorithm\n")
    rval <- optim(theta, obj, gr=NULL, ..., method = "L-BFGS-B", hessian = TRUE) ## , control=list(parscale=myparscale))
    #rval <- optim(theta, obj, gr=NULL, ..., method = "SANN", hessian = TRUE)
    if (rval$convergence>0) {
        print(rval)
        warnings("optimization not converged - using result of Nelder-Mead")
    } else {
      # print("THETA IS FROM SECOND METHOD")
      theta <- rval$par
    }
    cat("LLK at optimum:", -obj(theta, ...), "\n")
    checkKKT(theta, obj, ...)

    theta
}

getCondPostLambda <- function(rvar.res) {
    ## Calculate moments of conditional posterior distribution for lambda
    ## Arguments:
    ##  rvar.res - list with results from restricted VAR estimation
    ## Globals: priors, N
    # browser()
    post.var <- solve( diag(1/priors$lambda.sd^2) + rvar.res$inv.cov.mat )
    post.mean <- post.var %*% rvar.res$inv.cov.mat %*% rvar.res$lambda.hat

    out.matrix <- matrix(NA, N+(N*(2+N)), 5)
    rownames(out.matrix) <- c(sapply(1:N, function(x) paste('lam0_', as.character(x), sep='')), sapply( rep((1:N),N+2)*10 + as.numeric(gl(N+2, N)), function(x) paste('Lam1_', as.character(x), sep='')))
    colnames(out.matrix) <- c('mean', 't-stat', 'LB', 'UB', 'sig 5%')
    out.matrix[,1] <- post.mean
    out.matrix[,2] <- post.mean/sqrt(diag(post.var))
    out.matrix[,3] <- post.mean-1.96*sqrt(diag(post.var))
    out.matrix[,4] <- post.mean+1.96*sqrt(diag(post.var))
    out.matrix[,5] <- sign(out.matrix[,3])*sign(out.matrix[,4])==1
    print(round(out.matrix,digi=4))

    return(list(mean=post.mean, var=post.var))
}

############################################################################
############################################################################
### MLE

obj.mle <- function(theta, Y, W, mats, gamma) {
    ## objective function is sum of negative log likelihoods
    # browser()
    pars <- theta2pars(theta)
    # 
    # cat(pars$lamQ, "\n")
    
    ## check restrictions on parameter space
    valid <- TRUE
    ## diagonal elements of Sigma positive and bounded away from zero
    if (any(diag(pars$Sigma)<1e-7)) valid <- FALSE
    if (any(diag(pars$Sigma)>1)) valid <- FALSE
    ## eigenvalues of Phi.Q not explosive
    if (any(pars$lamQ>0)) valid <- FALSE

    if (any(abs(pars$dlamQ> 9999999))) valid <- FALSE
    ## eigenvalues sorted
    if (any(pars$dlamQ>0)) valid <- FALSE
    # valid <- TRUE
    ## if parameters satisfy restriction on param space
    if (valid) {
        ## evaluate likelihood function and return sum of negative logliks
        # browser()
        if (missing(gamma)) {
            warning("did not provide gamma")
            res.llk <- jsz.llk(Y, W, ZP= ZP, K1Q.X=diag(pars$lamQ), Sigma.cP=pars$Omega, mats=mats, dt=1)
        } else {
            # browser()
            res.llk <- jsz.llk(Y, W,ZP = ZP, K1Q.X=diag(pars$lamQ), Sigma.cP=pars$Omega, Sigma.cQ = pars$Omega[1:N,1:N], mats=mats, dt=1, restr=1, ind.restr=gamma)
        }

        return(sum(res.llk$llk))
    } else {
        ## else return penalty value
        return(1e6)
    }
}

scale.dlamQ <- 100
scale.Sigma <- 50000
scale.MEF <- 1

theta2pars <- function(theta) {
    ## convert theta vector to list of individual parameters
    ## Globals: uses N
    ## Q parameters
    # browser()
    #  
    pars <- list(dlamQ=theta[1:N]/scale.dlamQ)
    pars$lamQ=cumsum(pars$dlamQ)
    ## P-innovation covariance matrix
    # browser()
    pars$Sigma <- matrix(0,N_new,N_new)
    pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)] <- tail(theta,-N)/scale.Sigma
    pars$Sigma[(N+1):(N+2),(N+1):(N+2)] <- pars$Sigma[(N+1):(N+2),(N+1):(N+2)]/scale.MEF
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)
    # 
    return(pars)
}

pars2theta <- function(pars) {
    # 
    # browser()
    ## convert list of individual parameters to theta vector
    ## Globals: uses N
    dlamQ <- c(pars$lamQ[1],diff(pars$lamQ));
    #if (length(pars$lamQ)!=N) stop("lamQ has wrong length")
    Sigma.vec <- pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)]
    
    # 
    return(c(dlamQ*scale.dlamQ, Sigma.vec*scale.Sigma))
}

estML <- function(Y, W, mats, gamma) {
    ## Obtain MLE for affine Gaussian DTSM
    ## (plus additional model derived parameters used in MCMC sampler)
    ## Arguments:
    ##  gamma -- risk price specification: vector with length N+N^2 indicating which elements of lambda are unrestricted (1) or restricted to zero (0)
    ## Value:
    ##  pars -- list with parameter estimates

    getStartingValuesForMLE <- function(Sigma, Sigma.Q = 0) {
        ## obtain starting values for lamQ for MLE
        ##  -> random seeds for Q-eigenvalues
        ##
        ## Arguments: Sigma
        ## Value: list with starting values

        if (missing(Sigma))
            error("Sigma needs to be provided")
        nSeeds <- 100;  # how many random seeds to try
        best.llk <- Inf
        Omega <- Sigma %*% t(Sigma)
        #Omega_Q <- Sigma.Q %*% t(Sigma.Q)
        Omega_Q <- Omega[1:N, 1:N]  
        for (i in 1:nSeeds) {
          ##    lamQ <- -sort(abs(.01*rnorm(N)))
          (lamQ <- -abs(sort(runif(N, 0, .1))))
          res.llk <- jsz.llk(Y, W, ZP, K1Q.X=diag(lamQ), Sigma.cQ = Omega_Q, Sigma.cP=Omega, mats=mats, dt=1, restr=1, ind.restr=gamma)
          llk <- sum(res.llk$llk)
          if (llk<best.llk) {
            cat('Improved seed llk to ', llk, '\n')
            best.llk <- llk
            best.lamQ <- lamQ
          }
        }
        return(list(lamQ=best.lamQ, Sigma=Sigma, Sigma.Q = Sigma.Q))
    }
    # browser()
    cP <- Y %*% t(W[1:N,])
    cP.MEF <- Y %*% t(W)
    ZP <<- cbind(cP, MF)
    N <<- nrow(W)
    N_new <<- ncol(ZP)
    #N <<- N_new
    J <- ncol(Y)
    # if (missing(mats)||length(mats)!=J)
    #   stop("estML: mats needs to be provided and have length consistent with yield data")
    # if (missing(W)||ncol(W)!=J)
    #   stop("estML: W needs to be provided and have dimensions consistent with yield data")
    if (missing(gamma))
      gamma <- rep(1, N+(N*(2+N)))
    
    cat("*** MLE ***\n")
    
    ## (1) estimate VAR -- just to get Sigma.hat
    lm <- ar.ols(ZP, aic=FALSE, order.max = 1,intercept=TRUE, demean=FALSE)
    
    Omega.hat <- lm$var.pred
    Sigma.hat <- t(chol(Omega.hat))

    Sigma.hat.Q <- Sigma.hat[1:N,1:N]
    ## (2) numerical optimization
    pars.start <- getStartingValuesForMLE(Sigma.hat, Sigma.hat.Q)
    theta.start <- pars2theta(pars.start)

    theta <- getOptim(theta.start, obj.mle, Y, W, mats, gamma)
    pars <- theta2pars(theta)  ## lamQ, Omega
    ## value of likelihood
    #  
    # browser()
    res.llk <- jsz.llk(Y, W, ZP, K1Q.X=diag(pars$lamQ), Sigma.cQ = pars$Omega[1:N,1:N], Sigma.cP=pars$Omega, mats=mats, dt=1, restr=1, ind.restr=gamma)
    pars$kinfQ <- res.llk$K0Q.X[which(res.llk$K0Q.X!=0)]
    pars$mu <- res.llk$K0P.cP
    #  
    pars$Phi <- res.llk$K1P.cP + diag(N+2)
    pars$sige2 <- res.llk$sigma.e^2
    # browser()
    res.llkQ <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, pars$Omega[1:N,1:N], pars$sige2)
    pars$loads <- res.llkQ$loads
    # pars$loads$K0Q.cP <- pars$loads$K0Q.cP/1000
    # pars$loads$K0Q.cP <- pars$loads$K0Q.cP
    #  
    # browser()
    pars$lambda <- getRP(pars$loads, pars$mu, pars$Phi)$lambda
    
    ## remaining starting values needed
    pars$gamma <- gamma
    pars$lambda[pars$gamma==0] <- 0  ## make elements EXACTLY zero
    pars$OmegaInv <- solve(pars$Omega)
    
    pars$errors <- res.llkQ$errors
    pars$llkQ <- res.llkQ$llkQ
    
    pars$llkP <- getLlkP(t(ZP), pars$mu, pars$Phi, pars$Omega, pars$OmegaInv)
    # browser()
    return(pars)
}

############################################################################
############################################################################

drawLambda <- function(pars) {
  ## prices of risk
  ## lambda = (lambda_0, lambda_1)
  ## draw jointly using Gibbs step
  ##
  ## Arguments:
  ##  pars - list with current model parameters
  ##
  ## Value:
  ##  pars - list with updated parameters
  ##
  ## Globals: Y, W, priors
  # browser()
  cP <- Y %*% t(W)
  ZP <- cbind(cP, MF)
  # browser()    # 
  if (sum(pars$gamma)>0) {
    ## at least one parameter unrestricted
    crossterm <- cbind(diag(N), matrix(rep(0, (2*N)), nrow =N))
    phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
    rvar.res <- getLambdaRVAR(pars$gamma, ZP, pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
    if (sum(pars$gamma)>1) {
      ## more than one parameter unrestricted
      lambda.var <- solve(diag(1/priors$lambda.sd[pars$gamma==1]^2) + rvar.res$inv.cov.mat)
    } else {
      ## only one parameter unrestricted -- scalar
      lambda.var <- 1/(1/priors$lambda.sd[pars$gamma==1]^2 + rvar.res$inv.cov.mat)
    }
    lambda.mean <- lambda.var %*% rvar.res$inv.cov.mat %*% rvar.res$lambda.hat
    lambda.draw <- lambda.mean + t(chol(lambda.var)) %*% rnorm(sum(pars$gamma))
        pars$lambda <- rvar.res$R %*% lambda.draw
        # browser()
        rvar.res <- getLambdaRVAR(pars$gamma, ZP, pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
        # if(any(pars$lambda>1)){browser()}
    } else {
        pars$lambda <- numeric(N+(N*(2+N)))
    }
    # browser()
    Pdyn <- getPdyn(pars$loads, pars$lambda)
    pars$mu[1:N] <- Pdyn$mu; pars$Phi[1:N,] <- Pdyn$Phi
    # browser()
    # pars$llkP <- getLlkP(t(ZP), Pdyn$mu, Pdyn$Phi, pars$Omega, pars$OmegaInv)
    pars$llkP <- getLlkP(t(ZP), pars$mu, pars$Phi, pars$Omega, pars$OmegaInv)

    return(pars)
}

############################################################################
############################################################################

drawLambdaSSVS <- function(pars) {
    ## prices of risk
    ## lambda = (lambda_0, lambda_1)
    ## draw jointly using Gibbs step
    ## given gamma = current parameter restrictions
    ## hierarchical prior according to SSVS
    ##
    ## Arguments:
    ##  pars - list with current model parameters
    ##
    ## Value:
    ##  pars - list with updated parameters
    ##
    ## Globals: Y, W, N, alpha.lambda, m
    cP <- Y %*% t(W)
    crossterm <- cbind(diag(N), matrix(rep(0, (2*N)), nrow =N))
    phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
    rvar.res <- getLambdaRVAR(rep(1, N+(N*(2+N))), ZP, pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
    lambda.hat <- rvar.res$lambda.hat ## least-squares estimate

    ## prior variance
    ## (I) prior conditional independence, use tau0 and tau1 as prior standard deviation
    if ("tau0" %in% names(priors)) {
        ## choose high or low prior SD, depending on whether included or excluded
        tau <- priors$tau0*(1-pars$gamma) + priors$tau1*(pars$gamma)
        ## invert and put into matrix D^-1
        D.inv <- diag(1/tau)
        R <- diag( N+(N*(2+N)) ); ## prior conditional independence
        R.inv <- R
    } else {
        stop("need to focus on case with prior conditional independence")
        ## (II) g-prior
        g <- priors$c0*(1-pars$gamma) + priors$c1*(pars$gamma)
        D.inv <- diag(1/g)
        pars$R.inv <- rvar.res$inv.cov.mat
    }
    DRD.inv <- D.inv %*% R.inv %*% D.inv

    if (any(is.na(DRD.inv)))
        stop("some elements of DRD^-1 are NA!?\n")

    ## posterior variance
    lambda.var <- solve( rvar.res$inv.cov.mat + DRD.inv )
    ## posterior mean
    lambda.mean <- lambda.var %*% rvar.res$inv.cov.mat %*% lambda.hat
    ## draw from conditional posterior
    #lambda.new <- rvar.res$R %*% (lambda.mean + t(chol(lambda.var)) %*% rnorm(N*(N+1)))
    lambda.new <- lambda.mean + t(chol(lambda.var)) %*% rnorm(N+(N*(2+N)))
    
    Pdyn <- getPdyn(pars$loads, lambda.new)
    pars$lambda <- lambda.new
    pars$mu[1:N] <- Pdyn$mu; pars$Phi[1:N,] <- Pdyn$Phi
    backup.mu <- c(Pdyn$mu, pars$mu[(N+1):(N+2)])
    backup.phi <- rbind(Pdyn$Phi, pars$Phi[(N+1):(N+2),])
    pars$llkP <- getLlkP(t(ZP), backup.mu, backup.phi, pars$Omega, pars$OmegaInv)
   #  
    return(pars)
}

drawGammaSSVS <- function(pars) {
    ## draw vector of indicators using Gibbs step
    ## conditional posterior is Bernoulli
    ##  -> does not depend on the data
    ##
    ## Arguments:
    ##  pars - current model parameters
    ##
    ## Value:
    ##  new model parameters (only gamma changed)
    ##
    ## Globals: N, priors, m, priors, prob.gamma (probability of inclusion saved for tracking)
    ## Side effects: changes prob.gamma[m,]

    ## assumptions:
    ## - prior on gamma is independent Bernoulli
    ## - no other prior parameters depend on gamma


    for (i in sample.int(N+(N*(2+N)))) {
        ## sample elements of gamma consecutively, in random order
        if ("tau0" %in% names(priors)) {
            ## (I) lambda | gamma is independent (R = I)
            a <- priors$tau1[i]^-1*exp(-.5*pars$lambda[i]^2/priors$tau1[i]^2)*priors$p[i]
            b <- priors$tau0[i]^-1*exp(-.5*pars$lambda[i]^2/priors$tau0[i]^2)*(1-priors$p[i])
        } else {
            ## (II) lambda | gamma, Omega is multivariate normal
            stop("focus on prior conditional independence")
            gamma.1 <- pars$gamma; gamma.1[i] <- 1
            g1 <- priors$c0*(1-gamma.1) + priors$c1*(gamma.1)
            D1.inv <- diag( 1/g1 )
            D1RD1.inv <- D1.inv %*% pars$R.inv %*% D1.inv
            gamma.0 <- pars$gamma; gamma.0[i] <- 0
            g0 <- priors$c0*(1-gamma.0) + priors$c1*(gamma.0)
            D0.inv <- diag( 1/g0 )
            D0RD0.inv <- D0.inv %*% pars$R.inv %*% D0.inv
            a <- sqrt(det(D1RD1.inv))*exp(-.5*t(pars$lambda)%*%D1RD1.inv%*%pars$lambda)*priors$p[i]
            b <- sqrt(det(D0RD0.inv))*exp(-.5*t(pars$lambda)%*%D0RD0.inv%*%pars$lambda)*(1-priors$p[i])
        }
        prob.gamma[m,i] <<- a/(a+b)
        pars$gamma[i] <- runif(1)< (a/(a+b))
    }

    if (any(is.na(pars$gamma))) {
        cat("some gamma's are NA!?\n")
    }
    return(pars)
}

############################################################################
############################################################################

drawGVS <- function(pars) {
    ## draw pairs of (lambda_i, gamma_i)
    ## using Gibbs Variable Selection
    ##
    ## Arguments:
    ##  pars - current model parameters
    ##
    ## Value:
    ##  new model parameters
    ##
    ## Globals: Y, W, priors, prob.gamma, m (probability of inclusion saved for tracking)
    ## Side effects: changes prob.gamma[m,]

    ## assumptions:
    ## - prior on gamma is independent Bernoulli
    ## - no other prior parameters depend on gamma
    ## - conditional prior independence of lambda

    cP <- Y %*% t(W)
    ## 1. draw included elements from conditional posterior
    if (sum(pars$gamma)>0) {
        ## at least one parameter unrestricted
        # browser()
        phi.Q <- cbind(pars$loads$K1Q.cP, matrix(rep(0, (2*N)), nrow = N))
        rvar.res <- getLambdaRVAR(pars$gamma, ZP, pars$loads$K0Q.cP, phi.Q, pars$OmegaInv, crossterm)
        if (sum(pars$gamma)>1) {
            # browser()
            # more than one parameter unrestricted
            lambda.var <- solve(diag(1/priors$lambda.sd[pars$gamma==1]^2) + rvar.res$inv.cov.mat)
        } else {
            ## only one parameter unrestricted -- scalar
            lambda.var <- 1/(1/priors$lambda.sd[pars$gamma==1]^2 + rvar.res$inv.cov.mat)
        }
        lambda.mean <- lambda.var %*% rvar.res$inv.cov.mat %*% rvar.res$lambda.hat
        lambda.draw <- lambda.mean + t(chol(lambda.var)) %*% rnorm(sum(pars$gamma))
        pars$lambda <- rvar.res$R %*% lambda.draw
    } else {
        pars$lambda <- numeric(N+N*(N+2))
    }
    pars$lambda.all <- pars$lambda
    Pdyn <- getPdyn(pars$loads, pars$lambda)
    pars$mu[1:N] <- Pdyn$mu; pars$Phi[1:N,] <- Pdyn$Phi
    
    pars$llkP <- getLlkP(t(ZP), pars$mu, pars$Phi, pars$Omega, pars$OmegaInv)
    ## 2. draw excluded elements from pseudo-prior
    pars$lambda.all[pars$gamma==0] <- priors$lambda.pseudo.mean[pars$gamma==0] + priors$lambda.pseudo.sd[pars$gamma==0]*rnorm(sum(pars$gamma==0))

    ## 3. draw elements of gamma, in random order
    for (i in sample.int(N+(N*(2+N)))) {
    # for (i in sample.int(N+(N*(N)))) {
        ## (ii) draw gamma_i | lambda_i
        lambda.incl <- pars$lambda
        lambda.incl[i] <- pars$lambda.all[i]
        Pdyn.incl <- getPdyn(pars$loads, lambda.incl)
        backup.incl.mu <- c(Pdyn.incl$mu, pars$mu[(N+1):(N+2)])
        backup.incl.phi <- rbind(Pdyn.incl$Phi, pars$Phi[(N+1):(N+2),])
        
        lambda.excl <- pars$lambda
        lambda.excl[i] <- 0
        Pdyn.excl <- getPdyn(pars$loads, lambda.excl)
        backup.excl.mu <- c(Pdyn.excl$mu, pars$mu[(N+1):(N+2)])
        backup.excl.phi <- rbind(Pdyn.excl$Phi, pars$Phi[(N+1):(N+2),])
        

        ## log-likelihood given gamma_i = 1
        llkP.incl <- getLlkP(t(ZP), backup.incl.mu, backup.incl.phi, pars$Omega, pars$OmegaInv)
        ## log-likelihood given gamma_i = 0
        llkP.excl <- getLlkP(t(ZP), backup.excl.mu, backup.excl.phi, pars$Omega, pars$OmegaInv)

        ## prior of lambda_i (gamma_i=1) / pseudo-prior (gamma_i=0)
        prior.ratio <- priors$lambda.sd[i]^-1*exp(-.5*pars$lambda.all[i]^2/priors$lambda.sd[i]^2)/priors$lambda.pseudo.sd[i]^-1/exp(-.5*(pars$lambda.all[i]-priors$lambda.pseudo.mean[i])^2/priors$lambda.pseudo.sd[i]^2)

        ## P(gamma(i)=1|Y,gamma(-i),lambda,...) / P(gamma(i)=0|Y,gamma(-i),lambda,...)
        q <- exp(llkP.incl-llkP.excl)*prior.ratio*priors$p[i]/(1-priors$p[i])
        incl.prob <- ifelse(is.infinite(q), 1, q/(q+1))
        if(is.nan(q)){
          incl.prob <- 0
          # browser()
        }
        if (runif(1) < incl.prob) {
            ## include element i
            pars$gamma[i] <- 1
            pars$lambda <- lambda.incl
            pars$mu[1:N] <- Pdyn.incl$mu; pars$Phi[1:N,]<- Pdyn.incl$Phi
            pars$llkP <- llkP.incl
        } else {
            ## exclude element i
            pars$gamma[i] <- 0
            pars$lambda <- lambda.excl
            pars$mu[1:N] <- Pdyn.excl$mu; pars$Phi[1:N,] <- Pdyn.excl$Phi
            pars$llkP <- llkP.excl
        }
        # if (any(pars$Phi == 0)){
        #   browser()
        # }
        prob.gamma[m,i] <<- incl.prob  ## P(gamma_i=1 | Y, ...)
    }
    # 
    return(pars)
}


############################################################################
############################################################################

getLambda.i <- function(i, pars) {
    ## Globals: Y, W, N
    cP <- Y %*% t(W)
    T <- nrow(cP) - 1

    ## draw lambda[i], given all other values of lambda/gamma
    ## -> approach suggested in Geweke and Kuo/Mallick to draw individual elements of lambda
    j <- ((i-1) %% N)+1 ## which equation
    ## get error variance
    sigm2 <- pars$Omega[j,j] - t(pars$Omega[j,-j]) %*% pseudoinverse(pars$Omega[-j,-j]) %*% pars$Omega[j,-j]
    ## dependent variable in equation j
    ydat <- cP[2:(T+1),j] - rep(pars$loads$K0Q.cP[j],T) - cP[1:T,]%*%(pars$loads$K1Q.cP+diag(N))[j,]
    ## subtract out other terms that we condition on
    ind.row <- matrix(1:(N+N*(N+2)), N, N+3)[j,]
    ind.other <- setdiff(ind.row, i)
    lam.other <- pars$lambda[ind.other]
    xdat.row <- cbind(1, ZP[1:T,]) ## all regressors
    xdat.other <- xdat.row[,match(ind.other, ind.row)] ## only other regressors
    zdat <- ydat - xdat.other %*% lam.other
    ## OLS estimate
    xdat <- xdat.row[,match(i, ind.row)]
    b.ols <- crossprod(zdat, xdat)/crossprod(xdat)
    ## posterior mean and variance
    omega2 <- sigm2/crossprod(xdat)
    sig2.post <- 1/(1/omega2 + 1/priors$lambda.sd[i]^2)
    lami.post <- sig2.post*b.ols/omega2 ## since prior mean is zero
    ## draw from posterior
    return(list(lami = lami.post + sqrt(sig2.post)*rnorm(1), lami.post=lami.post, sig2.post=sig2.post))
}

drawJump <- function(pars) {
    ## RJMCMC -- jump step
    ## given gamma = current parameter restrictions
    ##
    ## Arguments:
    ##  pars - list with current model parameters
    ##
    ## Value:
    ##  pars - list with updated parameters
    ##
    ## Globals: Y, W, N, alpha.jump, m
    ## Side effects: alpha.jump[m]

    cP <- Y %*% t(W)

    ## draw proposed model
    p.null <- 0.25  ## probability of null move
    if (runif(1)<p.null) {
        ## within-model move
        pars <- drawLambda(pars)
    } else {
        ## jump to other model
        gamma.prop <- pars$gamma
        ind <- sample.int(N+(N*(N+2)),1)    ## pick one element randomly

        post.moments <- getLambda.i(ind, pars)
        u.mean <- post.moments$lami.post
        u.sd <- sqrt(post.moments$sig2.post)
        if (pars$gamma[ind]==0) {
            ## switch on -- dim(lambda) < dim(lambda')
            gamma.prop[ind] <- 1
            ## lambda' = g(lambda, u)  -- u is random scalar
            lambda.prop <- pars$lambda
            if (lambda.prop[ind]!=0)
                stop("should be zero")
            u <- u.mean + u.sd*rnorm(1)
            lambda.prop[ind] <- u

            ## ratio of priors is f(u) (other priors cancel out)
            ratio.prior <- priors$lambda.sd[ind]^-1*exp(-.5*u^2/priors$lambda.sd[ind]^2)

            ## ratio of proposals is 1/q(u) (jump proposal cancels out)
            ratio.prop <- 1/( u.sd^-1*exp(-.5*(u - u.mean)^2/u.sd^2) )
        } else {
            ## switch off -- dim(lambda) > dim(lambda')
            gamma.prop[ind] <- 0
            lambda.prop <- pars$lambda
            lambda.prop[ind] <- 0
            u.prime <- pars$lambda[ind]

            ## ratio of priors is 1/f(u) (other priors cancel out)
            ratio.prior <- 1/(priors$lambda.sd[ind]^-1*exp(-.5*u.prime^2/priors$lambda.sd[ind]^2))

            ## ratio of proposals is q(u) (jump proposal cancels out)
            ratio.prop <- u.sd^-1*exp(-.5*(u.prime - u.mean)^2/u.sd^2)
        }

        ## acceptance probability
        Pdyn.prop <- getPdyn(pars$loads, lambda.prop);
        backup.mu <- c(Pdyn.prop$mu, pars$mu[(N+1):(N+2)])
        backup.phi <- rbind(Pdyn.prop$Phi, pars$Phi[(N+1):(N+2),])
        llkP.prop <- getLlkP(t(ZP), backup.mu, backup.phi, pars$Omega, pars$OmegaInv)
        llr.P <- llkP.prop - pars$llkP
        alpha.jump[m] <<- min(exp(llr.P)*ratio.prior*ratio.prop, 1)

        if (runif(1)<alpha.jump[m]) {
            pars$gamma <- gamma.prop
            pars$lambda <- lambda.prop
            pars$mu[1:N] <- Pdyn.prop$mu; pars$Phi[1:N,] <- Pdyn.prop$Phi
            pars$llkP <- llkP.prop
        }
    }
    return(pars)
}

############################################################################
############################################################################

drawOmega <- function(pars) {
    ## draw Omega -- variance-covariance matrix
    ## draw Sigma using Independence Proposal
    ##
    ## Arguments:
    ##  pars - list with current model parameters
    ##
    ## Value:
    ##   pars - list with updated model parameters
    ##
    ## Globals: N, Y, W, mats, alpha.Omega, m
    ## Side effects: changes alpha.Omega[m]
    require(mvtnorm)  # for rmvt
    require(numDeriv) # for hessian
    # browser()
    #pars$kinfQ <- pars$kinfQ/1000
    
    cP <- Y %*% t(W)
    ltri <- lower.tri(matrix(NA, N+2, N+2), diag=TRUE)
    obj.Sigma <- function(vSigma) {
        ## value of neg. log cond. posterior -- due to flat prior this
        ## is just log likelihood (unless prior restrictions on
        ## eigenvalues are violated)
        # browser()
        Sigma <- matrix(0,N+2,N+2)
        Sigma[ltri] <- vSigma  ##/scale.Sigma
        Omega <- Sigma %*% t(Sigma)
        res.llkQ <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, Omega[1:N,1:N], pars$sige2)
        llkQ <- res.llkQ$llkQ
        Pdyn <- getPdyn(res.llkQ$loads, pars$lambda);
        backup.mu <- c(Pdyn$mu, pars$mu[(N+1):(N+2)])
        backup.phi <- rbind(Pdyn$Phi, pars$Phi[(N+1):(N+2),])
        llkP <- getLlkP(t(ZP), backup.mu, backup.phi, Omega, pseudoinverse(Omega))
        # browser()
        -(llkQ + llkP)
    }
    obj.Sigma1 <- function(vSigma) {
      ## value of neg. log cond. posterior -- due to flat prior this
      ## is just log likelihood (unless prior restrictions on
      ## eigenvalues are violated)
      # browser()
      Sigma <- matrix(0,N+2,N+2)
      Sigma1 <- matrix(0, N,N)
      ltri.1 <- ltri[1:N,1:N]
      Sigma1[ltri.1] <- vSigma
      
      Sigma2 <- matrix(0, 2,(N+2))
      ltri.2 <- ltri[(N+1):(N+2),]
      
      Sigma2[ltri.2] <- pars.mle$Sigma[(N+1):(N+2),][ltri.2]
      # browser()
      Sigma <- rbind(cbind(Sigma1, matrix(rep(0,N*2), nrow = N)), Sigma2)
      # Sigma[ltri[1:N,1:N]] <- vSigma  ##/scale.Sigma
      # Sigma[ltri[(N+1):(N+2),]] <- pars.mle$Sigma[(N+1):(N+2),]
      Omega <- Sigma %*% t(Sigma)
      res.llkQ <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, Omega[1:N,1:N], pars$sige2)
      llkQ <- res.llkQ$llkQ
      Pdyn <- getPdyn(res.llkQ$loads, pars$lambda);
      backup.mu <- c(Pdyn$mu, pars$mu[(N+1):(N+2)])
      backup.phi <- rbind(Pdyn$Phi, pars$Phi[(N+1):(N+2),])
      llkP <- getLlkP(t(ZP), backup.mu, backup.phi, Omega, pseudoinverse(Omega))
      # browser()
      -(llkQ + llkP)
    }
    obj.Sigma2 <- function(vSigma) {
      ## value of neg. log cond. posterior -- due to flat prior this
      ## is just log likelihood (unless prior restrictions on
      ## eigenvalues are violated)
      # browser()
      Sigma <- matrix(0,N+2,N+2)
      Sigma1 <- matrix(0, N,N)
      ltri.1 <- ltri[1:N,1:N]
      Sigma1[ltri.1] <- pars.mle$Sigma[1:N,1:N][ltri.1]
      
      Sigma2 <- matrix(0, 2,(N+2))
      ltri.2 <- ltri[(N+1):(N+2),]
      Sigma2[ltri.2] <- vSigma
      
      Sigma <- rbind(cbind(Sigma1, matrix(rep(0,N*2), nrow = N)), Sigma2)
      # Sigma2 <- matrix(0, 2, (N+2))
      # Sigma[ltri[1:N,1:N]] <-  pars.mle$Sigma[1:N,1:N]  ##/scale.Sigma
      # Sigma[ltri[(N+1):(N+2),]] <- vSigma
      Omega <- Sigma %*% t(Sigma)
      res.llkQ <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, Omega[1:N,1:N], pars$sige2)
      llkQ <- res.llkQ$llkQ
      Pdyn <- getPdyn(res.llkQ$loads, pars$lambda);
      backup.mu <- c(Pdyn$mu, pars$mu[(N+1):(N+2)])
      backup.phi <- rbind(Pdyn$Phi, pars$Phi[(N+1):(N+2),])
      llkP <- getLlkP(t(ZP), backup.mu, backup.phi, Omega, pseudoinverse(Omega))
      # browser()
      -(llkQ + llkP)
    }
    cP <- Y %*% t(W)

    ## independence proposal
    # browser()
    vSigma.mean <- as.numeric(pars.mle$Sigma[ltri])
    Sigma1.mean <- pars.mle$Sigma[1:N,1:N]
    vSigma1.mean <- as.numeric(Sigma1.mean[ltri[1:N,1:N]])
    
    Sigma2.mean <- pars.mle$Sigma[(N+1):(N+2),]
    vSigma2.mean <- as.numeric(Sigma2.mean[ltri[(N+1):(N+2),]])
    
    if (!("vSigma1.var" %in% names(pars)))
      pars$vSigma1.var <- makePD(pseudoinverse(hessian(obj.Sigma1, vSigma1.mean)))
    
    if (!("vSigma2.var" %in% names(pars)))
      pars$vSigma2.var <- makePD(pseudoinverse(hessian(obj.Sigma2, vSigma2.mean)))
  
    ## variance: use Hessian of cond. posterior
    # if (!("vSigma.var" %in% names(pars)))
    #     pars$vSigma.var <- makePD(pseudoinverse(hessian(obj.Sigma, vSigma.mean)))
    nu <- 5

    # scale.mat <- pars$vSigma.var*(nu-2)/nu
    
    scale.mat1 <- (pars$vSigma1.var*(nu-2)/nu)/Scale.Omega.1
    vSigma1.prop <- as.numeric(vSigma1.mean + rmvt(1, sigma=scale.mat1, df=nu))
    
    scale.mat2 <- (pars$vSigma2.var*(nu-2)/nu)/Scale.Omega.2

    vSigma2.prop <- as.numeric(vSigma2.mean + rmvt(1, sigma=scale.mat2, df=nu))

    vSigma.prop <- c(vSigma1.prop, vSigma2.prop)
    
    Sigma1.prop <- matrix(0, N,N)
    ltri.1 <- ltri[1:N,1:N]
    Sigma1.prop[ltri.1] <- vSigma1.prop
    
    Sigma2.prop <- matrix(0, 2,(N+2))
    ltri.2 <- ltri[(N+1):(N+2),]
    Sigma2.prop[ltri.2] <- vSigma2.prop
    
    Sigma.prop <- rbind(cbind(Sigma1.prop, matrix(rep(0,N*2), nrow = N)), Sigma2.prop)

    Omega.prop <- matrix(0, N+2, N+2)
  
    Omega.prop <- Sigma.prop %*% t(Sigma.prop)
    # browser()
    Omega.prop.inv <- solve(Omega.prop)

    ## calculate acceptance probability
    vSigma.current <- as.numeric(pars$Sigma[ltri])
    vSigma1.current <- pars$Sigma[1:N,1:N][ltri[1:N,1:N]]
    vSigma2.current <- pars$Sigma[(N+1):(N+2),][ltri[(N+1):(N+2),]]
    
    ## (1) prior - make sure it's a covariance matrix
    if (all(eigen(Omega.prop)$values > 0)) {
        ## (2) Q-likelihood
        llkQ.prop <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, Omega.prop[1:N,1:N], pars$sige2)
        llr.Q <- llkQ.prop$llkQ - pars$llkQ
        ## get P-dynamics implied by proposal
        Pdyn.prop <- getPdyn(llkQ.prop$loads, pars$lambda);
        backup.mu <- c(Pdyn.prop$mu, pars$mu[(N+1):(N+2)])
        backup.phi <- rbind(Pdyn.prop$Phi, pars$Phi[(N+1):(N+2),])
        ## (3) log-ratio P-likelihoods
        llkP.prop <- getLlkP(t(ZP), backup.mu, backup.phi, Omega.prop, Omega.prop.inv)
        llr.P <- llkP.prop - pars$llkP
        ## (4) log-ratio of proposals
        lr.prop1 <- dmvt(vSigma1.current, sigma = scale.mat1, df = nu, log = TRUE) -
              dmvt(vSigma1.prop, sigma = scale.mat1, df = nu, log = TRUE)
        lr.prop2 <- dmvt(vSigma2.current, sigma = scale.mat2, df = nu, log = TRUE) -
              dmvt(vSigma2.prop, sigma = scale.mat2, df = nu, log = TRUE)
        lr.prop <- lr.prop1 + lr.prop2
        alpha.Omega[m] <<- min(exp(llr.Q + llr.P + lr.prop), 1)
        a <- alpha.Omega[m]
    } else {
        alpha.Omega[m] <<- 0
    }
    if (runif(1)<alpha.Omega[m]) {
        pars$Omega <- Omega.prop
        pars$OmegaInv <- Omega.prop.inv
        pars$loads <- llkQ.prop$loads
        pars$errors <- llkQ.prop$errors
        pars$llkQ <- llkQ.prop$llkQ
        pars$mu[1:N] <- Pdyn.prop$mu; pars$Phi[1:N,] <- Pdyn.prop$Phi
        pars$llkP <- llkP.prop
    }
    pars
}

############################################################################
############################################################################

drawThetaQ <- function(pars) {
    ## draws kinfQ and lamQ using Independence Proposal
    ##
    ## Arguments:
    ##  pars - list of current model parameters
    ##
    ## Value:
    ##  pars - list of updated model parameters
    ##
    ## Globals: Y, W, mats, alpha.thetaQ, m
    require(mvtnorm)  # for rmvt
    require(numDeriv) # for hessian

    # browser()
    # 
    scale.kinfQ <- 1000

    obj.thetaQ <- function(thetaQ) {
        ## value of neg. log cond. posterior -- due to flat prior this
        ## is just log likelihood (unless prior restrictions on
        ## eigenvalues are violated)
        kinfQ <- thetaQ[1]/scale.kinfQ
        dlamQ <- thetaQ[2:(N+1)]
        lamQ <- cumsum(dlamQ)
        # 
        if (all(dlamQ<0, lamQ<priors$lamQ.max, lamQ>priors$lamQ.min)) {
            
            res.llkQ <- getLlkQ(Y, W, mats, lamQ, kinfQ, pars$Omega[1:N,1:N], pars$sige2)
            llkQ <- res.llkQ$llkQ
            Pdyn <- getPdyn(res.llkQ$loads, pars$lambda);
            backup.mu <- c(Pdyn$mu, pars$mu[(N+1):(N+2)])
            backup.phi <- rbind(Pdyn$Phi, pars$Phi[(N+1):(N+2),])
            llkP <- getLlkP(t(ZP), backup.mu, backup.phi, pars$Omega, pars$OmegaInv)
            -(llkQ + llkP)
        } else {
            1e6
        }
    }
    cP <- Y %*% t(W)

    ## independence proposal
    # dlamQ.mean <- c(pars.mle$lamQ[1], diff(pars.mle$lamQ)[1:N-1])
    dlamQ.mean <- c(pars.mle$lamQ[1], diff(pars.mle$lamQ))
    thetaQ.mean <- c(pars.mle$kinfQ * scale.kinfQ, dlamQ.mean)
    ## variance: use Hessian of cond. posterior
    # browser()
    if (!("thetaQ.var" %in% names(pars))){
        
        if (det(hessian(obj.thetaQ, thetaQ.mean)) > 10^20){
          pars$thetaQ.var <- makePD(solve(hessian(obj.thetaQ, thetaQ.mean)))
        } else {
          pars$thetaQ.var <- makePD(pseudoinverse(hessian(obj.thetaQ, thetaQ.mean)))
        }
    }    
    # 
    nu <- 5
    # scale.mat <- pars$thetaQ.var*(nu-2)/nu
    # thetaQ.prop <- as.numeric(thetaQ.mean + rmvt(1, sigma=scale.mat, df=nu))
    # browser()
    # scale.mat <- (pars$thetaQ.var[1:(N+1), 1:(N+1)]*(nu-2)/nu)
    scale.mat <- (pars$thetaQ.var[1:(N+1), 1:(N+1)]*(nu-2)/nu)/Scale.Theta
    # browser()
    thetaQ.prop <- as.numeric(thetaQ.mean[1:(N+1)] + rmvt(1, sigma = scale.mat, df= nu))

    # browser()
    kinfQ.prop <- thetaQ.prop[1]/(scale.kinfQ)
    # kinfQ.prop <- thetaQ.prop[1]/(scale.kinfQ*190)
    dlamQ.prop <- thetaQ.prop[2:(N+1)]
    # dlamQ.prop <- thetaQ.prop[2:(N+2+1)]
    lamQ.prop <- cumsum(dlamQ.prop)
    # browser()
    ## calculate acceptance probability
    dlamQ.current <- c(pars$lamQ[1], diff(pars$lamQ))
    thetaQ.current <- c(pars$kinfQ * scale.kinfQ, dlamQ.current)
    #browser()
    ## (1) prior
    if (all(dlamQ.prop<0, lamQ.prop<priors$lamQ.max, lamQ.prop>priors$lamQ.min)) {
        ## (2) Q-likelihood
        # browser()
        llkQ.prop <- getLlkQ(Y, W, mats, lamQ.prop[1:N], kinfQ.prop, pars$Omega[1:N,1:N], pars$sige2)
        llr.Q <- llkQ.prop$llkQ - pars$llkQ
        ## get P-dynamics implied by proposal
        Pdyn.prop <- getPdyn(llkQ.prop$loads, pars$lambda);
        backup.mu <- c(Pdyn.prop$mu, pars$mu[(N+1):(N+2)])
        backup.phi <- rbind(Pdyn.prop$Phi, pars$Phi[(N+1):(N+2),])
        
        ## (3) log-ratio P-likelihoods
        # llkP.prop <- getLlkP(t(ZP), Pdyn.prop$mu, Pdyn.prop$Phi, pars$Omega, pars$OmegaInv)
        llkP.prop <- getLlkP(t(ZP), backup.mu, backup.phi, pars$Omega, pars$OmegaInv)
        llr.P <- llkP.prop - pars$llkP
        ## (4) log-ratio of proposals
        
        lr.prop <- dmvt(thetaQ.current[1:(N+1)], sigma = scale.mat, df = nu, log = TRUE) -
            dmvt(thetaQ.prop, sigma = scale.mat, df = nu, log = TRUE)
        
        alpha.thetaQ[m] <<- min(exp(llr.Q + llr.P + lr.prop), 1)
        a <- alpha.thetaQ[m]
        # cat("exp value: ", alpha.thetaQ[m])
        # a2 <- lamQ.prop[1:N]
        # a3 <- kinfQ.prop
        # a4 <- pars$Omega[1:N,1:N]
        # a5 <- pars$sige2
        # a6 <- llkQ.prop$test.llkQ - pars$llkQ
        # a7 <- sum(llkQ.prop$errors)
        # a8 <- llkQ.prop$loads$AX
        # a9 <- llkQ.prop$loads$BX
        # browser()
        #llkQ.prop <- getLlkQ(Y, W, mats, lamQ.prop[1:N], kinfQ.prop, pars$Omega[1:N,1:N], pars$sige2)
    } else {
        # browser()
        alpha.thetaQ[m] <<- 0
    }
     
    if (runif(1) < alpha.thetaQ[m]) {
        pars$kinfQ <- thetaQ.prop[1]/scale.kinfQ
        pars$lamQ <- cumsum(thetaQ.prop[2:(1+N)])
        pars$loads <- llkQ.prop$loads
        pars$errors <- llkQ.prop$errors
        pars$llkQ <- llkQ.prop$llkQ
        pars$mu[1:N] <- Pdyn.prop$mu; pars$Phi[1:N,] <- Pdyn.prop$Phi
        pars$llkP <- llkP.prop
    }
    return(pars)
}

drawSige2 <- function(pars) {
    ## draw measurement error variance
    ## Gibbs step -- pooled linear regression, conjugate prior
    ##
    ## Arguments:
    ##  pars - current model parameters
    ##
    ## Value:
    ##  pars - updated model parameters
    ##
    ## Globals: Y, W
     
    J <- ncol(W)
    N <- nrow(W)

    ## prior: inverse gamma
    ## uninformative prior:
    alpha.0 <- 0
    delta.0 <- 0

    ssr <- sum(pars$errors^2)
    alpha.1 <- alpha.0 + (J-N)*nrow(Y)  ## only J-N independent linear combinations!
    delta.1 <- delta.0 + ssr

    pars$sige2 <- 1/rgamma(n=1,shape=alpha.1/2,rate=delta.1/2)
    pars$llkQ <- -.5*sum(pars$errors^2)/pars$sige2

    return(pars)
}

drawPhi <- function(pars){
  obj.Phi <- function(Phi){
    # browser()
    mu <- Phi[1:2]
    phi.vec <- matrix(Phi[(2+1):(2*(N+2+1))], nrow =2)
    res.llkQ <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, pars$Omega[1:N,1:N], pars$sige2)
    llkQ <- res.llkQ$llkQ
    Pdyn <- getPdyn(res.llkQ$loads, pars$lambda);
    backup.mu <- c(Pdyn$mu, mu)
    backup.phi <- rbind(Pdyn$Phi, phi.vec)
    llkP <- getLlkP(t(ZP), backup.mu, backup.phi, pars$Omega, pars$OmegaInv)
    -(llkQ + llkP)
  }
  # browser()
  ##independance proposal
  mu.mean <- pars.mle$mu[(N+1):(N+2)]
  phi.mean <- as.numeric(pars.mle$Phi[(N+1):(N+2),])
  par.mean <- c(mu.mean, phi.mean)
  
  if (!("Phi.var" %in% names(pars)))
      pars$Phi.var <- makePD(pseudoinverse(hessian(obj.Phi, par.mean)))
  nu <- 5
  scale.mat <- (pars$Phi.var*(nu-2)/nu)/Scale.Phi
  # scale.mat <- (pars$Phi.var*(nu-2)/nu)
  # browser()
  par.prop <- as.numeric(par.mean + rmvt(1, sigma = scale.mat, df = nu))
  # browser()
  mu.prop <- par.prop[1:2]
  phi.prop <- matrix(par.prop[(2+1):(2*(N+2+1))], nrow = 2)
  
  mu.current <- pars$mu[(N+1):(N+2)]
  phi.current <- pars$Phi[(N+1):(N+2),]
  par.current <-c(mu.current, as.numeric(phi.current))
  # browser()
  if (all(mu.prop<999999)){
    llkQ.prop <- getLlkQ(Y, W, mats, pars$lamQ[1:N], pars$kinfQ, pars$Omega[1:N,1:N], pars$sige2)
    llr.Q <- llkQ.prop$llkQ - pars$llkQ
    ## get P-dynamics implied by proposal
    Pdyn.prop <- getPdyn(llkQ.prop$loads, pars$lambda);
    backup.mu.prop <- c(Pdyn.prop$mu, mu.prop)
    # browser()
    backup.phi.prop <- rbind(Pdyn.prop$Phi, phi.prop)
    
    ## (3) log-ratio P-likelihoods
    # llkP.prop <- getLlkP(t(ZP), Pdyn.prop$mu, Pdyn.prop$Phi, pars$Omega, pars$OmegaInv)
    llkP.prop <- getLlkP(t(ZP), backup.mu.prop, backup.phi.prop, pars$Omega, pars$OmegaInv)
    llr.P <- llkP.prop - pars$llkP
    ## (4) log-ratio of proposals
    
    lr.prop <- dmvt(par.current, sigma = scale.mat, df = nu, log = TRUE) -
      dmvt(par.prop, sigma = scale.mat, df = nu, log = TRUE)
    
    alpha.par[m] <<- min(exp(llr.Q + llr.P + lr.prop), 1)
    a <- alpha.par[m]
    # browser()
  
  } else{
    alpha.par[m] <<- 0
  }
  a <- alpha.par[m]
  # browser()
  if (runif(1) < alpha.par[m]) {
    pars$loads <- llkQ.prop$loads
    pars$errors <- llkQ.prop$errors
    pars$llkQ <- llkQ.prop$llkQ
    pars$mu <- backup.mu.prop; pars$Phi <- backup.phi.prop
    # browser()
    pars$llkP <- llkP.prop
    
  }
  pars
}

estAllModels <- function(Lam0.free=TRUE, kinfQ=TRUE) {
    ## estimate ALL models using MLE
    ## - take all parameters other than gamma as given (MLE estimates)
    source("R/rrp_functions_MEF.r", local=TRUE)
    source("R/estimation_fns_MEF.r", local=TRUE)
    if (!kinfQ)
        source("R/rinfQ_fns.r")
    set.seed(4493)
    
    # browser()
    ## load data, init parameters -- creates global cariables
    init(N=3)
    cP <- Y %*% t(W)
    T <- nrow(cP)
    J <- length(mats)

    ## get MLE for unrestricted model
    ## (all parameters will be fixed at these values except for lambda/mu/Phi)
    if (kinfQ) {
        pars <- estML(Y, W, mats, rep(1, N+(N*(N+2))))
    } else {
        pars <- estML.rinfQ(Y, W, mats, rep(1, N+(N*(N+2))))
    }
    A <- pars$loads$AcP; B <- pars$loads$BcP
    
    # browser()

    Yhat <- rep(1,T)%*%A + cP%*%B
    cat("RMSE = ", 10000*n.per*sqrt( mean((Y-Yhat)^2) ),"\n")

    ## show persistence under Q
    ## (this will be the same for all models)
    PhiQ <- diag(N)+pars$loads$K1Q.cP
    cat("PhiQ[1,1] = ", PhiQ[1,1], ";",
    "maxev-Q = ", max(abs(eigen(PhiQ)$values)), ";",
    "IRF-Q = ",  irf.var1(PhiQ, 120)[120], "\n")
    if (!kinfQ)
        cat("rinfQ =", 1200*pars$rinfQ, "\n")

    cols <- 14
    # browser()
    if (Lam0.free) {
        K <- 2^18 ## all models
        models <- matrix(NA, K, cols)
        rownames(models) <- sapply(1:K, function(k) {
            gamma <- as.numeric(intToBits(k-1)[1:18])
            return(paste(gamma, collapse=""))
        })
    } else {
        K <- 2^15 ## only models with Lam0 = (1,1,1)
        models <- matrix(NA, K, cols)
        rownames(models) <- sapply(1:K, function(k) {
            gamma <- c(1,1,1, as.numeric(intToBits(k-1)[1:15]))
            return(paste(gamma, collapse=""))
        })
    }

    colnames(models) <- c("Phi11", "maxev-P", "IRF-P", "AIC", "BIC",
                          "E.r", "E.y", "sig.y", "sig.dy", "sig.dyrn", "sighat.dyrn",
                          "EP_1", "EP_2", "EP_3")
    for (k in 1:K) {
        # browser()
        if( k %% 500 == 0)
            cat("*** model", k, "out of", K, "\n")
        gamma <- as.numeric(strsplit(rownames(models)[k], "")[[1]])
        if (kinfQ) {
            # browser()
            res.llk <- jsz.llk(Y, W, ZP= ZP, K1Q.X=diag(pars$lamQ), 
                                Sigma.cP=pars$Omega, Sigma.cQ = pars$Omega[1:N,1:N], mats=mats, dt=1, restr = 1, ind.restr = gamma)
            # res.llk <- jsz.llk(Y, W, K1Q.X=diag(pars$lamQ),
            #                    Sigma.cP=pars$Omega, mats=mats, dt=1, restr=1, ind.restr=gamma)
        } else {
            res.llk <- jsz.llk.rinfQ(Y, W, K1Q.X=diag(pars$lamQ), rinfQ=pars$rinfQ,
                                     Sigma.cP=pars$Omega, mats=mats, dt=1, restr=1, ind.restr=gamma)
        }
        numparam <- 1+N+N*(N+1)/2+sum(gamma)
        models[k,4] <- 2*sum(res.llk$llk) + 2*numparam
        models[k,5] <- 2*sum(res.llk$llk) + numparam*log(nrow(Y))
        # browser()
        Phi <- res.llk$K1P.cP + diag(N+2)
        mu <- res.llk$K0P.cP
        models[k,1] <- Phi[1,1]
        models[k,2] <- max(abs(eigen(Phi[1:N,1:N])$values))
        ## continue only if stationary
        if (models[k,2]<1) {
            models[k,3] <- irf.var1(Phi, 120)[120]
            ## population moments
            # browser()
            EcP <- as.numeric(solve(diag(N) - Phi[1:N,1:N]) %*% mu[1:N])
            models[k,6] <- 1200*(res.llk$rho0.cP + crossprod(res.llk$rho1.cP, EcP))
            models[k,7] <- 1200*(A[J] + EcP%*%B[,J])
            VarcP <- matrix( solve(diag((N)^2) - kronecker(Phi[1:N,1:N], Phi[1:N,1:N]))%*%as.numeric(pars$Omega[1:N,1:N]), (N), (N))
            models[k,8] <- 1200*sqrt(t(B[,J]) %*% VarcP[1:N,1:N] %*% B[,J])
            VardcP <- (diag(N) - Phi[1:N,1:N]) %*% VarcP %*% t(diag(N) - Phi[1:N,1:N]) + pars$Omega[1:N,1:N]
            models[k,9] <- 1200*sqrt(t(B[,J]) %*% VardcP %*% B[,J])
            ## vol of risk-neutral rates
            loads.rn <- gaussian.loadings(mats, mu[1:N], Phi[1:N,1:N]-diag(N), pars$Omega[1:N,1:N], pars$loads$rho0, pars$loads$rho1)
            Brn <- loads.rn$B
            models[k,10] <- 1200*sqrt(t(Brn[,J]) %*% VardcP %*% Brn[,J])
            models[k,11] <- 1200*sqrt(t(Brn[,J]) %*% cov(cP[2:T,]-cP[1:(T-1),]) %*% Brn[,J])
            models[k, 12:14] <- EcP
        }
    }
    # browser()
    order.aic <- order(models[,4])
    order.bic <- order(models[,5])
    cat("Unrestricted and best restricted models:\n")
    print(round(models[c(K, order.bic[1:20]), 1:10], digi=4))

    ## mean risk factors
    tmp <- rbind(colMeans(cP),
                 models[c(K, order.bic[1:10]), 12:14])
    print(round(1200*tmp, digi=2))
    return(models)
}

drawPhi_Gibbs <- function(pars){
  T <- length(ZP[,1])
  Y.test <- as.numeric(ZP[2:T,])
  S_alt <- rbind(matrix(0, nrow = N+N*(N+2), ncol = 2+2*(N+2)), diag(2+2*(N+2)))
  Z_full <- t(cbind(rep(1,T-1), ZP[1:(T-1),]))
  
  V_hat <- t(S_alt) %*% kronecker((Z_full%*%t(Z_full)), pars$OmegaInv) %*% S_alt
  
  prior <- c(pars.mle$mu[(N+1):(N+2)], pars.mle$Phi[(N+1):(N+2),])
  
  lambda.var <- solve(diag(1/prior^2) + V_hat)
  lambda.hat <- solve(V_hat) %*% t(S_alt) %*% kronecker(Z_full, pars$OmegaInv) %*% Y.test
  lambda.mean <- lambda.var %*% V_hat %*% lambda.hat
  lambda.draw <- lambda.mean + t(chol(lambda.var))%*%rnorm(2+2*(N+2))
  
  new.mu <- lambda.draw[1:2]
  new.phi <- lambda.draw[3:(2+2*(N+2))]
  
  pars$mu[(N+1):(N+2)] <- new.mu
  pars$Phi[(N+1):(N+2),] <- new.phi
  
  pars$llkP <- getLlkP(t(ZP), pars$mu, pars$Phi, pars$Omega, pars$OmegaInv)
  # browser()
  return(pars)
  
}
