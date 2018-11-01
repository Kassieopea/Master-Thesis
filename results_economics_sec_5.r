## analyze estimation results and economic implications
## Tables 7-9
## Figures 1, 2

rm(list=ls())
source("R/jsz_fns_MEF.r")
source("R/rrp_functions_MEF.r")
source("R/analysis_fns_MEF.r")

library(corpcor)
library(mvtnorm)
library(MASS)

to.file <- TRUE
init(N=3)
# browser()
## load model estimates
models <- list(M0 = loadModel("M0", "estimates/THESIS_mcmc_M0_N3_20181030.RData"),
               M1 = loadModel("M1", "estimates/THESIS_mcmc_M1_N3_20181031.RData"),
               M2 = loadModel("M2", "estimates/THESIS_mcmc_M2_N3_20181031.RData"),
               M3 = loadModel("M3", "estimates/THESIS_mcmc_M3_N3_20181031.RData"),
               BMA = loadModel("BMA", "estimates/THESIS_gvs_N3_20181030.RData"))

###cross-sectional fit
for (model in models)
    analyzeFit(model)

printAvgFactors(models)
printAvgYieldCurves(models, sample.mean=FALSE)
#
cat("# Table 7 - Persistence and volatility\n")
# browser()
printPersVol(models, to.file)
#
## Figure 1: Term structure of volatility
plotVolatilities(models[c("M0", "BMA")], to.file)

# browser()
## Figure 2: risk-neutral yield and term premium
name.list <- c("_short", "_medium", "" )

for (i in 1:3){
  mats_plot <- mats[(1+3*(i-1))]
  name <- name.list[i]
  plotExpTP(models[c("M0", "BMA")], to.file, mats_plot, name)
}
# plotExpTP(models[c("M0", "BMA")], to.file)

cat("# Table 8 - Historical changes in long-term rates and expectations \n")
printHistChanges(models, to.file)

cat("# Table 9 - Return predictability\n")
mat.sel <- match(c(3,5,7,10)*n.per, mats)
for (m in seq_along(names(models))) {
    cat("calculating returns for", models[[m]]$name, "\n")
    models[[m]] <- calculateReturns(models[[m]])
}
printCP(models, to.file)

