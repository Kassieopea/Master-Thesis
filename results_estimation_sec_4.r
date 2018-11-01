## analyze results of model estimation and selection
## Tables 3, 4, 5, 6
rm(list=ls())
load("TEST_results.RData")
source("R/rrp_functions_MEF.r")
source("R/estimation_fns_MEF.r") # for estAllModels()
source("R/analysis_fns_MEF.r")
source("R/jsz_fns_MEF.r")

library(corpcor)


to.latex <- TRUE
filename.mcmc <- "estimates/THESIS_mcmc_M0_N3_20181030.RData"
# filename.mcmc <- "estimates/THESIS_mcmc_M0_N3_20180919.RData"
filename.gvs <- "estimates/THESIS_gvs_N3_20181030.RData"
filename.ssvs <- "estimates/THESIS_ssvs_N3_20181030.RData"
filename.rjmcmc <- "estimates/THESIS_rjmcmc_N3_20181030.RData"
# 
# # browser()
cat("# Table 3 - Estimates for unrestricted model\n")
# printParameterEstimates(filename.mcmc, to.latex)
# # 
cat("# Table 4 - Risk price restrictions\n")
# # # ## summary statistics for gamma
# printGammaStats(c(filename.ssvs, filename.gvs, filename.rjmcmc), to.latex)
# #
cat("# Table 5 - Posterior model probabilities\n")
all.models <- estAllModels() ## for AIC/BIC
save.image("TEST_results.RData")
# printModels(c(filename.ssvs, filename.gvs, filename.rjmcmc), all.models, to.latex)
printModels(c(filename.gvs, filename.gvs, filename.rjmcmc), all.models, to.latex)
# 
# #####posterior probability that level/slope/curve risk is priced
pricedRisks(filename.gvs)
# 
# cat("# Table 6 - Model selection and prior dispersion\n")
# # ## sensitivity to prior dispersion
printModelFreq("estimates/THESIS_gvs_10000_N3_20181030.RData")
printModelFreq("estimates/THESIS_gvs_1000_N3_20181030.RData")
printModelFreq("estimates/THESIS_gvs_N3_20181030.RData")
printModelFreq("estimates/THESIS_gvs_10_N3_20181030.RData")
