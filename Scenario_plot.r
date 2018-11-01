# for (i in 1:length(Scenarios[1,])){
#   plot(mats, Scenarios[i,])
#   lines(Scenarios[i,])
# }
setwd("C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017")
set.seed(4493)
source("R/jsz_fns.r")
source("R/rrp_functions.r")
source("R/analysis_fns.r")

library(corpcor)
library(mvtnorm)
library(MASS)
# models_base <- list(BMA_base = loadModel("BMA", "C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/estimates/THESIS_gvs_N3_20180919.RData"))

setwd("C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/R_MEF")
rm(list=ls())
source("R/jsz_fns_MEF.r")
source("R/rrp_functions_MEF.r")
source("R/analysis_fns_MEF.r")

library(corpcor)
library(mvtnorm)
library(MASS)


yield_curve <- c(-0.269, -0.184, 0.028, 0.338, 0.613, 0.8433, 0.935)
new_6m <- c(-0.271, -0.278, -0.271, -0.27, -0.269, -0.269)
new_5y <- c(0.3195, 0.504, 0.447, 0.385, 0.4041, 0.338)
new_10y <- c(0.913, 1.097, 1.0825, 0.959, 0.9897, 0.935)
new_data <- matrix(c(new_6m, new_5y, new_10y), nrow = 6)

name <- "GVS"
# name <- "unrestricted"

to.file <- TRUE
init(N=3)
# browser()
## load model estimates
filename_Base <- "C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/estimates/THESIS_gvs_N3_20180919.RData"
# filename_Base <- "C:/Users/wijsk/OneDrive/Documenten/R_Bauer2017/estimates/THESIS_mcmc_M0_N3_20180919.RData"

filename_1 <- "estimates/THESIS_gvs_N3_20181030.RData"
# filename_1 <- "estimates/THESIS_mcmc_M0_N3_20181030.RData"

quantiles <- c(0.1, 0.5, 0.9)
Scenarios <- matrix(0, nrow = length(mats), ncol = 2*length(quantiles)+1)
Forecasts <- array(0, c(length(mats), 12, 2*length(quantiles)))
Scenario_short <- matrix(0, nrow = length(mats), ncol = 2*length(quantiles))
Scenario_long <- matrix(0, nrow = length(mats), ncol = 2*length(quantiles))
Scenario_aver <- matrix(0, nrow = length(mats), ncol = 2*length(quantiles))

shock_curve <- matrix(0, nrow = length(mats), ncol = length(quantiles))
par_shock_lev <- matrix(0, nrow = 2*length(quantiles), 1)


for(i in 1:length(quantiles)){
  # browser()
  Scenario_list <- extremeScenarios(filename_1, 12, quantiles[i], 5,4)
  Scenario_short[,i] <- Scenario_list$score_scen
  Scenario_long[,i] <- Scenario_list$score_scen_long
  Scenario_aver[,i] <- Scenario_list$score_scen_aver
  Scenarios[,i] <- Scenario_list$scenario
  Forecasts[,,i] <- Scenario_list$scenario.indiv
  par_shock_lev[i] <- Scenario_list$par_shock
  
  
  
  Scenario_list <- extremeScenarios(filename_Base, 12, quantiles[i], 3,4 )
  Scenarios[,i+length(quantiles)] <- Scenario_list$scenario
  Forecasts[,,i+length(quantiles)] <- Scenario_list$scenario.indiv
  Scenario_short[,i+length(quantiles)] <- Scenario_list$score_scen
  Scenario_long[,i+length(quantiles)] <- Scenario_list$score_scen_long
  Scenario_aver[,i+length(quantiles)] <- Scenario_list$score_scen_aver
  par_shock_lev[i+length(quantiles)] <- Scenario_list$par_shock

  Delta.Y <- Y[2:120,] - Y[1:119,]
  sam.Mean <- colMeans(Delta.Y)
  sam.std <- apply(Delta.Y, 2, sd)
  Phi <- matrix(0, length(mats), 1)
  for (k in 1:length(mats)){
    VAR <- ar.ols(Delta.Y[,k], order = 1, aic = FALSE, demean = FALSE, intercept = TRUE)
    Phi[k] <- VAR$ar[,,]
  }
  # VAR <- ar.ols(Delta.Y, order = 1, aic = FALSE, demean = FALSE, intercept = TRUE)
  Mu.Year <- 12*sam.Mean
  Std.year <- matrix(0, length(mats), 1)
  for (k in 1:12){
    for (j in 1:12){
      Std.year <- Std.year + Phi^(abs(k-j))*sam.std^2
    }
  }
  # browser()
  shock_curve[,i]<-  Mu.Year + sqrt(Std.year)*qnorm(quantiles[i])
  # shock2 <- Mu.Year + sqrt(Std.year)*qnorm(1-quantiles[i]
  # shock_up = max(shock1, shock2)


    # distribution <- rnorm(50000, Mu.Year, Std.Year)

}

save.image(paste("Scenarios_data_",name, ".RData", sep = ""))
load(paste("Scenarios_data_",name, ".RData", sep = ""))

# browser()
Y_par <- matrix(NA, 6, length(mats))
Y_par[1:3,] <- shock_curve[1,]%*%t(rep(1,7)) + rep(1,3) %*% t(Y[120,])
Y_par[4:6,] <- shock_curve[7,]%*%t(rep(1,7)) + rep(1,3) %*% t(Y[120,])

Scen_par_short <- (Scenario_short[1,]-Y[120,1])%*%t(rep(1,7)) + rep(1,6) %*% t(Y[120,])
Scen_par_long <- (Scenario_long[7,]-Y[120,7])%*%t(rep(1,7)) + rep(1,6) %*% t(Y[120,])


# browser()
H <- 120
Scenarios[,7] <- Y[H,]




Scenarios_plot <- Scenarios*1200
# browser()
# filename <- "figures/Scenarios_unrestricted.eps"
filename <- paste("figures/Scenarios_",name , ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios_plot[,7], pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  lines(mats, Scenarios_plot[,k], col=colors[k], lwd=lwds[k], lty=ltys[k])
}
points(mats, yield_curve, col = 1, lty = 0, pch = 8)

# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs."), pch=1, lwd = lwds,col=c(2,4,1), lty = ltys,bg="white", cex = 1)
# title("Extreme Scenarios for the Unrestricted Model")
title(paste("Extreme scenarios for the ", name,  " model", sep = ""))
dev.off()


Scenarios_plot <- Scenario_short*1200
# filename <- "figures/Scenarios_short_unrestricted.eps"
filename <- paste("figures/Scenarios_short_", name, ".eps")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  lines(mats, Scenarios_plot[,k], col=colors[k], lwd=lwds[k], lty=ltys[k])
}

points(mats, yield_curve, col = 1, lty = 0, pch = 8)

# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs."), pch=1, lwd = lwds,col=c(2,4,1), lty = ltys,bg="white", cex = 1)
# title("Extreme Scenarios for the Unrestricted Model")
title(paste("Extreme scenarios for the ", name, " model"))
dev.off()

Forecasts.plot <- Forecasts*1200
Title.list <- c("6 months", "5 years", "10 years")
for (i in 1:3){
  # browser()
  index <- 1+3*(i-1)
  Forecasts.ind <- Forecasts.plot[index,,]
  filename <- paste("figures/Forecasts_", name, "_", toString(index), ".eps", sep = "")
  # filename <- paste("figures/Forecasts_unrestricted_", toString(index), ".eps", sep = "")
  postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
  cat("*** writing figure", filename, "\n")
  lwds <- c(2,2,2,2,2,2)
  ltys <- c(1,1,1,1,1,1)
  
  colors <- c(2,2,2,4,4,4)
  par(mar = c(4,4,2,1)+.1)
  par(mfrow = c(1, 1))
  yrange <- range(c(min(Forecasts.ind), max(Forecasts.ind)))
  plot(1:12, Forecasts.ind[,1], pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Month", col=colors[1], lwd=lwds[1], lty=ltys[1])
  for (k in 2:(length(quantiles)*2)){
    lines(1:12, Forecasts.ind[,k], col=colors[k], lwd=lwds[k], lty=ltys[k])
  }
  points(1:6, new_data[,i], col = 1, lty = 0, pch = 8)
  legend("topleft", legend = c("Ext. Model", "Base Model"), lwd = lwds, col=c(2,4),lty = ltys, bg="white",pch=1, cex = 1)
  title(paste("Quantiles of predictions for each month of 2018 \n with a maturity of ", Title.list[i], sep = ""))
  dev.off()
}


Scenarios_plot <- Scenario_long*1200
# filename <- "figures/Scenarios_long_unrestricted.eps"
filename <- paste("figures/Scenarios_long_", name, ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  lines(mats, Scenarios_plot[,k], col=colors[k], lwd=lwds[k], lty=ltys[k])
}

points(mats, yield_curve, col = 1, lty = 0, pch = 8)


# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs."), pch=1, lwd = lwds,col=c(2,4,1), lty = ltys,bg="white", cex = 1)
# title("Extreme Scenarios for the Unrestricted Model")
title(paste("Extreme Scenarios for the ", name, " Model"))
dev.off()


Scenarios_plot <- Scenario_aver*1200
# filename <- "figures/Scenarios_aver_unrestricted.eps"
filename <- paste("figures/Scenarios_aver_", name, ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  lines(mats, Scenarios_plot[,k], col=colors[k], lwd=lwds[k], lty=ltys[k])
}

points(mats, yield_curve, col = 1, lty = 0, pch = 8)


# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs."), pch=1, lwd = lwds,col=c(2,4,1), lty = ltys,bg="white", cex = 1)
# title("Extreme Scenarios for the Unrestricted Model")
title(paste("Extreme scenarios for the ", name, " model"))
dev.off()


Y_plot <- Y_par*1200
Scenarios_plot <- Scen_par_short*1200
# browser()
# filename <- "figures/Scenarios_par_short_unrestricted.eps"
filename <- paste("figures/Scenarios_par_short_", name, ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1,3)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  if (k <= length(quantiles)){
    lines(mats, Y_plot[k,], col = 3, lwd = lwds[k], lty = ltys[k])
  }
  lines(mats, Scenarios_plot[k,], col=colors[k], lwd=lwds[k], lty=ltys[k])
}
points(mats, yield_curve, col = 1, lty = 0, pch = 8)

# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs.", "AR(1) model"), pch=1, lwd = lwds,col=c(2,4,1,3), lty = ltys,bg="white", cex = 1)
# title("Parallel Scenarios for the Unrestricted Model")
title(paste("Parallel scenarios for the ", name, " model", sep = ""))
dev.off()



Scenarios_plot <- Scen_par_long*1200
# browser()
# filename <- "figures/Scenarios_par_long_unrestricted.eps"
filename <- paste("figures/Scenarios_par_long_", name, ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1,3)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  if (k > length(quantiles)){
    lines(mats, Y_plot[k,], col = 3, lwd = lwds[k], lty = ltys[k])
  }
  lines(mats, Scenarios_plot[k,], col=colors[k], lwd=lwds[k], lty=ltys[k])
}
points(mats, yield_curve, col = 1, lty = 0, pch = 8)

# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs.", "AR(1) model"), pch=1, lwd = lwds,col=c(2,4,1,3), lty = ltys,bg="white", cex = 1)
# title("Parallel Scenarios for the Unrestricted Model")
title(paste("Parallel Scenarios for the ", name, " Model", sep = ""))
dev.off()





Scenarios_plot <- Scen_par_long*1200
# browser()
# filename <- "figures/Scenarios_par_long_unrestricted.eps"
filename <- paste("figures/Scenarios_par_long_", name, ".eps", sep = "")
postscript(filename, width=6, height=5, horizontal=FALSE, pointsize=12)
cat("*** writing figure", filename, "\n")
lwds <- c(2,2,2,2,2,2,1)
ltys <- c(1,1,1,1,1,1,1)

colors <- c(2,2,2,4,4,4,1,3)
par(mar = c(4,4,2,1)+.1)
par(mfrow = c(1, 1))
yrange <- range(c(min(Scenarios_plot), max(Scenarios_plot)))
plot(mats, Scenarios[,7]*1200, pch = 1,type="l", ylab="Yield", ylim = yrange, xlab = "Maturity", col=colors[7], lwd=lwds[7], lty=ltys[7])
for (k in 1:(length(quantiles)*2)){
  if (k > length(quantiles)){
    lines(mats, Y_plot[k,], col = 3, lwd = lwds[k], lty = ltys[k])
  }
  lines(mats, Scenarios_plot[k,], col=colors[k], lwd=lwds[k], lty=ltys[k])
}
points(mats, yield_curve, col = 1, lty = 0, pch = 8)

# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot
legend("topleft", legend = c("Ext. Model", "Base Model", "Last Obs.", "AR(1) model"), pch=1, lwd = lwds,col=c(2,4,1,3), lty = ltys,bg="white", cex = 1)
# title("Parallel Scenarios for the Unrestricted Model")
title(paste("Parallel Scenarios for the ", name, " Model", sep = ""))
dev.off()
# yrange <- c(min(Scenarios_plot),(max(Scenarios_plot)+1))
# matplot(mats, t(Scenarios_plot), type = c("b"),pch=1,col = c(2,2,2,4,4,4,1), ylab = "Yields in Percentage points", xlab = "Horizon ", ylim = yrange) #plot

