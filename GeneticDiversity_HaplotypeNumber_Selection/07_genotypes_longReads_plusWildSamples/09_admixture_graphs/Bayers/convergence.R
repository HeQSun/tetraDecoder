library("coda") # We load the coda library

#This is our fraction of samples to discard as burn-in. This can be changed to
#decide what a good burn-in is
BurnInFraction = 0.35

#We format the plots
par(mfrow=c(3,3))

#We assume that the output of the 3 independently run chains are
#chain1.txt, chain2.txt, and chain3.txt
# chain2.MCMC.SampleC_Hap_3.txt
Chain1 = read.csv("chain1.MCMC.SampleC_Hap_3.txt")
Chain2 = read.csv("chain2.MCMC.SampleC_Hap_3.txt")
Chain3 = read.csv("chain3.MCMC.SampleC_Hap_3.txt")

#We extract the total branch length
Branches1 = Chain1[,which(names(Chain1) == "total_branch_length")]
Branches2 = Chain2[,which(names(Chain2) == "total_branch_length")]
Branches3 = Chain3[,which(names(Chain3) == "total_branch_length")]

#And our total branch length after a burn-in
Branches1Burn = Branches1[-(1:round(BurnInFraction*length(Branches1)))]
Branches2Burn = Branches2[-(1:round(BurnInFraction*length(Branches2)))]
Branches3Burn = Branches3[-(1:round(BurnInFraction*length(Branches3)))]

#We do the same thing for the posterior
Post1 = Chain1[,which(names(Chain1) == "posterior")]
Post2 = Chain2[,which(names(Chain2) == "posterior")]
Post3 = Chain3[,which(names(Chain3) == "posterior")]
Post1Burn = Post1[-(1:round(BurnInFraction*length(Post1)))]
Post2Burn = Post2[-(1:round(BurnInFraction*length(Post2)))]
Post3Burn = Post3[-(1:round(BurnInFraction*length(Post3)))]

#And for the number of admixture events
admixes1 = Chain1[,which(names(Chain1) == "no_admixes")]
admixes2 = Chain2[,which(names(Chain2) == "no_admixes")]
admixes3 = Chain3[,which(names(Chain3) == "no_admixes")]
admixes1Burn = admixes1[-(1:round(BurnInFraction*length(admixes1)))]
admixes2Burn = admixes2[-(1:round(BurnInFraction*length(admixes2)))]
admixes3Burn = admixes3[-(1:round(BurnInFraction*length(admixes3)))]


#We plot our trace plots (these include all samples,
#including those in the burn-in period)
posteriormax = max(c(Post1, Post2, Post3))
posteriormin =  min(c(Post1[-c(1:1000)], Post2[-c(1:1000)], Post3[-c(1:1000)]))
plot(Post1, type = "l", col ="steelblue4", xlab = "", ylab = "log-Posterior", main = "Chain 1", ylim = c(posteriormin, posteriormax))
plot(Post2, type = "l", col ="steelblue4", xlab = "", ylab = "", main = "Chain 2", ylim = c(posteriormin, posteriormax))
plot(Post3, type = "l", col ="steelblue4", xlab = "", ylab = "", main = "Chain 3", ylim = c(posteriormin, posteriormax))

plot(Branches1, type = "l", col ="limegreen",  xlab = "", ylab = "Total Branch Length")
plot(Branches2, type = "l", col ="limegreen",  xlab = "", ylab = "")
plot(Branches3, type = "l", col ="limegreen",  xlab = "", ylab = "")

plot(admixes1, pch = 19, col ="orange", xlab = "Iteration", ylab = "Number of Admixture Events", cex = 0.48)
plot(admixes2, pch = 19, col ="orange", xlab = "Iteration", ylab = "", cex = 0.48)
plot(admixes3, pch = 19, col ="orange", xlab = "Iteration", ylab = "", cex = 0.48)

#We reformat the plotting space for the Gelman-Rubin Plots
par(mfrow=c(1,1))

#This computes and plots the median of the scale reduction factor for total branch length
#after the burn-in period for a cumulative increasing number of samples
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(Branches1Burn)/50*i))
  x = c(x, (round(length(Branches1Burn)/50*i)))
 a = gelman.diag(mcmc.list(list(mcmc(Branches1Burn[seq]), mcmc(Branches2Burn[seq]), mcmc(Branches3Burn[seq]))), autoburnin = F)
y = c(y, (a[[1]][1]) )
 }
plot(x, y, type = "l", ylab = "Gelman-Rubin Scale Reduction Factor", main = "Gelman-Rubin Convergence Diagnostics",
     xlab = paste0("Cumulative Total Samples (After Burn-in of ",toString(round(BurnInFraction*length(admixes1))) , " Samples)"),
     col ="steelblue4", lwd = 2, ylim = c(0.98, 3.2))

#We do the same thing for the number of admixture events
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(admixes1Burn)/50*i))
  x = c(x, (round(length(admixes1Burn)/50*i)))
  a = gelman.diag(mcmc.list(list(mcmc(admixes1Burn[seq]), mcmc(admixes2Burn[seq]), mcmc(admixes3Burn[seq]))), autoburnin = F)
  y = c(y, (a[[1]][1]) )
}

lines(x, y, col ="orange", lwd = 2)

#And the posterior
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(Post1)/50*i))
  x = c(x, (round(length(Post1)/50*i)))
  a = gelman.diag(mcmc.list(list(mcmc(Post1Burn[seq]), mcmc(Post2Burn[seq]), mcmc(Post3Burn[seq]))), autoburnin = F)
  y = c(y, (a[[1]][1]))
}

lines(x, y, col ="limegreen", lwd = 2)

#We plot a horizontal line at 1 for reference as well as a legend
lines(c(-100, 1e10), c(1,1), col = "red", lty = 3)

legend("topright", legend=c("Posterior", "Total Branch Length", "Number of Admixtures"),
       col=c("steelblue4", "limegreen", "orange"), lwd = 2)

#We can also plot the plots of the autocorrelation after the burn-in
par(mfcol=c(2,1))
acf(Post1Burn,lag.max = 500 , ci = -1, col ="steelblue4", ylim = c(-0.2,1), main = "Autocorrelation of Posterior (Chain 1)", ylab = "Autocorrelation")
acf(Branches1Burn,lag.max = 500 , ci = -1, col ="limegreen", ylim = c(-0.2,1), main = "Autocorrelation of Total Branch Length (Chain 1)", ylab = "Autocorrelation")




