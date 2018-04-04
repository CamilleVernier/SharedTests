### Inference with the summary likelihood method 
### under a population size change (PSC) model.
###
### Alexandre Gouy - alexandre.gouy@iee.unibe.ch
###
### This script allows to perform simulations using the software 
### IBDsim and to infer the model parameters using summary likelihood

# Required: this R script, "PSC_simulation.R", and
# the IBDsim sources which can be found here: 
# http://www1.montpellier.inra.fr/CBGP/software/ibdsim/download.html
# (Leblois et al. 2009 Mol. Ecol. Resources)

# Important note:
# This script is designed to perform one simulation at each call of IBDsim,
# which is not optimal in terms of computation times. In our study, 
# to save computation time, we generated the whole distribution of 
# summary statistics with one call of IBDsim (i.e. 1000 simulations
# for one set of parameters). But this procedure involves manual 
# steps, especially during the iteration procedure, to project and 
# add the new summary statistics.
# Therefore, we preferred to provide a script which is functional 
# for any user, even if it is substantially slower.

library(Infusion)
library(caret)

setwd("/Users/raph/Documents/Taf/EtudiantsPostdocVisiteurs+Stages+Theses/_CamilleVernier/GitHub/SharedTests/") #pour Raph

# If the user wants to run this script on a Linux platform, the IBDsim sources
# must be compiled thanks to the script "compile.sh" provided with the sources
# on IBDsim website.
# The g++ compiler must be installed on the user's system.
IBDSimExec<<-"./IBDSimMac" # if ran on Linux


# We load a function to run PSC simulations and return the 6 summary 
# statistics used in this study:
source("PSC_simulation.R") 
#IBDSimExec<-"./ibdsimV2.0-win-i386.exe" #IBDsim executable if ran on Windows


# We first define the parameter grid.
# The number of sets is equal to the number of empirical distribution 
# desired for the first iteration (here, 100).
parsp <- init_grid(lower=c(log10theta=-1,log10a=-3,log10tau=-4),
                   upper=c(log10theta=1.5,log10a=3,log10tau=1),
                   nUnique=100)

# We compute a pseudo-observed set of summary statistics with
# with arbitrary parameters values, e.g. theta = 1, a = 1, tau = -1
sobs <- IBDSim_wrapper(log10theta=1,log10a =-1,log10tau=-2)

# We generate the first set of empirical distributions
# Note: To avoid these computationnaly intensive first steps, 
# an initial set of projected summary statistics has been provided
# (you can skip the next lines and load the simulations l.68)
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper",par.grid=parsp,nRealizations=10, verbose=FALSE)# nRealizations = 1

# We project the summary statistics with neural networks.
# Statistics to be projected:
allstats <- c("H","VarK","K","M","varH","f") 

thetaproj <- project("log10theta",stats=allstats,data=simuls,method="neuralNet")
aproj <- project("log10a",stats=allstats,data=simuls,method="neuralNet")
tauproj <- project("log10tau",stats=allstats,data=simuls,method="neuralNet")

# We apply projections on simulated statistics:
corrSimuls <- project(simuls,projectors=list("THETA"=thetaproj,"A"=aproj,"TAU"=tauproj))
corrSobs <- project(sobs,projectors=list("THETA"=thetaproj,"A"=aproj,"TAU"=tauproj))

# load("PSC_simulations.rda") # loads the projected summary statistics
# We then infer the summary-likelihood surface:
densv <- infer_logLs(corrSimuls,stat.obs=corrSobs)
slik <- infer_surface(densv)
slik <- MSL(slik)

# We iterate to improve parameter estimation:
nIterations <- 3
for(i in 1:nIterations) {
  slik <- refine(slik)
}

# The slik object contains the ML estimates and 95% CIs:
summary(slik)

# We can also get the x% CIs (here, 90% CIs):
confint(slik,"log10theta",level=0.9)
confint(slik,"log10a",level=0.9)
confint(slik,"log10tau",level=0.9)
