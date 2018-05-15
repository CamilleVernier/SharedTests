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

rm(list = ls())
options(error=recover) 

library(Infusion)
library(caret)

setwd(dir = "/home/cvernier/Infusion")
deb <- Sys.time()
# setwd(dir = "/home/vernierc/Documents/GitCamille/Version locale/")
#setwd("/Users/raph/Documents/Taf/EtudiantsPostdocVisiteurs+Stages+Theses/_CamilleVernier/GitHub/SharedTests/") #pour Raph

# We load a function to run PSC simulations and return the 6 summary 
# statistics used in this study:
source("PSC_simulation_c.R") 
#IBDSimExec<-"./ibdsimV2.0-win-i386.exe" #IBDsim executable if ran on Windows

# If the user wants to run this script on a Linux platform, the IBDsim sources
# must be compiled thanks to the script "compile.sh" provided with the sources
# on IBDsim website.
# The g++ compiler must be installed on the user's system.
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


if(get_os()=="linux") {
  IBDSimExec<-"../IBDSim"
}  else if(get_os()=="osx") {
    IBDSimExec<-"../IBDSimMac"
    } else {IBDSimExec<-"../ibdsimv2.0_Win7.exe"}


gr <- 300
parsp <- init_grid(lower=c(log10theta=-2.5,log10thetaanc=-1.0,log10tau=-1),
                   upper=c(log10theta=1.5,log10thetaanc=3.5,log10tau=1),
                   nUnique=gr)


#10^parsp[,1]/mu
#10^parsp[,2]/mu
#10^parsp[,3]*10^parsp[,1]/mu

message("SIM: Nact(min/max)=",paste(min(10^parsp[,1]/mu),max(10^parsp[,1]/mu)," "),
        "\nT(min/max)=",paste(min(10^parsp[,3]*10^parsp[,1]/mu),max(10^parsp[,3]*10^parsp[,1]/mu)," "),
        "\nNanc(min/max)=",paste(min(10^parsp[,2]/mu),max(10^parsp[,2]/mu)," "))



# We compute a pseudo-observed set of summary statistics with
# with arbitrary parameters values, e.g. theta = 3, a = 3, tau = -1
#Expansion
log10thetaOBS <- -0.398
log10thetaancOBS <- 1.30
log10tauOBS <- 0.097


#Grille
hist(parsp$log10theta, breaks = 20)
abline(v = log10thetaOBS, col="red")
hist(parsp$log10thetaanc, breaks=20)
abline(v = log10thetaancOBS, col="red")
hist(parsp$log10tau, breaks=20)
abline(v = log10tauOBS, col="red")


sobs <- IBDSim_wrapper(log10theta=log10thetaOBS,log10thetaanc =log10thetaancOBS,log10tau=log10tauOBS,execName=IBDSimExec)
sobs

message("OBS: Nact=",10^log10thetaOBS/mu,
        "; T=",10^log10tauOBS*10^log10thetaOBS/mu,
        "; Nanc=",10^log10thetaancOBS/mu)

nR <- 1000
et <- proc.time()
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper",par.grid=parsp, nRealizations = nR, nb_cores = NULL)# nRealizations = 1
time <- proc.time()-et
time

# We project the summary statistics with neural networks.
# Statistics to be projected:
allstats <- c("H","VarK","K","M","varH","f")

thetaproj <- project("log10theta",stats=allstats,data=simuls,method="neuralNet")
thetaancproj <- project("log10thetaanc",stats=allstats,data=simuls,method="neuralNet")
tauproj <- project("log10tau",stats=allstats,data=simuls,method="neuralNet")

# We apply projections on simulated statistics:
corrSimuls <- project(simuls,projectors=list("THETA"=thetaproj,"THETAANC"=thetaancproj,"TAU"=tauproj))
corrSobs <- project(sobs,projectors=list("THETA"=thetaproj,"THETAANC"=thetaancproj,"TAU"=tauproj))

#load("PSC_simulations.rda") # loads the projected summary statistics
# We then infer the summary-likelihood surface:
densv <- infer_logLs(corrSimuls,stat.obs=corrSobs)
slik <- infer_surface(densv)
slik <- MSL(slik)
plot(slik, filled=TRUE) # seems to be valid ony for analyses with 2 parameters


# We iterate to improve parameter estimation:
nIterations <- 3
for(i in 1:nIterations) {
  slik <- refine(slik)
}
#Error in .calc_all_slices(object, fittedPars, color.palette) : 
# objet 'plot.axes' introuvable

# The slik object contains the ML estimates and 95% CIs:
summary(slik)

# We can also get the x% CIs (here, 90% CIs):
confint(slik,"log10theta",level=0.9)
confint(slik,"log10thetaanc",level=0.9)
confint(slik,"log10tau",level=0.9)


name  <- NULL
if (log10thetaancOBS>log10thetaOBS)
{name <- "contr"
}else{
  name <- "exp"}

setwd("/work/cvernier/")
save.image(paste ("/work/cvernier/", 
                  name, "th",log10thetaOBS, "thanc", log10thetaancOBS, "tau", log10tauOBS, 
                  "g", gr,"nR", nR, "nloc", nloc, "smplsize", samplesize, "mu", mu, deb,"sim.RData", sep="_"))

png(file = paste ("/work/cvernier/", deb, "plot1D.png", sep="_"))
plot1Dprof(slik)
dev.off()

png(file = paste ("/work/cvernier/", deb, "plot2D.png", sep="_"))
plot2Dprof(slik)
dev.off()
