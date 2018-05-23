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
nbcores <- 20
gr <- 300
nR <- 1000
n_loc=40 # number of loci
smplsize=100 # number of individuals to simulate
Mu=5e-2
#Expansion
log10thetaOBS <- - 1.5
log10thetaancOBS <- 1
log10tauOBS <- -0.2
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

log10thetaOBS
log10thetaancOBS
log10tauOBS

gr_theta <- c(-3, 2)
gr_thetaanc <- c(-1, 3)
gr_tau <- c(-2,2)
parsp <- init_grid(lower=c(log10theta=gr_theta[1],log10thetaanc=gr_thetaanc[1],log10tau=gr_tau[1]),
                   upper=c(log10theta=gr_theta[2],log10thetaanc=gr_thetaanc[2],log10tau=gr_tau[2]),
                   nUnique=gr)



#10^parsp[,1]/mu
#10^parsp[,2]/mu
#10^parsp[,3]*10^parsp[,1]/mu

message("SIM: Nact(min/max)=",paste(min(10^parsp[,1]/Mu),max(10^parsp[,1]/Mu)," "),
        "\nT(min/max)=",paste(min(10^parsp[,3]*10^parsp[,1]/Mu),max(10^parsp[,3]*10^parsp[,1]/Mu)," "),
        "\nNanc(min/max)=",paste(min(10^parsp[,2]/Mu),max(10^parsp[,2]/Mu)," "))



# We compute a pseudo-observed set of summary statistics with
# with arbitrary parameters values, e.g. theta = 3, a = 3, tau = -1



#Grille
op <- par(no.readonly = TRUE) # the whole list of settable "par" 's to be reset afterwards
par(mfrow=c(2, 2)) # les modifs de "par" pour les graphiques ci dessous
 hist(parsp$log10theta, breaks = 20)
 abline(v = log10thetaOBS, col="red")
 hist(parsp$log10thetaanc, breaks=20)
 abline(v = log10thetaancOBS, col="red")
 hist(parsp$log10tau, breaks=20)
 abline(v = log10tauOBS, col="red") 
par(mfrow=c(1, 1))
par(op) #et a la fin on remet "par" comme avant


sobs <- IBDSim_wrapper(log10theta=log10thetaOBS,log10thetaanc =log10thetaancOBS,log10tau=log10tauOBS, mu = Mu, sampleSize = smplsize, nloc=n_loc,execName=IBDSimExec)
sobs

message("OBS: Nact=",10^log10thetaOBS/Mu,
        "; T=",10^log10tauOBS*10^log10thetaOBS/Mu,
        "; Nanc=",10^log10thetaancOBS/Mu)


et <- proc.time()
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper",par.grid=parsp, nRealizations = nR, nb_cores = nbcores)# nRealizations = 1
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
densv <- infer_logLs(corrSimuls,stat.obs=corrSobs, nb_cores=nbcores)
slik <- infer_surface(densv)
slik <- MSL(slik)
plot(slik, filled=TRUE) # seems to be valid ony for analyses with 2 parameters


# We iterate to improve parameter estimation:
nIterations <- 3
ntimes_iter <- 1
nIterations_total <- ntimes_iter*nIterations
for (j in 1:ntimes_iter) 
{
  for(i in 1:nIterations) 
    {
  name_slik <- paste("slik", i, sep="")
  assign(name_slik, refine(slik, nb_cores=7))
    }
}

# Sauvegarder slik tous les 3 refine

# The slik object contains the ML estimates and 95% CIs:
summary(slik)

# We can also get the x% CIs (here, 90% CIs):
confint_theta <- confint(slik,"log10theta",level=0.9)
confint_thetaanc <- confint(slik,"log10thetaanc",level=0.9)
confint_tau <- confint(slik,"log10tau",level=0.9)


########################################################"

# Calcul durée du script
fin <- Sys.time()
duree <- fin - deb

# Type de simulation constraction, expansion, constant)
name  <- NULL
if (log10thetaancOBS>log10thetaOBS)
{name <- "contr"
}else if (log10thetaancOBS<log10thetaOBS){
  name <- "exp"
}else{
  name <- "cst"
}


setwd("/work/cvernier/")
exist <- exists("slik")

# Classement de la table selon si slik a fonctionné ou non
if (exist == "TRUE")
{
  chemin <- "/work/cvernier/ok/"
}else{
  chemin <- "/work/cvernier/bug/"
}

# Sauvegarde données Rdata
save.image(paste (chemin, 
                  name, "log10theta =",log10thetaOBS, "log10thetaanc =", log10thetaancOBS, "log10tau =", log10tauOBS, 
                  "taille grille =", gr,"nRealizations =", nR, "nloc =", n_loc, "smplsize =", smplsize, "mu =", Mu, 
                  deb,".RData", sep=" "))


# Sauvegarde des paramètres dans un fichier txt
dat <- capture.output(slik)
sobs_dat <- capture.output(sobs)
write(paste(name, "\nlog10theta =",log10thetaOBS, " Bornes theta = [", gr_theta[1], ",", gr_theta[2], "]", 
            "\nlog10thetaanc =", log10thetaancOBS, " Bornes thetaanc =[", gr_thetaanc[1], ",", gr_thetaanc[2], "]", 
            "\nlog10tau =", log10tauOBS, " Bornes tau =[", gr_tau[1], ",", gr_tau[2], "]",
            "\n\nOBS: Nact=",10^log10thetaOBS/Mu,"; T=",10^log10tauOBS*10^log10thetaOBS/Mu,"; Nanc=",10^log10thetaancOBS/Mu,
            "\n\n", sobs_dat[1], "\n", sobs_dat[2],
            "\n\nSIM: Nact(min/max)=",paste(min(10^parsp[,1]/Mu),max(10^parsp[,1]/Mu)," "),
            "\nT(min/max)=",paste(min(10^parsp[,3]*10^parsp[,1]/Mu),max(10^parsp[,3]*10^parsp[,1]/Mu)," "),
            "\nNanc(min/max)=",paste(min(10^parsp[,2]/Mu),max(10^parsp[,2]/Mu)," "),
            "\n\ntaille grille =", gr, "\nnRealizations =", nR, "\nnloc =", n_loc, "\nsmplsize =", smplsize, "\nmu =", Mu,  
            "\nNb refine =", nIterations_total, 
            "\n\nConfint theta = [", confint_theta$lowerpar[1], ",", confint_theta$upperpar[1], "]",
            "\nConfint theta anc = [", confint_thetaanc$lowerpar[2], ",", confint_thetaanc$upperpar[2], "]",
            "\nConfint tau = [", confint_tau$lowerpar[3], ",", confint_tau$upperpar[3], "]",
            "\n\n", dat[1], "\n", dat[2], "\n", dat[3], "\n", dat[4],"\n", dat[5], "\n", dat[6], "\n", dat[7], "\n", dat[8],
            "\n","\n", capture.output(duree), "\n", deb, sep=" "), 
              file=paste(deb,".txt"))
  
# Sauvegarde des graphiques (en cours)
pdf(file = paste (chemin, deb, "plot.pdf", sep="_"))
 op <- par(no.readonly = TRUE) # the whole list of settable "par" 's to be reset afterwards
 par(mfrow=c(2, 2)) # les modifs de "par" pour les graphiques ci dessous
  hist(parsp$log10theta, breaks = 20)
  abline(v = log10thetaOBS, col="red")
  hist(parsp$log10thetaanc, breaks=20)
  abline(v = log10thetaancOBS, col="red")
  hist(parsp$log10tau, breaks=20)
  abline(v = log10tauOBS, col="red") 
 par(mfrow=c(1, 1))
 par(op) #et a la fin on remet "par" comme avant
  plot(slik)
  plot1Dprof(slik)
  plot(thetaproj)
  plot(thetaancproj)
  plot(tauproj)
dev.off()


