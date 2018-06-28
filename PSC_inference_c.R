### Inference with the summary likelihood method 


rm(list = ls())

library(Infusion)
library(caret)


#setwd(dir = "/work/cvernier/Infusion")
deb <- Sys.time()
# setwd(dir = "/home/vernierc/Documents/GitCamille/Version locale/")
#setwd("/Users/raph/Documents/Taf/EtudiantsPostdocVisiteurs+Stages+Theses/_CamilleVernier/GitHub/SharedTests/") #pour Raph

# We load a function to run PSC simulations and return the 6 summary 
# statistics used in this study:
source("PSC_simulation_c.R") 
nbcores <- 25
gr <- 300
nR <- 1000
n_loc=100 # number of loci
smplsize=200 # number of individuals to simulate
Mu=5e-4
#Expansion
log10thetaOBS <- -0.398
log10thetaancOBS <- 1.3
log10tauOBS <- 0.097


if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("./bug/dump",deb, "__grille=",gr, "__nloc=",n_loc, "__smpsize=", smplsize,"__log10theta=",log10thetaOBS, "__log10thetaanc=", log10thetaancOBS, "__log10tau=", log10tauOBS, sep=""), TRUE)))
}
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
  } else {
    os <- .Platform$OS.type
    if (grepl("darwin", R.version$os))
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

gr_theta <- c(-3, 4)
gr_thetaanc <- c(-2, 3.5)
gr_tau <- c(-3,3)
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


sobs <- IBDSim_wrapper(log10theta=log10thetaOBS,log10thetaanc =log10thetaancOBS,log10tau=log10tauOBS, mu = Mu, 
                       sampleSize = smplsize, nloc=n_loc,execName=IBDSimExec)
sobs

message("OBS: Nact=",10^log10thetaOBS/Mu,
        "; T=",10^log10tauOBS*10^log10thetaOBS/Mu,
        "; Nanc=",10^log10thetaancOBS/Mu)


et <- Sys.time()
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper",par.grid=parsp, nRealizations = nR, nb_cores = nbcores)
time <- Sys.time()-et
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

# We can also get the x% CIs (here, 90% CIs):
confint_theta0 <- confint(slik,"log10theta",level=0.9)
confint_thetaanc0 <- confint(slik,"log10thetaanc",level=0.9)
confint_tau0 <- confint(slik,"log10tau",level=0.9)

# Type de simulation constraction, expansion, constant)
name  <- NULL
if (log10thetaancOBS>log10thetaOBS)
{name <- "contr"
}else if (log10thetaancOBS<log10thetaOBS){
  name <- "exp"
}else{
  name <- "cst"
}



exist <- exists("slik")
exist


# Test grille : vérification des bornes
if (log10thetaOBS <= gr_theta[1]+1 | log10thetaOBS >= gr_theta[2]-1 
    | log10thetaancOBS <= gr_thetaanc[1]+1 | log10thetaancOBS >= gr_thetaanc[2]-1 
    |log10tauOBS <= gr_tau[1]+1 | log10tauOBS >= gr_tau[2]-1) 
{
  war_gr <- "Grille mal définie"
  wg <- "gm"
}else{
  war_gr <- "Grille OK"
  wg <- "go"
}


# Classement de la table selon si slik a fonctionné ou non et sleon la grille
if (exist == "TRUE")
{
  if (war_gr == "Grille OK") 
  {
    chemin <- "./ok/"
  }else{
    chemin <- "./bug/"
  }
  
}else{    # Si slik n'existe pas, on attribue des valeur nulles ou NA à slike et 
          # aux intervales de confiance pour que le script puisse s'exécuter jusqu'au bout
  chemin <- "./bug/"
  slik <- matrix(0, nrow = 8, ncol=1)
  confint_tau <- 0
  confint_tau$lowerpar <- 0
  confint_tau$lowerpar[3] <-NA
  confint_tau$upperpar[3] <-NA
  confint_theta <- 0
  confint_theta$lowerpar <- 0
  confint_theta$lowerpar[1] <- NA  
  confint_theta$upperpar[1] <- NA
  confint_thetaanc <- 0
  confint_thetaanc$lowerpar <- 0
  confint_thetaanc$lowerpar[2] <-NA
  confint_thetaanc$upperpar[2] <-NA
}




# Itérations (par groupe de trois)
nIterations <- 3
ntimes_iter <- 2
nIterations_total <- ntimes_iter*nIterations
slik0 <- slik
out <- capture.output(summary(slik0))

#########################################################
# Fichier txt qui résume les résultats et les données
sobs_dat <- capture.output(sobs)
file_name <- paste(chemin, "test_refine_", name, wg, deb, ".txt", sep="_")
write(paste(name, deb, "\n\nlog10theta =",log10thetaOBS, " Bornes theta = [", gr_theta[1], ",", gr_theta[2], "]", 
            "\nlog10thetaanc =", log10thetaancOBS, " Bornes thetaanc =[", gr_thetaanc[1], ",", gr_thetaanc[2], "]", 
            "\nlog10tau =", log10tauOBS, " Bornes tau =[", gr_tau[1], ",", gr_tau[2], "]","\n", war_gr,
            "\n\nOBS: Nact=",10^log10thetaOBS/Mu,"; T=",10^log10tauOBS*10^log10thetaOBS/Mu,"; Nanc=",10^log10thetaancOBS/Mu,
            "\n\n", sobs_dat[1], "\n", sobs_dat[2],
            "\n\nSIM: Nact(min/max)=",paste(min(10^parsp[,1]/Mu),max(10^parsp[,1]/Mu)," "),
            "\nT(min/max)=",paste(min(10^parsp[,3]*10^parsp[,1]/Mu),max(10^parsp[,3]*10^parsp[,1]/Mu)," "),
            "\nNanc(min/max)=",paste(min(10^parsp[,2]/Mu),max(10^parsp[,2]/Mu)," "),
            "\n\ntaille grille =", gr, "\nnRealizations =", nR, "\nnloc =", n_loc, "\nsmplsize =", smplsize, "\nmu =", Mu,  
            "\nNb refine =", nIterations_total, sep=" "), file=file_name)
cat("\n\n", file=file_name, append=TRUE)
cat(out, file=file_name, sep="\n", append=TRUE)
cat("\n\n", file=file_name, append=TRUE)
#######################################################################
# Refine sauvegardés à chaque itérations (et on en sauvegarde un sur trois dans le fichier texte)
for (j in 1:ntimes_iter)
{
  for(i in 1:nIterations)
    {
    slik <- refine(slik, nb_cores=nbcores)
    }
  name_slik <- paste("slik", j, sep="")
  name_out <- paste("out", i, sep="")
  assign(name_slik, slik)
  test_assign <- assign(name_slik, slik)
  test_out <- assign(name_out, capture.output(summary(test_assign)))
  cat(test_out, file=file_name, sep= "\n", append=TRUE)
  cat("\n\n", file=file_name, sep= "\n", append=TRUE)
}
# Sauvegarder slik tous les 3 refine

confint_theta <- confint(slik,"log10theta",level=0.9)
confint_thetaanc <- confint(slik,"log10thetaanc",level=0.9)
confint_tau <- confint(slik,"log10tau",level=0.9)

write(paste("Intervalles à 90% :","\nConfint theta = [", confint_theta$lowerpar[1], ",", confint_theta$upperpar[1], "]",
      "\nConfint theta anc = [", confint_thetaanc$lowerpar[2], ",", confint_thetaanc$upperpar[2], "]",
      "\nConfint tau = [", confint_tau$lowerpar[3], ",", confint_tau$upperpar[3], "]",sep=" "), file=file_name, append = TRUE)


########################################################"

# Calcul durée du script
fin <- Sys.time()
duree <- fin - deb
write(paste("\n", capture.output(duree), "\n", deb), sep="", file=file_name, append=TRUE)

# Conversion de la date deb de manière à remplacer les espaces par des "_" (pour ouvrir la base plus tard sans problème)
deb <- gsub(" ", "_", deb)

# Sauvegarde données Rdata
save.image(paste (chemin, "test_",
                  name,"_", wg,"_", deb, "_log10theta_=_",log10thetaOBS, "_log10thetaanc_=_", log10thetaancOBS, "_log10tau_=_", log10tauOBS, 
                  "_taille_grille_=_", gr,"_nRealizations_=_", nR, "_nloc_=_", n_loc, "_smplsize_=_", smplsize, "_mu_=_", Mu, 
                  ".RData", sep=""))


# Sauvegarde des graphiques (en cours)
pdf(file = paste (chemin, name, wg, deb, "plot.pdf", sep="_"))
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
dev.off()

dim(parsp)
h.simuls <- NULL
m.simuls <- NULL
k.simuls <- NULL
vk.simuls <- NULL
vh.simuls <- NULL
f.simuls <- NULL
for (i in 1:dim(parsp)[1]) 
{
  h.simuls <- c(h.simuls, simuls[[i]][,1])
  vk.simuls <- c(vk.simuls, simuls[[i]][,2])
  m.simuls <- c(m.simuls, simuls[[i]][,3]) 
  k.simuls <- c(k.simuls, simuls[[i]][,4]) 
  vh.simuls <- c(vh.simuls, simuls[[i]][,5]) 
  f.simuls <- c(f.simuls, simuls[[i]][,6]) 
}
hist(h.simuls, main="Histogramme des valeurs de H")
abline(v=sobs[1], col="red") # valeur observée
abline(v=median(h.simuls), col="blue") # mediane simulations
abline(v=mean(h.simuls), col="green") # mediane simulations
legend("topleft", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)

hist(vh.simuls)
abline(v=sobs[5], col="red")
abline(v=median(vh.simuls), col="blue") # mediane simulations
abline(v=mean(vh.simuls), col="green") # mediane simulations
legend("topright", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)

hist(k.simuls, main="Histogramme des valeurs de K")
abline(v=sobs[4], col="red") # valeur observée
abline(v=median(k.simuls), col="blue") # mediane simulations
abline(v=mean(k.simuls), col="green") # mediane simulations
legend("topright", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)

hist(vk.simuls)
abline(v=sobs[2], col="red")
abline(v=median(vk.simuls), col="blue") # mediane simulations
abline(v=mean(vk.simuls), col="green") # mediane simulations
legend("topright", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)

hist(m.simuls)
abline(v=sobs[3], col="red")
abline(v=median(m.simuls), col="blue") # mediane simulations
abline(v=mean(m.simuls), col="green") # mediane simulations
legend("topright", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)

hist(f.simuls)
abline(v=sobs[6], col="red")
abline(v=median(f.simuls), col="blue") # mediane simulations
abline(v=mean(f.simuls), col="green") # moyenne simulations
legend("topleft", legend=c("Valeur observée", "Médiane des simulations","Moyenne des simulations"),
       col=c("red", "blue", "green"), lty=1)