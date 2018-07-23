####################################
##            ABC RF              ## 
####################################

rm(list = ls())

# setwd(dir = "/work/cvernier/ABCRF")
# setwd(dir="/Users/raph/Downloads/++Ajeter/Camille/")
# setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

source("IBD_simulation_Camille.R") 

if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("dump",gsub(" ", "_", Sys.time()), sep=""), TRUE)))
}
# options(warn=0)
# options(error = quote(dump.frames(paste("/home/vernierc/Documents/GitCamille/SharedTests/Bugs/dump"
#                                        ,gsub(" ", "_", Sys.time()), sep=""), TRUE)))


######################### LIBRARIES ######################### 

library(abcrf)
library(KScorrect) # pour les log uniformes
library(doParallel) 


############################### PARAMETRES ############################## 
#txtStart("./temp50.txt") #capture input and output in a txt file
deb <- Sys.time()
print(deb)

IBDSimExec<-"../IBDSim"
#IBDSimExec<-"/home/vernierc/Documents/GitCamille/SharedTests/IBDSim"

###  Arguments nbcores (à tester)
# args = commandArgs(trailingOnly=TRUE)
# for(i in 1:length(args)){
#   eval(parse(text=args[[i]]))
# }
# print(paste(nbcores))

nbcores <- 4
gr <- 200
nR <- 1

latt=c(40,40)
sample=c(15,15)
minsample=c(15,15)

d=1
n_tree = 500
nb_sim=100 #table de ref
n_sim=1 #IBDSim
n_loc=20
Mu=5e-4

dist_max=20

g_obs <- 0.575
m_obs <- 0.25
log10_g_obs <- log10(g_obs)
log10_m_obs <- log10(m_obs)
log_10=TRUE
habitatsize_obs=70

IBDSimExec<-"../IBDSim"

######################### PRIORS ######################### 

# Prior (non informatif) pour les paramètres
sim_g <- rlunif(n = nb_sim, min = 0.01, max = 2)
sim_m <- rlunif(nb_sim, 0.01, 2)
sim_habitatsize <- runif(nb_sim, 1, 200)

######################### STATS RESUMEES #########################

sobs <- IBDSim_wrapper_IBD(g=g_obs, m=m_obs, habitatSize=habitatsize_obs, mu=Mu, nloc=n_loc,
                           lattice=latt, samp=sample, min_sample=minsample, D=d, nsim=n_sim,
                           dist_max=20,log10=log_10, execName=IBDSimExec)
# sobs <- IBDSim_wrapper_IBD(g=log10_g_obs, m=log10_m_obs, habitatSize=habitatsize_obs, mu=Mu,
#                            log10=log_10, nloc=n_loc, lattice=latt, samp=sample, min_sample=minsample,
#                            D=d, nsim=n_sim,dist_max=20,execName=IBDSimExec)
sobs


######################### TABLE DE REFERENCE #########################

# Table de référence avec les valeurs des paramètres et les stats résumées
ncores <- nbcores 

cl <- makeCluster(ncores) # Tu 'initialises' ton cluster
registerDoParallel(cl)   # Obligatoire avant de faire les choses en parall?lisation
sim_concat <- foreach(i=1:nb_sim, .combine="rbind") %dopar% 
{
  return(iter_ibdsim(i, sim_g, sim_m, sim_habitatsize, execName=IBDSimExec)) 
}
stopCluster(cl)    # Obligatoire apr?s la parall?lisation
colnames(sim_concat)[1:3] <- c("g", "m", "habitat_size")

sim_concat
sim_concat2 <- sim_concat
rownames(sim_concat2) <- NULL


######################### DATA FRAME #########################

ref_table<- data.frame(sim_concat2)
# ref_table$param <- as.vector(ref_table[,1:3])
# ref_table$sumsta <- as.vector(ref_table[,4:9])
# ref_table2 <- ref_table[,-c(1:9)]
# data2 <- data.frame(ref_table)
data2 <- data.frame(ref_table$g, ref_table[,-c(1:3)])


######################### RANDOM FOREST #########################
# Construction de Ntree arbres CART qui prédissent les valeurs des paramètes d'après S(x)

model_RF <- regAbcrf(formula = data2$ref_table.g~.,
                     data    = data2,
                     ntree   = n_tree,
                     paral   = T)

plot(model_RF,training=data2)

predictOOB(model_RF, data2)
err.regAbcrf(model_RF, data2)
sobs_df <- t(as.data.frame(sobs))
# pred <- predict(model_RF, obs = sobs_df, training=data2, paral = TRUE) # bug, à revoir
# Nb arbres: 500

fin <- Sys.time()
duree <- fin - deb
