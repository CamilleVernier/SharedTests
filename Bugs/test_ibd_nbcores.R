rm(list = ls())


source("IBD_simulation_Camille.R") 

library(Infusion)
library(caret)
library(doSNOW)


deb <- Sys.time()
print(deb)

IBDSimExec<-"../IBDSim"

nbcores <- 2


gr <- 200
nR <- 1

latt=c(40,40)
sample=c(15,15)
minsample=c(15,15)

d=1
n_sim=1
n_loc=20
Mu=5e-4

dist_max=20

g_obs <- 0.575
m_obs <- 0.25
log10_g_obs <- log10(g_obs)
log10_m_obs <- log10(m_obs)
log_10=FALSE
habitatsize_obs=70


g_grille <- c(0, 1)
m_grille <- c(0, 1)
log10_g_grille <- c(-2, 0)
log10_m_grille <- c(-2, 0)
#habitatsize_grille <-c(2,100)
habitatsize_grille <-c(16,200)


Infusion.options(nb_cores = c(param=nbcores)) #parallelisation
Infusion.options(nRealizations=c(as_one=nR))

parsp <- init_grid(lower=c(g=g_grille[1],m=m_grille[1], habitatSize=habitatsize_grille[1]),
                   upper=c(g=g_grille[2],m=m_grille[2], habitatSize=habitatsize_grille[2]),
                   nUnique=gr)

# parsp <- init_grid(lower=c(g=log10_g_grille[1],m=log10_m_grille[1], habitatSize=habitatsize_grille[1]),
#                    upper=c(g=log10_g_grille[2],m=log10_m_grille[2], habitatSize=habitatsize_grille[2]),
#                    nUnique=gr)

parsp2 <- unique(parsp)


# sobs <- IBDSim_wrapper_IBD(habitatsize_obs = habitatsize_obsObs, g = g_obs, m = m_obs, execName=IBDSimExec)

sobs <- IBDSim_wrapper_IBD(g=g_obs, m=m_obs, habitatSize=habitatsize_obs, mu=Mu, nloc=n_loc,
                           lattice=latt, samp=sample, min_sample=minsample, D=d, nsim=n_sim,
                           dist_max=20,log10=log_10, execName=IBDSimExec)
# sobs <- IBDSim_wrapper_IBD(g=log10_g_obs, m=log10_m_obs, habitatSize=habitatsize_obs, mu=Mu,
#                            log10=log_10, nloc=n_loc, lattice=latt, samp=sample, min_sample=minsample,
#                            D=d, nsim=n_sim,dist_max=20,execName=IBDSimExec)
sobs



simtable <- add_reftable(Simulate="IBDSim_wrapper_IBD", nloc=n_loc, par.grid=parsp2,
                         lattice=latt, samp=sample, min_sample=minsample, D=d, log10=log_10,
                         dist_max=20,execName=IBDSimExec,
                         nb_cores=c(param=nbcores))
#simtable

print(42)
sys <- Sys.time()
save.image(file=paste(sys,".Rdata", sep=""))