rm(list = ls())

library(Infusion)
library(caret)
library(doSNOW)

#setwd(dir = "/work/cvernier/Infusion")
setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

deb <- Sys.time()

source("IBD_simulation_FR.R") 


########################################################"
if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("./bug/dump",Sys.time(), "__grille=",gr, "__nloc=",n_loc, sep=""), TRUE)))
}
options(warn=2)

# get_os <- function(){
#   sysinf <- Sys.info()
#   if (!is.null(sysinf)){
#     os <- sysinf['sysname']
#     if (os == 'Darwin')
#       os <- "osx"
#   } else {
#     os <- .Platform$OS.type
#     if (grepl("darwin", R.version$os))
#       os <- "osx"
#     if (grepl("linux-gnu", R.version$os))
#       os <- "linux"
#   }
#   tolower(os)
# }
# if(get_os()=="linux") {
#   IBDSimExec<-"../IBDSim"
# }  else if(get_os()=="osx") {
#   IBDSimExec<-"../IBDSimMac"
# } else {IBDSimExec<-"../ibdsimv2.0_Win7.exe"}

IBDSimExec<-"../IBDSim"
#IBDSimExec<-"/home/vernierc/Documents/GitCamille/SharedTests/IBDSim"

########################################################
nbcores <- 6
gr <- 500
nR <- 1
## attention a changer ces valeures de parametres comme valeurs par défault dans le wrapper
# n_loc=20 # number of loci
# Mu=5e-3
# latt <- c(20,20)
# sample <- c(10,10)
# min <- c(5,5)

############################ ECHELLE NON LOG
# gshape <- 0.25
# em_rate <- 0.45
# 
# 
# gr_gshape <- c(0, 1)
# gr_emigrate <- c(0, 1)
# 
# parsp <- init_grid(lower=c(g_shape=gr_gshape[1],m=gr_emigrate[1]),
#                    upper=c(g_shape=gr_gshape[2],m=gr_emigrate[2]),
#                    nUnique=gr)
# 
# parsp2 <- unique(parsp)

############################# ECHELLE LOG10
log10_g_OBS <- log10(0.25)
log10_m_OBS <- log10(0.45)
log10_g_OBS
log10_m_OBS

gr_gshape <- c(-2, 0)
gr_emigrate <- c(-2, 0)

parsp <- init_grid(lower=c(g_shape=gr_gshape[1],m=gr_emigrate[1]),
                   upper=c(g_shape=gr_gshape[2],m=gr_emigrate[2]),
                   nUnique=gr)

parsp2 <- unique(parsp)

sobs <- IBDSim_wrapper_IBD(g_shape = gshape, m = em_rate, execName=IBDSimExec)
sobs

simtable <- add_reftable(Simulate="IBDSim_wrapper_IBD",par.grid=parsp2, nb_cores=c(param=nbcores))

dens <- infer_SLik_joint(simtable,stat.obs=sobs)

slik_j <- MSL(dens)
plot(slik_j)

slik_j <- refine(slik_j, maxit=5, nb_cores=nbcores)
slik_j2 <- refine(slik_j, nb_cores=c(param=nbcores))

Rmixmod::plotCluster(slik_j$jointdens,slik_j$logLs,variable1 = "g_shape",variable2="m")

Infusion.options(nRealizations=c(as_one=nR))
#Infusion.options(nRealizations=nR)

et <- Sys.time()
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper_IBD", par.grid=parsp,nRealizations=c(as_one=nR), nb_cores = nbcores, env=list2env(list(IBDSim_wrapper_IBD=IBDSim_wrapper_IBD)))
time <- Sys.time()-et
time







# for (i in 1:dim(parsp)[1]) {
#   setwd(dir="/Users/raph/Downloads/++Ajeter/Camille/")
#   message(paste("iter=",i,"/",dim(parsp)[1]))
#   print(parsp[i,])
#   for (j in 1:nR){
#     message(paste("Rep=",j,"/",nR))
#     res <- IBDSim_wrapper_IBD(lattice=latt,samp=sample,min_sample=min,nloc=n_loc,
#                               mu=Mu, g_shape = parsp[i,1], m = parsp[i,2], execName=IBDSimExec,nsim=1)
#     print(res)
#   }
# }