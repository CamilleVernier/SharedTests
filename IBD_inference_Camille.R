rm(list = ls())

# if (interactive()) {options(error=recover)} else {
#   options(echo = FALSE)
#   options(error = quote(dump.frames(paste("dump",Sys.time(), "__grille=",gr, "__nloc=",n_loc, sep=""), TRUE)))
# }
options(warn=2)
options(error = quote(dump.frames(paste("dump",Sys.time(), "__grille=",gr, "__nloc=",n_loc, sep=""), TRUE)))




############################### FCT GET OS ############################### 

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




######################### LIBRARIES, INFUSION OPTION, IBDSIM EXEC, SOURCE ######################### 

IBDSimExec<-"../IBDSim"
#IBDSimExec<-"/home/vernierc/Documents/GitCamille/SharedTests/IBDSim"

library(Infusion)
library(caret)
library(doSNOW)


Infusion.options(nb_cores = c(param=nbcores)) #parallelisation


# setwd(dir = "/work/cvernier/Infusion")
# setwd(dir="/Users/raph/Downloads/++Ajeter/Camille/")
setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

source("IBD_simulation_Camille.R") 



############################### PARAMETRES OBS ############################## 

deb <- Sys.time()

nbcores <- 6
gr <- 500
nR <- 1

g_obs <- 0.575
m_obs <- 0.25

######################### ECHELLE NON LOG ######################### 

gshape <- 0.575
em_rate <- 0.25


gr_gshape <- c(0, 1)
gr_emigrate <- c(0, 1)
# habitatSize_bounds <-c(16,400)

parsp <- init_grid(lower=c(g_shape=gr_gshape[1],m=gr_emigrate[1]),
                   upper=c(g_shape=gr_gshape[2],m=gr_emigrate[2]),
                   nUnique=gr)

parsp2 <- unique(parsp)

parsp <- init_grid(lower=c(g_shape=gr_gshape[1],m=gr_emigrate[1],habitatSize=habitatSize_bounds[1]),
                   upper=c(g_shape=gr_gshape[2],m=gr_emigrate[2],habitatSize=habitatSize_bounds[2]),
                   nUnique=gr)
parsp2 <- unique(parsp)

sobs <- IBDSim_wrapper_IBD(habitatSize = habitatSizeObs, g_shape = gshape, m = em_rate, execName=IBDSimExec)
sobs
############################# ECHELLE LOG10 ############################# 
#  gshape <- 0.25
#  em_rate <- 0.45
# 
# # gr_gshape <- c(0, 1)
# # gr_emigrate <- c(0, 1)
# 
# log10gshape <- log10(0.25)
# log10em_rate <- log10(0.45)
# 
# log10gr_gshape <- c(-2, 0)
# log10gr_emigrate <- c(-2, 0)
# 
# 
# 
# parsp <- init_grid(lower=c(log10gshape=log10gr_gshape[1],log10m=log10gr_emigrate[1]),
#                    upper=c(log10gshape=log10gr_gshape[2],log10m=log10gr_emigrate[2]),
#                    nUnique=gr)
# 
# parsp2 <- unique(parsp)
##############################################"

# En cas de bug, retour au dossier
# setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

sobs <- IBDSim_wrapper_IBD(execName=IBDSimExec)
sobs

simtable <- add_reftable(Simulate="IBDSim_wrapper_IBD",par.grid=parsp2, nb_cores=c(param=nbcores))

dens <- infer_SLik_joint(simtable,stat.obs=sobs)

slik_j <- MSL(dens)
plot(slik_j)

slik_j <- refine(slik_j, maxit=5, nb_cores=c(param=nbcores))
slik_j2 <- refine(slik_j, nb_cores=c(param=nbcores))

Rmixmod::plotCluster(slik_j$jointdens,slik_j$logLs,variable1 = "g_shape",variable2="m")

Infusion.options(nRealizations=c(as_one=nR))
#Infusion.options(nRealizations=nR)



###################################### PROJECTIONS ###################################### 

allstats <- c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept")

gproj <- project("g_shape",stats=allstats,data=simtable,method="neuralNet")
mproj <- project("emmig_rate",stats=allstats,data=simtable,method="neuralNet")
habsizeproj <- project("habitat_size",stats=allstats,data=simtable,method="neuralNet")

# We apply projections on simulated statistics:
corrSimuls <- project(simtable,projectors=list("g"=gproj,"m"=mproj,"habsize"=habsizeproj))
corrSobs <- project(sobs,projectors=list("g"=gproj,"m"=mproj,"habsize"=habsizeproj))


###################################### ADD SIMULATION ###################################### 
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