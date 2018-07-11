rm(list = ls())
if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("./bug/dump",deb, "__grille=",gr, "__nloc=",n_loc, sep=""), TRUE)))
  options(warn=1)
}
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

########################################################
library(Infusion)
library(caret)
library(doSNOW)

#setwd(dir = "/work/cvernier/Infusion")
#setwd(dir="/home/vernierc/Documents/GitCamille/Version locale/")
setwd(dir="/Users/raph/Downloads/++Ajeter/Camille/")
source("IBD_simulation_Raph.R") 
# if(get_os()=="linux") {
#   IBDSimExec<-"../IBDSim"
# }  else if(get_os()=="osx") {
#   IBDSimExec<-"../IBDSimMac"
# } else {IBDSimExec<-"../ibdsimv2.0_Win7.exe"}
IBDSimExec<-"/Users/raph/Downloads/++Ajeter/Camille/IBDSim"

########################################################
nbcores <- 5
gr <- 500
nR <- 1

# gshape <- 0.25
# em_rate <- 0.45
g_shape_obs <- 0.575
em_rate_obs <- 0.25
hab_size_obs <- 70

g_gshape_bounds <- c(0, 1)
em_rate_bounds <- c(0, 1)
hab_size_bounds <-c(16,400)

deb <- Sys.time()

parsp <- init_grid(lower=c(g_shape=g_gshape_bounds[1],em_rate=em_rate_bounds[1],hab_size=hab_size_bounds[1]),
                   upper=c(g_shape=g_gshape_bounds[2],em_rate=em_rate_bounds[2],hab_size=hab_size_bounds[2]),
                   nUnique=gr)
parsp2 <- unique(parsp)

sobs <- IBDSim_wrapper_IBD(hab_size = hab_size_obs, g_shape = g_shape_obs, em_rate = em_rate_obs, execName=IBDSimExec)
sobs

###############
Infusion.options(nb_cores = c(param=nbcores))
rm(simtable,slik_j,dens)
simtable <- add_reftable(Simulate="IBDSim_wrapper_IBD",execName=IBDSimExec, par.grid=parsp2, nb_cores = c(param=nbcores))

allstats <- c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept")

gproj <- project("g_shape",stats=allstats,data=simtable,method="neuralNet")
mproj <- project("em_rate",stats=allstats,data=simtable,method="neuralNet")
habsizeproj <- project("hab_size",stats=allstats,data=simtable,method="neuralNet")

# We apply projections on simulated statistics:
corrSimuls <- project(simtable,projectors=list("g_shape_proj"=gproj,"em_rate_proj"=mproj,"hab_size_proj"=habsizeproj))
corrSobs <- project(sobs,projectors=list("g_shape_proj"=gproj,"em_rate_proj"=mproj,"hab_size_proj"=habsizeproj))

#dens <- infer_SLik_joint(simtable,stat.obs=drop(sobs))
dens <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs)

slik_j <-MSL(dens)
# The slik object contains the ML estimates and 95% CIs:
summary(slik_j)
plot(slik_j)
slik_j <- refine(slik_j,maxit=6,nb_cores = c(param=nbcores))
plot(slik_j)

time <- Sys.time() - deb
time

#######
rm(simtable,slik_j,dens)
simtable <- add_reftable(Simulate="IBDSim_wrapper_IBD_Git",par.grid=parsp2, nb_cores = c(param=nbcores))

dens <- infer_SLik_joint(simtable,stat.obs=drop(sobs))

slik_j <-MSL(dens)
plot(slik_j)
slik_j <- refine(slik_j,maxit=12,nb_cores = c(param=nbcores))
plot(slik_j)

Rmixmod::plotCluster(slik_j$jointdens,slik_j$logLs,variable1="m",variable2="g_shape")
#Infusion.options(nRealizations=c(as_one=nR))
#Infusion.options(nRealizations=nR)

# et <- Sys.time()
# simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper_IBD", par.grid=parsp,nRealizations=c(as_one=nR), nb_cores = nbcores, env=list2env(list(IBDSim_wrapper_IBD=IBDSim_wrapper_IBD)))
# time <- Sys.time()-et
# time


# sobs1 <- IBDSim_wrapper_IBD(g_shape = gshape, m = em_rate, execName="/Users/raph/Downloads/++Ajeter/Camille/IBDSim")
# sobs1
# 
# sobs2 <- IBDSim_wrapper_IBD(g_shape = gshape, m = em_rate, execName="/Users/raph/Downloads/++Ajeter/Camille/IBDSimGit")
# sobs2



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