rm(list = ls())

library(Infusion)
library(caret)
library(doSNOW)

#setwd(dir = "/work/cvernier/Infusion")
setwd(dir="/home/vernierc/Documents/GitCamille/Version locale/")

deb <- Sys.time()

source("IBD_simulation.R") 


########################################################"
if (interactive()) {options(error=recover)} else {
  options(echo = FALSE)
  options(error = quote(dump.frames(paste("./bug/dump",deb, "__grille=",gr, "__nloc=",n_loc, sep=""), TRUE)))
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


if(get_os()=="linux") {
  IBDSimExec<-"../IBDSim"
}  else if(get_os()=="osx") {
  IBDSimExec<-"../IBDSimMac"
} else {IBDSimExec<-"../ibdsimv2.0_Win7.exe"}

########################################################
nbcores <- 7
gr <- 300
nR <- 1000
n_loc=100 # number of loci
Mu=5e-4
latt <- c(100,100)
sample <- c(10,10)
min <- c(45,45)
gshape <- 0.25
em_rate <- 0.45


gr_gshape <- c(-1, 1)
gr_emigrate <- c(0, 3)

parsp <- init_grid(lower=c(gr_g=gr_gshape[1],gr_em=gr_emigrate[1]),
                   upper=c(gr_g=gr_gshape[2],gr_em=gr_emigrate[2]),
                   nUnique=gr)


sobs <- IBDSim_wrapper_IBD(lattice=latt,samp=sample,min_sample=min,nloc=n_loc, 
                            mu=Mu, g_shape = gshape, m = em_rate, execName="../IBDSim")
sobs

Infusion.options(nRealizations=c(as_one=500))
et <- Sys.time()
simuls <- add_simulation(NULL,Simulate="IBDSim_wrapper_IBD",par.grid=parsp, nRealizations = nR, nb_cores = 7, env=list2env(list(IBDSim_wrapper_IBD=IBDSim_wrapper_IBD)))
time <- Sys.time()-et
time