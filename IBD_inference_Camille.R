rm(list = ls())

# setwd(dir = "/work/cvernier/Infusion")
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




######################### LIBRARIES ######################### 

.libPaths("/save/cvernier/R/x86_64-pc-linux-gnu-library/3.3")
loc <- "/save/cvernier/R/x86_64-pc-linux-gnu-library/3.3"
library(Infusion, lib.loc = loc)
library(caret, lib.loc = loc)
library(doSNOW, lib.loc = loc)
#library(TeachingDemos)


############################### PARAMETRES ############################## 

#txtStart("./temp50.txt") #capture input and output in a txt file
deb <- Sys.time()
print(deb)

IBDSimExec<-"../IBDSim"
#IBDSimExec<-"/home/vernierc/Documents/GitCamille/SharedTests/IBDSim"

args = commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
print(paste(nbcores))

gr <- 500
nR <- 1

latt=c(70,70)
sample=c(15,15)
minsample=c(15,15)

d=1
n_sim=1
n_loc=20
Mu=5e-4

dist_max=20

g_obs <- 0.558615
m_obs <- 0.25
log10_g_obs <- log10(g_obs)
log10_m_obs <- log10(m_obs)
log_10=FALSE
habitatsize_obs=70

######################### INFUSION OPTIONS #########################

Infusion.options(nb_cores = c(param=nbcores)) #parallelisation
Infusion.options(nRealizations=c(as_one=nR))
#Infusion.options(nRealizations=nR)

######################### ECHELLE NON LOG #########################

g_grille <- c(0, 1)
m_grille <- c(0, 1)
log10_g_grille <- c(-2, 0)
log10_m_grille <- c(-2, 0)
#habitatsize_grille <-c(2,100)
habitatsize_grille <-c(16,200)

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
                         env=list2env(list(IBDSim_wrapper_IBD=IBDSim_wrapper_IBD)),
                         nb_cores=c(param=nbcores))

fin_tps_sim <- Sys.time()
duree_sim <-  difftime(fin_tps_sim,deb, units = c("hours"))
print(duree_sim)

##################### SANS PROJECTION ##################### 

# dens <- infer_SLik_joint(simtable,stat.obs=sobs)
# 
# slik_j <- MSL(dens)
# plot(slik_j)

###############################################################   


#   slik_j <- refine(slik_j, maxit=3, nb_cores=c(param=nbcores))
# plot(slik_j)

#Rmixmod::plotCluster(slik_j$jointdens,slik_j$logLs,variable1 = "g",variable2="m")


# setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")
# source("IBD_simulation_Camille.R")

###################################### PROJECTIONS ######################################
deb2 <- Sys.time()

qnames<-NULL
for (i in 1:sample[1]) 
{
  name <- paste("Qr", i-1, sep="")
  qnames <- c(qnames,name)
}

allstats <- c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample",
              "var_nballele", qnames, "ar_slope","ar_intercept","er_slope","er_intercept")

gproj <- project("g",stats=allstats,data=simtable,method="neuralNet")
mproj <- project("m",stats=allstats,data=simtable,method="neuralNet")
habsizeproj <- project("habitatSize",stats=allstats,data=simtable,method="neuralNet")

# We apply projections on simulated statistics:
# corrSimuls <- project(simtable,projectors=list("g"=gproj,"m"=mproj,"habsize"=habsizeproj))
# corrSobs <- project(sobs,projectors=list("g"=gproj,"m"=mproj,"habsize"=habsizeproj))
corrSimuls <- project(simtable,projectors=list("g_p"=gproj,"m_p"=mproj, "h_p"=habsizeproj))
corrSobs <- project(sobs,projectors=list("g_p"=gproj,"m_p"=mproj,"h_p"=habsizeproj))

dens_proj <- infer_SLik_joint(corrSimuls,stat.obs=corrSobs)

slik_j_proj <- MSL(dens_proj)
plot(slik_j_proj, filled=TRUE)
#slik_j <- refine(slik_j, maxit=3, nb_cores=c(param=nbcores))

#pval <- dchisq(2*(slik_j$MSL$maxlogL-predict(slik_j, newdata=c(g_obs,m_obs,habitatsize_obs))[1]), df=2)

###################################### SAUVEGARDE SLIK ######################################
#file_name_slik <- "Refine.txt"
file_name_slik_proj <- "Refine_proj.txt"
niterations <- 3
ntimes_iter <- 3
niterations_total <- ntimes_iter*niterations
#slik0 <- slik_j
slik0_proj <- slik_j_proj

#out_slik0 <- capture.output(summary(slik0))
# write(cat(out_slik0, file=file_name_slik, sep= "\n"), append=FALSE)
# write(cat("\n\n", file=file_name_slik, sep= "\n"), append=TRUE)

out_slik_proj <- capture.output(summary(slik0_proj))
# write(cat(out_slik_proj, file=file_name_slik, sep= "\n"), append=FALSE)
# write(cat("\n\n", file=file_name_slik, sep= "\n"), append=TRUE)

for (j in 1:ntimes_iter)
{
  for(i in 1:niterations)
  {
    #slik_j <- refine(slik_j, nb_cores=nbcores)
    slik_j_proj <- refine(slik_j_proj, nb_cores=nbcores)
  }
  #name_slik <- paste("slik", j, sep="")
  #name_out <- paste("out", i, sep="")
  #assign(name_slik, slik_j)
  #test_assign <- assign(name_slik, slik_j)
  #test_out <- assign(name_out, capture.output(summary(test_assign)))
  # write(cat(test_out, file=file_name_slik, sep= "\n"), append=TRUE)
  # write(cat("\n\n", file=file_name_slik, sep= "\n"), append=TRUE)

  name_slik_proj <- paste("slik_proj", j, sep="")
  name_out_proj <- paste("out_proj", i, sep="")
  assign(name_slik_proj, slik_j_proj)
  test_assign_proj <- assign(name_slik_proj, slik_j_proj)
  test_out_proj <- assign(name_out_proj, capture.output(summary(test_assign_proj)))
  write(cat(test_out_proj, file=file_name_slik_proj, sep= "\n"), append=TRUE)
  write(cat("\n\n", file=file_name_slik_proj, sep= "\n"), append=TRUE)
}

# tmp <- tempfile(pattern="Ana", tmpdir= ".",fileext=".txt")
# out_slik <- capture.output(summary(slik_j))
# cat(out_slik, file=tmp, sep="\n", append=TRUE)
# out_slik_proj <- capture.output(summary(slik_j_proj))
# cat(out_slik_proj, file=tmp, sep="\n", append=TRUE)

fin_10_refine <- Sys.time()
duree_10_refine <- difftime(fin_10_refine,deb2, units = c("hours"))
print(duree_10_refine)

deb_txt <- gsub(" ", "_", deb)

# write(cat(slik_j$MSL$MSLE,"\n", file=paste(deb,"Resultats.txt", sep="")))
# write(cat(slik_j$lower,"\n", file=paste(deb,"Resultats.txt", sep="")), append=TRUE)
# write(cat(slik_j$upper,"\n", file=paste(deb,"Resultats.txt", sep="")), append=TRUE)
write(c(slik_j_proj$MSL$MSLE,"\n",slik_j_proj$CIobject$CIs$g$interval, "\n", slik_j_proj$CIobject$CIs$m$interval, "\n", slik_j_proj$CIobject$CIs$habitatSize$interval,"\n", duree_sim,"\n", duree_10_refine), file="Resultats.txt")

save.image(file=paste(gsub(" ", "_", deb_txt),"grille=",gr,".Rdata", sep=""))




#txtStop()
###################################### ADD SIMULATION ######################################
# et <- Sys.time()
# simuls2 <- add_simulation(NULL,Simulate="IBDSim_wrapper_IBD", par.grid=parsp,nRealizations=c(as_one=nR),
#                          nb_cores = 7, env=list2env(list(IBDSim_wrapper_IBD=IBDSim_wrapper_IBD)))
# time <- Sys.time()-et
# time


############################# ECHELLE LOG10 #############################
#  g_obs <- 0.25
#  m_obs <- 0.45
#
# # g_grille <- c(0, 1)
# # m_grille <- c(0, 1)
#
# log10g_obs <- log10(0.25)
# log10m_obs <- log10(0.45)
#
# log10g_grille <- c(-2, 0)
# log10m_grille <- c(-2, 0)
#
#
#
# parsp <- init_grid(lower=c(log10g_obs=log10g_grille[1],log10m=log10m_grille[1]),
#                    upper=c(log10g_obs=log10g_grille[2],log10m=log10m_grille[2]),
#                    nUnique=gr)
#
# parsp2 <- unique(parsp)
##############################################"




# for (i in 1:dim(parsp)[1]) {
#   setwd(dir="/Users/raph/Downloads/++Ajeter/Camille/")
#   message(paste("iter=",i,"/",dim(parsp)[1]))
#   print(parsp[i,])
#   for (j in 1:nR){
#     message(paste("Rep=",j,"/",nR))
#     res <- IBDSim_wrapper_IBD(lattice=latt,samp=sample,min_sample=min,nloc=n_loc,
#                               mu=Mu, g = parsp[i,1], m = parsp[i,2], execName=IBDSimExec,nsim=1)
#     print(res)
#   }
# }