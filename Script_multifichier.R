######################################################
############    SCRIPT R MULTIFICHIERS    ############
######################################################

# création d'un dossier "ANALYSES" sur le cluster
# copie de IBD_sim et des scripts R de simulation et d'inférence
# script d'inference sans valeur de paramètres: on les remplit dans ce script?
# dossier [cas1] avec 100 sous dossiers [Ana_i]

# setwd(dir="/Users/verni/Documents/M2_Biostat/Stages/Scripts_R/")
# setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

##################  PARAMETRES OBSERVES   ##################  
g_obs <- 0.575
m_obs <- 0.25
h_obs <- 70

##################  SOUMISSIONSUR LE CLUSTER   ##################
bug_simu <- 0

args = commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

name_folder <- paste("Cas2_avec_Q", sep="_")
dir.create(name_folder)
setwd(dir = name_folder)
file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")


# lancement jobs
for(i in 1:100)
  {
  name_folder2 <- paste("Ana",i, sep="_")
  dir.create(name_folder2)
  setwd(dir = name_folder2)
  file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")
  system(paste("qsub -q workq -pe parallel_smp ", nbcores, " -b y \"R CMD BATCH --no-save --no-restore '--args nbcores=", nbcores,"' IBD_inference_Camille.R testibd",i,".Rout\"", sep=""))
  setwd(dir = "../")
}


  

##################  VERIFICATION DES FICHIERS ET RECUPERATION DES RESULTATS  ##################  
# setwd(dir = "./Cas2_avec_Q")
# bug_simu <- 0
# 
# 
# deb <- Sys.time()
# res_tab<-c("g","m","habitatsize")
# ci_tab_g<-c("lower_g","uper_g")
# ci_tab_m<-c("lower_m","upper_m")
# ci_tab_h<-c("lower_h","upper_h")
# 
# for (i in 1:100)
# {
#   name_folder <-  paste("./Ana", i, sep="_") 
#   setwd(dir = name_folder)
#   if (file.exists("Resultats.txt")==FALSE)
#   {
#     bug_simu <- bug_simu + 1
#   }else{
#     if(file.info("Resultats.txt")$size==0)
#     {
#       bug_simu <- bug_simu + 1
#     }
#   }
#   ##################    ANALYSE  ##################
#   res <- read.table("Resultats.txt", sep=" ")
#   res_tab <- rbind(res_tab, res[1:3,])
#   ci_tab_g <- rbind(ci_tab_g, res[4:5,])
#   ci_tab_m <- rbind(ci_tab_m, res[6:7,])
#   ci_tab_h <- rbind(ci_tab_h, res[8:9,])
#   setwd(dir = "../")
# }
# 
# 
# biais_g <- mean(as.numeric(res_tab[2:101,1])) - g_obs
# biais_m <- mean(as.numeric(res_tab[2:101,2])) - m_obs
# biais_h <- mean(as.numeric(res_tab[2:101,3])) - h_hobs
# 
# biais_g
# biais_m
# biais_h
# 
# MSE_g <- mean((as.numeric(res_tab[2:101,1]) - g_obs)^2)
# MSE_m <- mean((as.numeric(res_tab[2:101,2]) - m_obs)^2)
# MSE_h <- mean((as.numeric(res_tab[2:101,3]) - h_obs)^2)
# 
# nb_NA_ci_g <- sum(is.na(ci_tab_g)==TRUE)
# nb_NA_ci_m <- sum(is.na(ci_tab_m)==TRUE)
# nb_NA_ci_h <- sum(is.na(ci_tab_h)==TRUE)
# 
# # cover_ci_g <- sum(as.numeric(ci_tab_g[2:101,1])<=g_obs & as.numeric(ci_tab_g[2:101,2])>=g_obs & is.na(as.numeric(ci_tab_g[2:101,]))==FALSE)
# # cover_ci_m <- sum(is.na(ci_tab_m)==TRUE)
# # cover_ci_h <- sum(is.na(ci_tab_h)==TRUE)
# 
# save.image(file=paste(gsub(" ", "_", deb),".Rdata", sep=""))
# 



















#  system(paste("qsub -q workq -o out",i,".txt -pe parallel_smp ", nbcores, 
#            " -b y Rscript ./IBD_inference_Camille.R ", nbcores,sep=""))
# R CMD BATCH --no-save --no-restore '--args a=1 b=c(2,5,6)'
# test.R test.out &

# system(paste("qsub -q workq -pe parallel_smp ", nbcores, " -b y /usr/bin/Rscript args.R args",
#              i,".out ", nbcores, sep=""))
# paste("qsub -q workq -pe parallel_smp ", nbcores, "R --vanilla --args ", nbcores,
#       " --redirOut < IBD_inference_Camille.R >  ana",i,".out", sep="")
# "R --vanilla --args --redirOut <  scriptToRun.R >  scriptToRun.out"
# paste("qsub -q workq -pe parallel_smp", nbcores, "-b y /usr/bin/Rscript IBD_inference_Camille.R", nbcores)
# -b y /usr/bin/Rscript scriptToRun.R param1 param2
# appeler directement R sans passer par du bash

# copie de IBD_inference/IBD_simulation #ou bien appel dans un fichier externe
# appel d'IBDSim dans un unique fichier en dehors?
# faire un qsub sur le cluster
# unlink?

# name_file <-  paste("Ana", i, ".txt", sep="_")
# file.create(name_file)
##################  VERIFICATION DES FICHIERS   ##################  

# for (i in 1:100)
# {
#   name_file <-  paste("Ana", i, ".txt", sep="_")
#   if (file.exists(name_file)==FALSE)
#   {
#     bug_simu <- bug_simu + 1
#   }
#   else
#   {
#     if(file.info(name_file)$size==0)
#     {
#       bug_simu <- bug_simu + 1
#     }
#   }
# }

##################  RECUPERATION DES RESULTATS   ##################  
# res_tab <- t(c(0,0,0))
# res_g <- t(c(0,0,0))
# for (j in 1:100) 
#   {
#   setwd(paste("./Ana",j, sep=""))
#   res <- read.table("Resultats.txt", sep=" ")
#   res <- res[,-4]
#   res_g <- cbind(res_g, res[1])
#   res_tab <- cbind(res_tab,res)
#   setwd("../")
#   }
# res_tab <- res_tab[,-c(1,2,3)]
# res_g <- res_g[,-c(1,2,3)]
# bias_g <- mean(as.numeric(res_g[1,])) - g_obs

