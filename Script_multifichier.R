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


##################  ANALYSES SUR LE CLUSTER   ##################
bug_simu <- 0

args = commandArgs(trailingOnly=TRUE)
nbcores <- args[1]
print(paste("nbcores=",nbcores))
# création dossiers
# for (j in 1:3) 
#   {
#   
#   }
name_folder <- paste("Cas1", sep="_")
dir.create(name_folder)
setwd(dir = name_folder)

# lancement jobs
for(i in 1:100)
  {
  file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim","../args.R"), "./")
  system(paste("qsub -q workq -pe parallel_smp ", nbcores, " -b y /usr/bin/Rscript ./args.R ", nbcores,
               sep=""))
}

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
res_tab <- t(c(0,0,0))
res_g <- t(c(0,0,0))
for (j in 1:100) 
  {
  setwd(paste("./Ana",j, sep=""))
  res <- read.table("Resultats.txt", sep=" ")
  res <- res[,-4]
  res_g <- cbind(res_g, res[1])
  res_tab <- cbind(res_tab,res)
  setwd("../")
  }
res_tab <- res_tab[,-c(1,2,3)]
res_g <- res_g[,-c(1,2,3)]
bias_g <- mean(as.numeric(res_g[1,])) - g_obs
  
##################  ANALYSE DES RESULTATS   ##################  

