######################################################
############    SCRIPT R MULTIFICHIERS    ############
######################################################

# création d'un dossier "ANALYSES" sur le cluster
# copie de IBD_sim et des scripts R de simulation et d'inférence
# script d'inference sans valeur de paramètres: on les remplit dans ce script?
# dossier [cas1] avec 100 sous dossiers [Ana_i]

# setwd(dir="/Users/verni/Documents/M2_Biostat/Stages/Scripts_R/")

# library(rslurm)

##################  PARAMETRES OBSERVES   ##################  


##################  ANALYSES SUR LE CLUSTER   ##################
bug_simu <- 0

# création dossiers
name_folder <- paste("Cas2", sep="_")
dir.create(name_folder)
setwd(dir = name_folder)

# lancement jobs
for(i in 1:2)
  {
  # name_file <-  paste("Ana", i, ".txt", sep="_")
  # file.create(name_file)
  file.copy(c("../multi_test.sh","../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")
  system("qsub -q workq -pe parallel_smp 10 multi_test.sh") 
  }



# copie de IBD_inference/IBD_simulation #ou bien appel dans un fichier externe
# appel d'IBDSim dans un unique fichier en dehors?
# faire un qsub sur le cluster
# unlink?


##################  VERIFICATION DES FICHIERS   ##################  

# for (i in 1:2)
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

##################  ANALYSE DES RESULTATS   ##################  





