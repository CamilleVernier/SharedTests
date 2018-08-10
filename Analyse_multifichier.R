##################  VERIFICATION DES FICHIERS ET RECUPERATION DES RESULTATS  ##################  
setwd(dir = "./Cas_avec_Q_grille=300_2018-08-09_16:49:21")
bug_simu <- 0
g_obs <- 0.558615
m_obs <- 0.25
h_obs <- 70

deb <- Sys.time()
res_tab<-c("g","m","habitatsize")
ci_tab_g<-c("lower_g","uper_g")
ci_tab_m<-c("lower_m","upper_m")
ci_tab_h<-c("lower_h","upper_h")
temps_tab_sim <- c("sim","10 refine", "total")

for (i in 1:100)
{
  name_folder <-  paste("./Ana", i, sep="_") 
  setwd(dir = name_folder)
  if (file.exists("Resultats.txt")==FALSE)
  {
    bug_simu <- bug_simu + 1
  }else{
    if(file.info("Resultats.txt")$size==0)
    {
      bug_simu <- bug_simu + 1
    }
  }

  ##################    ANALYSE  ##################
  res <- read.table("Resultats.txt", sep=" ")
  res_tab <- rbind(res_tab, res[1:3,])
  ci_tab_g <- rbind(ci_tab_g, res[4:5,])
  ci_tab_m <- rbind(ci_tab_m, res[6:7,])
  ci_tab_h <- rbind(ci_tab_h, res[8:9,])
  temps_tab_sim <- rbind(temps_tab_sim,c(res[10,],res[11,],res[10,]+res[11,]))
  setwd(dir = "../")
}


mean_g <- mean(as.numeric(res_tab[2:101,1]))
mean_m <- mean(as.numeric(res_tab[2:101,2]))
mean_h <- mean(as.numeric(res_tab[2:101,3]))

med_g <- median(as.numeric(res_tab[2:101,1]))
med_m <- median(as.numeric(res_tab[2:101,2]))
med_h <- median(as.numeric(res_tab[2:101,3]))


biais_g <- mean_g - g_obs
biais_m <- mean_m - m_obs
biais_h <- mean_h - h_obs

biais_g
biais_m
biais_h

biais_relatif_g <- (mean_g - g_obs)/g_obs
biais_relatif_m <- (mean_m - m_obs)/m_obs
biais_relatif_h <- (mean_h - h_obs)/h_obs


MSE_g <- mean((as.numeric(res_tab[2:101,1]) - g_obs)^2)
MSE_m <- mean((as.numeric(res_tab[2:101,2]) - m_obs)^2)
MSE_h <- mean((as.numeric(res_tab[2:101,3]) - h_obs)^2)

MSE_relatif_g <- (mean((as.numeric(res_tab[2:101,1]) - g_obs)^2))/(g_obs^2)
MSE_relatif_m <- (mean((as.numeric(res_tab[2:101,2]) - m_obs)^2))/(m_obs^2)
MSE_relatif_h <- (mean((as.numeric(res_tab[2:101,3]) - h_obs)^2))/(h_obs^2)

nb_NA_ci_g <- sum(is.na(ci_tab_g)==TRUE)
nb_NA_ci_m <- sum(is.na(ci_tab_m)==TRUE)
nb_NA_ci_h <- sum(is.na(ci_tab_h)==TRUE)

nb_NA_ci_low_g <- sum(is.na(ci_tab_g[2:101,1])==TRUE)
nb_NA_ci_low_m <- sum(is.na(ci_tab_m[2:101,1])==TRUE)
nb_NA_ci_low_h <- sum(is.na(ci_tab_h[2:101,1])==TRUE)

nb_NA_ci_up_g <- sum(is.na(ci_tab_g[2:101,2])==TRUE)
nb_NA_ci_up_m <- sum(is.na(ci_tab_m[2:101,2])==TRUE)
nb_NA_ci_up_h <- sum(is.na(ci_tab_h[2:101,2])==TRUE)

for (i in 2:101)
{
  if(is.na(ci_tab_g[i,1]==TRUE))
  {
    ci_tab_g[i,1] <- 0
  }
  if(is.na(ci_tab_m[i,1]==TRUE))
  {
    ci_tab_m[i,1] <- 0
  }
  if(is.na(ci_tab_h[i,1]==TRUE))
  {
    ci_tab_h[i,1] <- 16
  }
  if(is.na(ci_tab_g[i,2]==TRUE))
  {
    ci_tab_g[i,2] <- 1
  }
  if(is.na(ci_tab_m[i,2]==TRUE))
  {
    ci_tab_m[i,2] <- 1
  }
  if(is.na(ci_tab_h[i,2]==TRUE))
  {
    ci_tab_h[i,2] <- 200
  }
}

cover_ci_g <- sum(as.numeric(ci_tab_g[2:101,1])<=g_obs & as.numeric(ci_tab_g[2:101,2])>=g_obs)
cover_ci_m <- sum(as.numeric(ci_tab_m[2:101,1])<=m_obs & as.numeric(ci_tab_m[2:101,2])>=m_obs)
cover_ci_h <- sum(as.numeric(ci_tab_h[2:101,1])<=h_obs & as.numeric(ci_tab_h[2:101,2])>=h_obs)

tableau <- cbind(c(mean_g, med_g, biais_g, biais_relatif_g, MSE_g, MSE_relatif_g, nb_NA_ci_g,
                   nb_NA_ci_low_g, nb_NA_ci_up_g, cover_ci_g),
                 c(mean_m,med_m,biais_m,  biais_relatif_m, MSE_m, MSE_relatif_m, nb_NA_ci_m,
                   nb_NA_ci_low_m, nb_NA_ci_up_m, cover_ci_m),
                 c(mean_h, med_h,biais_h, biais_relatif_h, MSE_h, MSE_relatif_h, nb_NA_ci_h,
                   nb_NA_ci_low_h, nb_NA_ci_up_h, cover_ci_h))

colnames(tableau) <- c("g","m","h")
rownames(tableau) <- c("mean", "mediane","biais", "biais_relatif", "MSE", "MSE_relatif", "nb_NA_ci",
                       "nb_NA_ci_low", "nb_NA_ci_up", "cover_ci")

write.table(tableau,paste(gsub(" ", "_", deb),"Resume_final.txt", sep="_"))

temps_sim_mean <- mean(as.numeric(temps_tab_sim[2:101,1]))
temps_refine_mean <- mean(as.numeric(temps_tab_sim[2:101,2]))
temps_mean <- mean(as.numeric(temps_tab_sim[2:101,3]))

write.table(c(temps_sim_mean, temps_refine_mean, temps_mean),paste(gsub(" ", "_", deb),"Temps.txt", sep="_"))


save.image(file=paste(gsub(" ", "_", deb),".Rdata", sep=""))


