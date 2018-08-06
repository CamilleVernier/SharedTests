##################  VERIFICATION DES FICHIERS ET RECUPERATION DES RESULTATS  ##################  
setwd(dir = "./Cas2_avec_Q")
bug_simu <- 0
g_obs <- 0.575
m_obs <- 0.25
h_obs <- 70

deb <- Sys.time()
res_tab<-c("g","m","habitatsize")
ci_tab_g<-c("lower_g","uper_g")
ci_tab_m<-c("lower_m","upper_m")
ci_tab_h<-c("lower_h","upper_h")

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
  setwd(dir = "../")
}


biais_g <- mean(as.numeric(res_tab[2:101,1])) - g_obs
biais_m <- mean(as.numeric(res_tab[2:101,2])) - m_obs
biais_h <- mean(as.numeric(res_tab[2:101,3])) - h_obs

biais_g
biais_m
biais_h

MSE_g <- mean((as.numeric(res_tab[2:101,1]) - g_obs)^2)
MSE_m <- mean((as.numeric(res_tab[2:101,2]) - m_obs)^2)
MSE_h <- mean((as.numeric(res_tab[2:101,3]) - h_obs)^2)

nb_NA_ci_g <- sum(is.na(ci_tab_g)==TRUE)
nb_NA_ci_m <- sum(is.na(ci_tab_m)==TRUE)
nb_NA_ci_h <- sum(is.na(ci_tab_h)==TRUE)

# cover_ci_g <- sum(as.numeric(ci_tab_g[2:101,1])<=g_obs & as.numeric(ci_tab_g[2:101,2])>=g_obs & is.na(as.numeric(ci_tab_g[2:101,]))==FALSE)
# cover_ci_m <- sum(is.na(ci_tab_m)==TRUE)
# cover_ci_h <- sum(is.na(ci_tab_h)==TRUE)

save.image(file=paste(gsub(" ", "_", deb),".Rdata", sep=""))


