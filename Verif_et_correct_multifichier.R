
setwd(dir = "./Cas_avec_Q_2018-08-08_18:51:28")
bug_simu <- 0

args = commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")

for (i in 1:100)
{
  name_folder <-  paste("./Ana", i, sep="_") 
  setwd(dir = name_folder)
  if (file.exists("Resultats.txt")==FALSE)
  {
    bug_simu <- bug_simu + 1
    print(getwd())
    #file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")
    system(paste("qsub -q workq -pe parallel_smp ", nbcores, " -b y \"R CMD BATCH --no-save --no-restore '--args nbcores=", nbcores,"' IBD_inference_Camille.R testibd_correct667",i,".Rout\"", sep=""))
  }else{
    if(file.info("Resultats.txt")$size==0)
    {
      bug_simu <- bug_simu + 1
      print(getwd())
      #file.copy(c("../IBD_inference_Camille.R", "../IBD_simulation_Camille.R","../IBDSim"), "./")
      system(paste("qsub -q workq -pe parallel_smp ", nbcores, " -b y \"R CMD BATCH --no-save --no-restore '--args nbcores=", nbcores,"' IBD_inference_Camille.R testibd_correct667",i,".Rout\"", sep=""))
    }
  }
  setwd("../")
}


