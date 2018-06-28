### This script generates the six statistics used in this study from 
# a PSC model with a given set of parameter, by calling IBDsim from R.
# IBDsim v2.0 reference: Leblois et al. 2009 MolEcol Ressources

IBDSim_wrapper<-function(log10theta,
                 log10a,
                 log10tau,
                 nsim=1,
                 nloc=40, # number of loci
                 sampleSize=100, # number of individuals to simulate
                 mu=5e-2, # mutation rate
                 execName="../IBDSim"){ # executable name
  #conversion from the scaled parameters:
  # créer dossier simul + chiffre aleatoire
  # setwd dans le dossier
 
  # name <- paste("sim",sample(1:10000, 1), sep="_")
  curDir<-getwd()
  dir1 <- tempfile(pattern = "sim", tmpdir =curDir )
  dir.create(dir1)
  #message("Directory created:",dir1)
  setwd(dir = dir1)
  #file.copy("../IBDSim", dir1)
  file.copy("/home/vernierc/Documents/GitCamille/SharedTests/IbdSettings.txt", dir1)
  #file.copy("/Users/raph/Documents/Taf/EtudiantsPostdocVisiteurs+Stages+Theses/_CamilleVernier/GitHub/SharedTests/IbdSettings.txt", dir1) #pour Raph
  
  
  #message("Files in the currect directory: ", paste(list.files(path = ".", all.files = TRUE,recursive = TRUE, include.dirs = TRUE, no.. = TRUE)," "))
  

  a <- 10^log10a
  N0<-(10^log10theta)/mu
  t<-(10^log10tau)*N0
  if (t<1) t=1
  Seed <- 1234567 + sample(1:10000, 1)

  # we write the input file for IBDsim:
  write.table(paste("%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%
                    Run_Number=",nsim,"
                    Migraine_Settings=F
                    Ploidy=Diploid
                    Random_Seeds=",Seed,"

                    %%%%% MARKERS PARAMETER S%%%%%%%%%%%%%%%
                    Locus_Number=",nloc,"
                    Mutation_Rate=",format(mu, scientific = FALSE),"
                    Mutation_Model=GSM
                    Allelic_Upper_Bound=200
                    %%%%%%%% DEMOGRAPHIC OPTIONS %%%%%%%%%%%%%
                    %% LATTICE
                    Lattice_SizeX=1
                    Lattice_SizeY=1
                    Ind_Per_Pop=",N0,"
                    %% SAMPLE
                    Sample_SizeX=1
                    Sample_SizeY=1
                    Min_Sample_CoordinateX=1
                    Min_Sample_CoordinateY=1
                    Ind_Per_Pop_Sampled=",sampleSize,"
                    %% DISPERSAL
                    Dispersal_Distribution=b
                    TotalEmigrationRate=0
                    Continuous_Deme_Size_Variation=Exponential

                    NewDemographicPhaseAt=",t,"
                    Continuous_Deme_Size_Variation=None
                    Ind_Per_Pop=",format(N0*a, scientific = FALSE),"
                    Genepop=false
                    DiagnosticTables=[Hexp,Fis,Allelic_Variance,Iterative_Statistics]",sep=""),
              file="IbdSettings.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


  # We call IBDsim through the system command,
  # it will read the "IbdSettings.txt" file and run the simulations:
  system(execName, ignore.stdout = TRUE)

  # We then load the simulated summary statistics and only keep the 6 used in this study:
  sumstats<-read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1)
  
  dat<-sumstats[,(dim(sumstats)[2]-4):dim(sumstats)[2]][c(1,3,5,2)] # H, K, M and varK
  colnames(dat)<-c("H","M","K","VarK") # H = Hexp ?
  dat$varH<-apply(sumstats[,(nloc+1):(2*nloc)],1,var) # variance of H ##modifié Hexp?
  hobs<-read.table("Various_Statistics_postdisp.txt",sep="",skip=9,nrows=1)[,5]
  hexp<-read.table("Various_Statistics_postdisp.txt",sep="",skip=10,nrows=1)[,10]
  dat$f<-1-(hobs/hexp) # F  # différent de 0

  setwd("../")
  unlink(dir1, recursive =TRUE) 
  
  
  dat<-unlist(dat) # required to be read by Infusion
  return(dat) # returns the 6 summary statistics

}

