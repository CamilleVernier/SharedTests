### This script generates the six statistics used in this study from 
# a PSC model with a given set of parameter, by calling IBDsim from R.
# IBDsim v2.0 reference: Leblois et al. 2009 MolEcol Ressources

IBDSim_wrapper<-function(log10theta,
                 log10a,
                 log10tau,
                 nsim=1,
                 nloc=20, # number of loci
                 sampleSize=10, # number of individuals to simulate
                 mu=5e-2, # mutation rate
                 execName="./IBDSim"){ # executable name
  #conversion from the scaled parameters:
  a <- 10^log10a
  N0<-(10^log10theta)/mu
  t<-(10^log10tau)*N0
  
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
  system(execName,ignore.stdout=TRUE,ignore.stderr=TRUE)
  
  # We then load the simulated summary statistics and only keep the 6 used in this study:
  sumstats<-read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1)
  dat<-sumstats[,(dim(sumstats)[2]-3):dim(sumstats)[2]] # H, K, M and varK
  colnames(dat)<-c("H","VarK","K","M")
  dat$varH<-apply(sumstats[,1:20],1,var) #variance of H
  hobs<-read.table("Various_Statistics_postdisp.txt",sep="",skip=9,nrows=1)[,5]
  hexp<-read.table("Various_Statistics_postdisp.txt",sep="",skip=10,nrows=1)[,10] 
  dat$f<-1-(hobs/hexp) # F

  dat<-unlist(dat) # required to be read by Infusion
  return(dat) # returns the 6 summary statistics
}

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}