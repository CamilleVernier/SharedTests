### This script generates the six statistics used in this study from 
# a PSC model with a given set of parameter, by calling IBDsim from R.
# IBDsim v2.0 reference: Leblois et al. 2009 MolEcol Ressources

IBDSim_wrapper_IBD <-function(
                         lattice=c(20,20),
                         samp=c(10,10),
                         min_sample=c(5,5),
                         D=1,
                         nsim=1,
                         nloc=40, # number of loci
                         mu=5e-3, # mutation rate
                         log10gshape=log10(0.25),# moche
                         log10m=log10(0.45), # très moche
                         g_shape=0.25, # geometric shape
                         m=0.45, # total emigration rate
                         dist_max=10,
                         execName="../IBDSim"){ # executable name

  curDir<-getwd()
  dir1 <- tempfile(pattern = "sim", tmpdir =curDir )
  dir.create(dir1)
  setwd(dir = dir1)
  #message("Files in the currect directory: ", paste(list.files(path = ".", all.files = TRUE,recursive = TRUE, include.dirs = TRUE, no.. = TRUE)," "))
  
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
                    Min_Allele_Number=2

                    %%%%%%%% DEMOGRAPHIC OPTIONS %%%%%%%%%%%%%
                    %% LATTICE
                    Lattice_SizeX=",lattice[1],"
                    Lattice_SizeY=",lattice[2],"
                    Ind_Per_Pop=",D,"

                    %% SAMPLE
                    Sample_SizeX=",samp[1],"
                    Sample_SizeY=",samp[2],"
                    Min_Sample_CoordinateX=",min_sample[1],"
                    Min_Sample_CoordinateY=",min_sample[2],"
                    Ind_Per_Pop_Sampled=1

                    %% DISPERSAL
                    Dispersal_Distribution=g
                    Geometric_Shape=",g_shape,"
                    Total_Emigration_Rate=",m,"
                    MinDistReg=0.000001
                    Dist_max=",dist_max,"
                    
                    DiagnosticTables=Hexp,Fis,Iterative_Statistics,arRegression,erRegression,Iterative_Identity_Probability",sep=""),
              file="IbdSettings.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  system(execName, ignore.stdout = TRUE)
  #sumstats_name <- as.matrix(read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="", nrows=1))
  #sumstats_name2 <- as.vector(sumstats_name[1,])
  
  tries <- 1
  while (tries<4) {
    sumstats<-try(read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1))
    if(inherits(sumstats,"try-error")) {
      message(paste("Reading simulated summary statistics failed after ", tries, " tries."))
      tries <- tries + 1;
      Seed <- ( 1234567 + sample(1:10000, 1) );
      
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
                        Min_Allele_Number=2
                        
                        %%%%%%%% DEMOGRAPHIC OPTIONS %%%%%%%%%%%%%
                        %% LATTICE
                        Lattice_SizeX=",lattice[1],"
                        Lattice_SizeY=",lattice[2],"
                        Ind_Per_Pop=",D,"
                        
                        %% SAMPLE
                        Sample_SizeX=",samp[1],"
                        Sample_SizeY=",samp[2],"
                        Min_Sample_CoordinateX=",min_sample[1],"
                        Min_Sample_CoordinateY=",min_sample[2],"
                        Ind_Per_Pop_Sampled=1
                        
                        %% DISPERSAL
                        Dispersal_Distribution=g
                        Geometric_Shape=",g_shape,"
                        Total_Emigration_Rate=",m,"
                        MinDistReg=0.000001
                        Dist_max=",dist_max,"
                        
                        DiagnosticTables=Hexp,Fis,Iterative_Statistics,arRegression,erRegression,Iterative_Identity_Probability",sep=""),
                  file="IbdSettings.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
      system(execName, ignore.stdout = TRUE)
    } else {
      dat<-sumstats[,(dim(sumstats)[2]-7-samp[1]):dim(sumstats)[2]] #on prend les stats qui nous intéressent
      #colnames(dat) <- sumstats_name[,(dim(sumstats)[2]-8):dim(sumstats)[2]]
      qnames<-NULL
      for (i in 1:samp[1]) 
      {
        name <- paste("Qr(", i-1, ")", sep="")
        qnames <- c(qnames,name)
      }
      
      colnames(dat)<-c("Hobs_moy","Hexp_moy","fis_moy","nb_allele_moyTotalSample",qnames
                       , "ar_slope","ar_intercept","er_slope","er_intercept")
      
      hobs<-read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1)[,1:nloc]
      hexp<-read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1)[,(nloc+1):(2*nloc)]  
      dat$varHobs<-apply(hobs,1,var)
      dat$varHexp<-apply(hexp,1,var)
      hobs_moy <- as.numeric(dat[,1])
      hexp_moy <- as.numeric(dat[,2])
      dat$fis <- 1-(hobs_moy/hexp_moy)
      nballele <- read.table("Iterative_Statistics_postdisp_PerLocus.txt",sep="",skip=1)[,(3*nloc+1):(4*nloc)]  
      dat$var_nballele <- apply(nballele,1,var)
      #dat$Dsigma2 <- D*((m*(1+g_shape))/(1-g_shape)^2)
      data2 <- dat[,c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept")]  
      
      setwd("../")
      if(tries<3) unlink(dir1, recursive =TRUE) 
      
      if(dim(data2)[1]>1) {data2<-t(as.matrix(data2))} else data2<-unlist(data2) # required to be read by Infusion
      
      break
      
    }
  }
  
  return(data2)
}
