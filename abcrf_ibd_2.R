
rm(list = ls())
deb <- Sys.time()

###########################    IBDSim_wrapper_IBD   ###########################    

IBDSim_wrapper_IBD <- function(
  datafile=FALSE,
  lattice=c(70,70),
  samp=c(15,15),
  min_sample=c(15,15),
  D=1,
  nsim=1,
  nloc=20, # number of loci
  mu=5e-4, # mutation rate
  log10=FALSE,
  g=0.575, # geometric shape
  m=0.25, # total emigration rate
  habitatSize=NULL,
  dist_max=20,
  execName="../IBDSim"){ # executable name
  
  curDir<-getwd()
  dir1 <- tempfile(pattern = "sim", tmpdir =curDir )
  dir.create(dir1)
  setwd(dir = dir1)
  #message("Files in the currect directory: ", paste(list.files(path = ".", all.files = TRUE,recursive = TRUE, include.dirs = TRUE, no.. = TRUE)," "))
  
  
  Seed <- 1234567 + sample(1:10000, 1)
  data_file_name <- paste("Data_File",1234567 + sample(1:10000, 1), sep="")
  
  if ( ! is.null(habitatSize) ) {
    habitatSize <- floor(habitatSize);
    lattice <- c(habitatSize,habitatSize);
    min_sample=c( floor( floor(lattice[1]/2) - floor(samp[1]/2) ) ,floor( floor(lattice[2]/2) - floor(samp[2]/2) ) );
  }
  
  if( log10==TRUE )
  {
    g <- 10^g
    m <- 10^m 
  }
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
                    Geometric_Shape=",g,"
                    Total_Emigration_Rate=",m,"
                    MinDistReg=0.000001
                    Dist_max=",dist_max,"
                    
                    .txt_extension=true
                    Data_File_Name=",data_file_name,"
                    DiagnosticTables=Hexp,Fis,Iterative_Statistics,arRegression,erRegression,Iterative_Identity_Probability
                    noSS=T",sep=""),
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
                        Geometric_Shape=",g,"
                        Total_Emigration_Rate=",m,"
                        MinDistReg=0.000001
                        Dist_max=",dist_max,"
                        
                        Data_File_Name=",data_file_name,"
                        .txt_extension=true
                        DiagnosticTables=Hexp,Fis,Iterative_Statistics,arRegression,erRegression,Iterative_Identity_Probability
                        noSS=T",sep=""),
                  file="IbdSettings.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
      system(execName, ignore.stdout = TRUE)
    } else {
      dat<-sumstats[,(dim(sumstats)[2]-7-samp[1]):dim(sumstats)[2]] #on prend les stats qui nous int??ressent
      #colnames(dat) <- sumstats_name[,(dim(sumstats)[2]-8):dim(sumstats)[2]]
      qnames<-NULL
      for (i in 1:samp[1]) 
      {
        name <- paste("Qr", i-1, sep="")
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
      dat$Dsigma2 <- D*((m*(1+g))/(1-g)^2)
      data2 <- dat[,c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept",qnames)]  
      
      setwd("../")
      if(datafile=TRUE){file.copy(paste(data_file_name,".txt", sep=""), "./")}
      if(tries<3) unlink(dir1, recursive =TRUE) 
      
      if(dim(data2)[1]>1) {data2<-t(as.matrix(data2))} else data2<-unlist(data2) # required to be read by Infusion
      
      break
      
    }
  }
  return(data2)
}


#setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

###########################    DONNEES OBSERVEES   ###########################    

theta1_obs <- 0.575
theta2_obs <- 0.25
theta3_obs <- 70
s_obs <- IBDSim_wrapper_IBD(g=theta1_obs,m=theta2_obs,habitatSize=theta3_obs)
s_obs <- as.data.frame(t(s_obs))
# sobs.test <- matrix(s_obs,12,1)

###########################    CODE TABLE REFERENCE   ###########################    

# N <- 50 # taille table ref
# N <- 10000
# theta <- matrix(0,N,3) # matrice des parametres
# s <- matrix(0,N,12) # matrice des stats résumées
# 
# deb<- Sys.time()
# for (i in 1:N)
# {
#   theta[i,1] <- runif(1)
#   theta[i,2] <- runif(1)
#   theta[i,3] <- runif(1,16,100)
#   s[i,] <- IBDSim_wrapper_IBD(g=theta[i,1],m=theta[i,2],
#                              habitatSize=theta[i,3])
# print(i)  
# }
# fin <- Sys.time()
# duree <- fin - deb
# 
# training1 <- data.frame(theta1=theta[,1],s=s)
# colnames(training1) = c("g","Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept")
# training2 <- data.frame(theta2=theta[,2],s=s)
# training3 <- data.frame(theta3=theta[,3],s=s)
# 
# 
# ###########################    CODE TABLE TEST   ###########################    
# 
# # Ntest <- 10
# # Ntest <- 100
# 
# theta.test <- matrix(0,Ntest,3)
# s.test <- matrix(0,Ntest,12)
# for (i in 1:Ntest)
# {
#   theta.test[i,1] <- runif(1)
#   theta.test[i,2] <- runif(1)
#   theta.test[i,3] <- runif(1,16,100)
#   s.test[i,] <- IBDSim_wrapper_IBD(g=theta.test[i,1],m=theta.test[i,2],
#                                    habitatSize=theta.test[i,3])
#   print(i)  
# }
# 
# test <- data.frame(s=s.test)
# colnames(test) = c("Hobs_moy","varHobs","Hexp_moy","varHexp","fis_moy","fis","nb_allele_moyTotalSample","var_nballele","ar_slope","ar_intercept","er_slope","er_intercept")

###########################    TABLE REF PARALLELISATION   ###########################    

#N=100
 N <- 10000
nb_sum_stats <- length(s_obs)
theta <- matrix(0,N,3) # matrice des parametres
s <- matrix(0,N,nb_sum_stats) # matrice des stats résumées

#nbcores <-6
args = commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
print(paste(nbcores))


#setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")
for (i in 1:N)
{
  theta[i,1] <- runif(1)
  theta[i,2] <- runif(1)
  theta[i,3] <- runif(1,16,100)
}

library(doParallel)

cl <- makeCluster(nbcores) # Initialisationcluster
registerDoParallel(cl)   # Obligatoire avant la parallelisation

deb2 <- Sys.time()
train <- foreach(i=1:N, .combine="rbind") %dopar% 
{
  s[i,] <- IBDSim_wrapper_IBD(g=theta[i,1],m=theta[i,2],habitatSize=theta[i,3])
}
stopCluster(cl)    # Obligatoire après la parallelisation


fin2 <- Sys.time()
duree2 <- fin2 - deb2


train <- cbind(theta,train)
colnames(train)[1:3] <- c("g","m","habitatsize")
head(train) 
train <- as.data.frame(train)
train1 <- train[,-c(2,3)]
train2 <- train[,-c(1,3)]
train3 <- train[,-c(1,2)]

len <- length(train)

train_without_Q <- train[,-c(15:len)]
train1_without_Q <- train_without_Q[,-c(2,3)]
train2_without_Q <- train_without_Q[,-c(1,3)]
train3_without_Q <- train_without_Q[,-c(1,2)]

train_without_ar_er <- train[,-c(12:15)]
train1_without_ar_er <- train_without_ar_er[,-c(2,3)]
train2_without_ar_er <- train_without_ar_er[,-c(1,3)]
train3_without_ar_er <- train_without_ar_er[,-c(1,2)]

train_without_ar_er_Q <- train[,-c(12:len)]
train1_without_ar_er_Q <- train_without_ar_er_Q[,-c(2,3)]
train2_without_ar_er_Q <- train_without_ar_er_Q[,-c(1,3)]
train3_without_ar_er_Q <- train_without_ar_er_Q[,-c(1,2)]


###########################    TABLE TEST PARALLELISATION   ###########################    

# Ntest <- 10
Ntest <- 100

theta.test <- matrix(c(theta1_obs,theta2_obs, theta3_obs),Ntest,3, byrow=TRUE)
s.test <- matrix(0,Ntest,nb_sum_stats)

library(doParallel)

cl <- makeCluster(nbcores) # Initialisationcluster
registerDoParallel(cl)   # Obligatoire avant la parallelisation
test <- foreach(i=1:Ntest, .combine="rbind") %dopar% 
{
  s.test[i,] <- IBDSim_wrapper_IBD(datafile = TRUE, g=theta.test[i,1],m=theta.test[i,2],habitatSize=theta.test[i,3])
}
stopCluster(cl) 

test <- as.data.frame(test)

len_test <- length(test)
test_without_Q <- test[,-c(13:len_test)]
test_without_ar_er <- test[,-c(9:12)]
test_without_ar_er_Q <- test[,-c(9:len_test)]

###########################    RANDOM FOREST   ###########################    

library(abcrf)

#### Paramètre g ####
model1 <- regAbcrf(g~.,data=train1,ntree=500)
pred.theta1 <- predict(model1,test,train1)

# MSE1 <- mean(pred.theta1$expectation)- theta.test[,1]
# MSE1 

MSE1 <- mean((pred.theta1$expectation - theta.test[,1])^2)

plot(theta.test[,1],pred.theta1$med)
abline(a=0,b=1,col="red")

model1_without_Q <- regAbcrf(g~.,data=train1_without_Q,ntree=500)
pred.theta1_without_Q <- predict(model1_without_Q,test_without_Q,train1_without_Q)
MSE1_without_Q <- mean((theta.test[,1] - pred.theta1_without_Q$expectation)^2)
MSE1_without_Q

model1_without_ar_er <- regAbcrf(g~.,data=train1_without_ar_er,ntree=500)
pred.theta1_without_ar_er <- predict(model1_without_ar_er,test_without_ar_er,train1_without_ar_er)
#pred_theta1_without_ar_er_obs <- predict(model1_without_ar_er,s_obs,train1_without_ar_er)
MSE1_without_ar_er <- mean((theta.test[,1] - pred.theta1_without_ar_er$expectation)^2)
MSE1_without_ar_er

model1_without_ar_er_Q <- regAbcrf(g~.,data=train1_without_ar_er_Q,ntree=500)
pred.theta1_without_ar_er_Q <- predict(model1_without_ar_er_Q,test_without_ar_er_Q,train1_without_ar_er_Q)
#pred_theta1_without_ar_er_obs <- predict(model1_without_ar_er,s_obs,train1_without_ar_er)
MSE1_without_ar_er_Q <- mean((theta.test[,1] - pred.theta1_without_ar_er_Q$expectation)^2)
MSE1_without_ar_er_Q




#### Paramètre m ####
model2 <- regAbcrf(m~.,data=train2,ntree=500)
pred.theta2 <- predict(model2,test,train2)

plot(theta.test[,2],pred.theta2$med)
abline(a=0,b=1,col="red")

MSE2 <- mean((theta.test[,2] - pred.theta2$expectation)^2)
MSE2 

model2_without_Q <- regAbcrf(m~.,data=train2_without_Q,ntree=500)
pred.theta2_without_Q <- predict(model2_without_Q,test_without_Q,train2_without_Q)
pred_theta2_without_Q_obs <- predict(model2_without_Q,s_obs,train2_without_Q)
MSE2_without_Q <- mean((theta.test[,2] - pred.theta2_without_Q$expectation)^2)
MSE2_without_Q

model2_without_ar_er <- regAbcrf(m~.,data=train2_without_ar_er,ntree=500)
pred.theta2_without_ar_er <- predict(model2_without_ar_er,test_without_ar_er,train2_without_ar_er)
#pred_theta2_without_ar_er_obs <- predict(model2_without_ar_er,s_obs,train2_without_ar_er)
MSE2_without_ar_er <- mean((theta.test[,2] - pred.theta2_without_ar_er$expectation)^2)
MSE2_without_ar_er

model2_without_ar_er_Q <- regAbcrf(m~.,data=train2_without_ar_er_Q,ntree=500)
pred.theta2_without_ar_er_Q <- predict(model2_without_ar_er_Q,test_without_ar_er_Q,train2_without_ar_er_Q)
#pred_theta2_without_ar_er_obs <- predict(model2_without_ar_er,s_obs,train2_without_ar_er)
MSE2_without_ar_er_Q <- mean((theta.test[,2] - pred.theta2_without_ar_er_Q$expectation)^2)
MSE2_without_ar_er_Q


#### Paramètre habitatsize ####
model3 <- regAbcrf(habitatsize~.,data=train3,ntree=500)
pred.theta3 <- predict(model3,test,train3)

plot(theta.test[,3],pred.theta3$med)
abline(a=0,b=1,col="red")

MSE3 <- mean((theta.test[,3] - pred.theta3$expectation)^2)
MSE3 

model3_without_Q <- regAbcrf(habitatsize~.,data=train3_without_Q,ntree=500)
pred.theta3_without_Q <- predict(model3_without_Q,test_without_Q,train3_without_Q)
pred_theta3_without_Q_obs <- predict(model3_without_Q,s_obs,train3_without_Q)
MSE3_without_Q <- mean((theta.test[,3] - pred.theta3_without_Q$expectation)^2)
MSE3_without_Q

model3_without_ar_er <- regAbcrf(habitatsize~.,data=train3_without_ar_er,ntree=500)
pred.theta3_without_ar_er <- predict(model3_without_ar_er,test_without_ar_er,train3_without_ar_er)
#pred_theta3_without_ar_er_obs <- predict(model3_without_ar_er,s_obs,train3_without_ar_er)
MSE3_without_ar_er <- mean((theta.test[,3] - pred.theta3_without_ar_er$expectation)^2)
MSE3_without_ar_er

model3_without_ar_er_Q <- regAbcrf(habitatsize~.,data=train3_without_ar_er_Q,ntree=500)
pred.theta3_without_ar_er_Q <- predict(model3_without_ar_er_Q,test_without_ar_er_Q,train3_without_ar_er_Q)
#pred_theta3_without_ar_er_obs <- predict(model3_without_ar_er,s_obs,train3_without_ar_er)
MSE3_without_ar_er_Q <- mean((theta.test[,3] - pred.theta3_without_ar_er_Q$expectation)^2)
MSE3_without_ar_er_Q

########################    Biais    ########################


biais1 <- mean(pred.theta1$expectation) - theta1_obs
biais1_without_Q <- mean(pred.theta1_without_Q$expectation) - theta1_obs
biais1_without_ar_er <- mean(pred.theta1_without_ar_er$expectation) - theta1_obs
biais1_without_ar_er_Q <- mean(pred.theta1_without_ar_er_Q$expectation) - theta1_obs

biais2 <- mean(pred.theta2$expectation) - theta2_obs
biais2_without_Q <- mean(pred.theta2_without_Q$expectation) - theta2_obs
biais2_without_ar_er <- mean(pred.theta2_without_ar_er$expectation) - theta2_obs
biais2_without_ar_er_Q <- mean(pred.theta2_without_ar_er_Q$expectation) - theta2_obs


biais3 <- mean(pred.theta3$expectation) - theta3_obs
biais3_without_Q <- mean(pred.theta3_without_Q$expectation) - theta3_obs
biais3_without_ar_er <- mean(pred.theta3_without_ar_er$expectation) - theta3_obs
biais3_without_ar_er_Q <- mean(pred.theta3_without_ar_er_Q$expectation) - theta3_obs


# biais1 <- mean(pred.theta1$expectation - theta[,1])
# biais1_without_Q <- mean(pred.theta1_without_Q$expectation - theta[,1])
# biais1_without_ar_er <- mean(pred.theta1_without_ar_er$expectation - theta[,1])
# biais1_without_ar_er_Q <- mean(pred.theta1_without_ar_er_Q$expectation - theta[,1])
# 
# biais2 <- mean(pred.theta2$expectation - theta[,2])
# biais2_without_Q <- mean(pred.theta2_without_Q$expectation - theta[,2])
# biais2_without_ar_er <- mean(pred.theta2_without_ar_er$expectation - theta[,2])
# biais2_without_ar_er_Q <- mean(pred.theta2_without_ar_er_Q$expectation - theta[,2])
# 
# 
# biais3 <- mean(pred.theta3$expectation - theta[,3])
# biais3_without_Q <- mean(pred.theta3_without_Q$expectation - theta[,3])
# biais3_without_ar_er <- mean(pred.theta3_without_ar_er$expectation - theta[,3])
# biais3_without_ar_er_Q <- mean(pred.theta3_without_ar_er_Q$expectation - theta[,3])

########################    NMAE    ########################
for (i in 1:3) {
  name <- paste("NMAE",i,sep="")
  tr_n <- paste("train",i, sep="")
  tr <- eval(parse(text=tr_n))
  pr_n <- paste("pred.theta",i,"$expectation", sep="")
  pr <- eval(parse(text=pr_n))
  nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
  assign(name, nmae)
  
  name <- paste("NMAE",i,"_without_Q",sep="")
  tr_n <- paste("train",i,"_without_Q", sep="")
  tr <- eval(parse(text=tr_n))
  pr_n <- paste("pred.theta",i,"_without_Q$expectation", sep="")
  pr <- eval(parse(text=pr_n))
  nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
  assign(name, nmae)
  
  name <- paste("NMAE",i,"_without_ar_er",sep="")
  tr_n <- paste("train",i,"_without_ar_er", sep="")
  tr <- eval(parse(text=tr_n))
  pr_n <- paste("pred.theta",i,"_without_ar_er$expectation", sep="")
  pr <- eval(parse(text=pr_n))
  nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
  assign(name, nmae)
  
  name <- paste("NMAE",i,"_without_ar_er_Q",sep="")
  tr_n <- paste("train",i,"_without_ar_er_Q", sep="")
  tr <- eval(parse(text=tr_n))
  pr_n <- paste("pred.theta",i,"_without_ar_er_Q$expectation", sep="")
  pr <- eval(parse(text=pr_n))
  nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
  assign(name, nmae)
}
# NMAE <- mean(abs((train3_without_Q[,1]-pred.theta3_without_Q$expectation)/train3_without_Q[,1]))

########################    PREDICT OOB    ########################

oob1 <- predictOOB(model1, train1)
oob1_without_Q <- predictOOB(model1_without_Q, train1_without_Q, paral=FALSE)
oob1_without_ar_er <- predictOOB(model1_without_ar_er, train1_without_ar_er, paral = FALSE)
oob1_without_ar_er_Q <- predictOOB(model1_without_ar_er_Q, train1_without_ar_er_Q, paral=FALSE)

OOB_results1 <- matrix(c(oob1$MSE,oob1_without_Q$MSE, oob1_without_ar_er$MSE, oob1_without_ar_er_Q$MSE,
                    oob1$NMAE,oob1_without_Q$NMAE, oob1_without_ar_er$NMAE, oob1_without_ar_er_Q$NMAE,
                    oob1$med_MSE,oob1_without_Q$med_MSE, oob1_without_ar_er$med_MSE, oob1_without_ar_er_Q$med_MSE,
                    oob1$med_NMAE,oob1_without_Q$med_NMAE, oob1_without_ar_er$med_NMAE, oob1_without_ar_er_Q$med_NMAE),4,4, byrow=FALSE)
colnames(OOB_results1) <- c("All","Without Qij", "Without ar,er", "Without ar,er,Qij")
rownames(OOB_results1) <- c("MSE", "NMAE", "med MSE", "med NMAE")


oob2 <- predictOOB(model2, train2, paral = FALSE)
oob2_without_Q <- predictOOB(model2_without_Q, train2_without_Q, paral=FALSE)
oob2_without_ar_er <- predictOOB(model2_without_ar_er, train2_without_ar_er, paral = FALSE)
oob2_without_ar_er_Q <- predictOOB(model2_without_ar_er_Q, train2_without_ar_er_Q, paral=FALSE)

OOB_results2 <- matrix(c(oob2$MSE,oob2_without_Q$MSE, oob2_without_ar_er$MSE, oob2_without_ar_er_Q$MSE,
                         oob2$NMAE,oob2_without_Q$NMAE, oob2_without_ar_er$NMAE, oob2_without_ar_er_Q$NMAE,
                         oob2$med_MSE,oob2_without_Q$med_MSE, oob2_without_ar_er$med_MSE, oob2_without_ar_er_Q$med_MSE,
                         oob2$med_NMAE,oob2_without_Q$med_NMAE, oob2_without_ar_er$med_NMAE, oob2_without_ar_er_Q$med_NMAE),4,4, byrow=FALSE)
colnames(OOB_results2) <- c("All","Without Qij", "Without ar,er", "Without ar,er,Qij")
rownames(OOB_results2) <- c("MSE", "NMAE", "med MSE", "med NMAE")



oob3 <- predictOOB(model3, train3, paral = FALSE)
oob3_without_Q <- predictOOB(model3_without_Q, train3_without_Q, paral=FALSE)
oob3_without_ar_er <- predictOOB(model3_without_ar_er, train3_without_ar_er, paral = FALSE)
oob3_without_ar_er_Q <- predictOOB(model3_without_ar_er_Q, train3_without_ar_er_Q, paral=FALSE)

OOB_results3 <- matrix(c(oob3$MSE,oob3_without_Q$MSE, oob3_without_ar_er$MSE, oob3_without_ar_er_Q$MSE,
                         oob3$NMAE,oob3_without_Q$NMAE, oob3_without_ar_er$NMAE, oob3_without_ar_er_Q$NMAE,
                         oob3$med_MSE,oob3_without_Q$med_MSE, oob3_without_ar_er$med_MSE, oob3_without_ar_er_Q$med_MSE,
                         oob3$med_NMAE,oob3_without_Q$med_NMAE, oob3_without_ar_er$med_NMAE, oob3_without_ar_er_Q$med_NMAE),4,4, byrow=FALSE)
colnames(OOB_results3) <- c("All","Without Qij", "Without ar,er", "Without ar,er,Qij")
rownames(OOB_results3) <- c("OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")



# densityPlot(model1,test[1:100,],train1)
densityPlot(model1,test[10,],train1)
plot(model1)
plot(model2)
plot(model3)
err.regAbcrf(model1, train1, paral=F)

op <- par(no.readonly = TRUE) 
par(mfrow=c(2, 2)) 
plot(model1, main="Estimation du paramètre g")
plot(model2, main="Estimation du paramètre m")
plot(model3, main="Estimation du paramètre habitatsize")
par(mfrow=c(1, 1))
par(op)
plot(model1_without_Q, main="Estimation du paramètre g")


#########################   RESULTS   ######################### 

results1 <- rbind(c(mean(pred.theta1$expectation), mean(pred.theta1_without_Q$expectation), mean(pred.theta1_without_ar_er$expectation), mean(pred.theta1_without_ar_er_Q$expectation)),
                  c(biais1, biais1_without_Q, biais1_without_ar_er, biais1_without_ar_er_Q),
                  c(var(pred.theta1$expectation), var(pred.theta1_without_Q$expectation), var(pred.theta1_without_ar_er$expectation), var(pred.theta1_without_ar_er_Q$expectation)),
                  c(MSE1, MSE1_without_Q, MSE1_without_ar_er,MSE1_without_ar_er_Q),                  
                  c(NMAE1, NMAE1_without_Q, NMAE1_without_ar_er,NMAE1_without_ar_er_Q),
                  OOB_results1)
rownames(results1) <- c("Valeur moyenne","Biais","Variance", "MSE","NMAE", "OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")

results2 <- rbind(c(mean(pred.theta2$expectation), mean(pred.theta2_without_Q$expectation), mean(pred.theta2_without_ar_er$expectation), mean(pred.theta2_without_ar_er_Q$expectation)),
                  c(biais2, biais2_without_Q, biais2_without_ar_er, biais2_without_ar_er_Q),
                  c(var(pred.theta2$expectation), var(pred.theta2_without_Q$expectation), var(pred.theta2_without_ar_er$expectation), var(pred.theta2_without_ar_er_Q$expectation)),
                  c(MSE2, MSE2_without_Q, MSE2_without_ar_er,MSE2_without_ar_er_Q),
                  c(NMAE2, NMAE2_without_Q, NMAE2_without_ar_er,NMAE2_without_ar_er_Q),
                  OOB_results2)
rownames(results2) <- c("Valeur moyenne","biais", "Variance","MSE", "NMAE","OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")


results3 <- rbind(c(mean(pred.theta3$expectation), mean(pred.theta3_without_Q$expectation), mean(pred.theta3_without_ar_er$expectation), mean(pred.theta3_without_ar_er_Q$expectation)),
                  c(biais3, biais3_without_Q, biais3_without_ar_er, biais3_without_ar_er_Q),
                  c(var(pred.theta3$expectation), var(pred.theta3_without_Q$expectation), var(pred.theta3_without_ar_er$expectation), var(pred.theta3_without_ar_er_Q$expectation)),
                  c(MSE3, MSE3_without_Q, MSE3_without_ar_er,MSE3_without_ar_er_Q),
                  c(NMAE3, NMAE3_without_Q, NMAE3_without_ar_er,NMAE3_without_ar_er_Q),
                  OOB_results3)
rownames(results3) <- c("Valeur moyenne","biais", "Variance","MSE", "NMAE","OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")

#########################   PREDICT OBS   ######################### 

pred_theta1_obs <- predict(model1,s_obs,train1)
pred_theta2_obs <- predict(model2,s_obs,train2)
pred_theta3_obs <- predict(model3,s_obs,train3)

######################  Sauvegardes Rdata  ###################### 

write.csv(results1, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results1_100_meme_para_gmh.csv", sep=""))
write.csv(results2, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results2_100_meme_para_gmh.csv", sep=""))
write.csv(results3, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results3_100_meme_para_gmh.csv", sep=""))

save.image(file=paste("/home/vernierc/Documents/Logiciels/ABCRF/Rdata/",gsub(" ", "_", deb),".Rdata", sep=""))
print(deb)
fin <- Sys.time()
duree <- fin - deb
#calcul de l'OOB error








# 
# 
# 
# 
# 
# ######################  Juste g et m  ###################### 
# 
# theta1_obs <- 0.575
# theta2_obs <- 0.25
# 
# s_obs <- IBDSim_wrapper_IBD(g=theta1_obs,m=theta2_obs)
# s_obs <- as.data.frame(t(s_obs))
# 
# #setwd(dir="/home/vernierc/Documents/Logiciels/ABCRF/")
# 
# N <- 10000
# nb_sum_stats <- length(s_obs)
# theta <- matrix(0,N,2) # matrice des parametres
# s <- matrix(0,N,nb_sum_stats) # matrice des stats résumées
# 
# #nbcores <-6
# 
# #setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")
# for (i in 1:N)
# {
#   theta[i,1] <- runif(1)
#   theta[i,2] <- runif(1)
# }
# 
# library(doParallel)
# 
# cl <- makeCluster(nbcores) 
# registerDoParallel(cl)   
# 
# train <- foreach(i=1:N, .combine="rbind") %dopar% 
# {
#   s[i,] <- IBDSim_wrapper_IBD(g=theta[i,1],m=theta[i,2])
# }
# stopCluster(cl)    
# 
# 
# train <- cbind(theta,train)
# colnames(train)[1:2] <- c("g","m")
# head(train) 
# train <- as.data.frame(train)
# train1 <- train[,-2]
# train2 <- train[,-1]
# 
# len <- length(train)
# 
# train_without_Q <- train[,-c(14:len)]
# train1_without_Q <- train_without_Q[,-2]
# train2_without_Q <- train_without_Q[,-1]
# 
# train_without_ar_er <- train[,-c(11:14)]
# train1_without_ar_er <- train_without_ar_er[,-2]
# train2_without_ar_er <- train_without_ar_er[,-1]
# 
# train_without_ar_er_Q <- train[,-c(11:len)]
# train1_without_ar_er_Q <- train_without_ar_er_Q[,-2]
# train2_without_ar_er_Q <- train_without_ar_er_Q[,-1]
# 
# 
# Ntest <- 100
# 
# theta.test <- matrix(c(theta1_obs,theta2_obs),Ntest,2, byrow=TRUE)
# s.test <- matrix(0,Ntest,nb_sum_stats)
# 
# library(doParallel)
# 
# cl <- makeCluster(nbcores) 
# registerDoParallel(cl)   
# test <- foreach(i=1:Ntest, .combine="rbind") %dopar% 
# {
#   s.test[i,] <- IBDSim_wrapper_IBD(g=theta.test[i,1],m=theta.test[i,2])
# }
# stopCluster(cl) 
# 
# test <- as.data.frame(test)
# 
# len_test <- length(test)
# test_without_Q <- test[,-c(13:len_test)]
# test_without_ar_er <- test[,-c(9:12)]
# test_without_ar_er_Q <- test[,-c(9:len_test)]
#  ###########
# 
# library(abcrf)
# 
# model1 <- regAbcrf(g~.,data=train1,ntree=500)
# pred.theta1 <- predict(model1,test,train1)
# 
# MSE1 <- mean((pred.theta1$expectation - theta.test[,1])^2)
# 
# plot(theta.test[,1],pred.theta1$med)
# abline(a=0,b=1,col="red")
# 
# model1_without_Q <- regAbcrf(g~.,data=train1_without_Q,ntree=500)
# pred.theta1_without_Q <- predict(model1_without_Q,test_without_Q,train1_without_Q)
# MSE1_without_Q <- mean((theta.test[,1] - pred.theta1_without_Q$expectation)^2)
# MSE1_without_Q
# 
# 
# model1_without_ar_er <- regAbcrf(g~.,data=dr,ntree=500)
# pred.theta1_without_ar_er <- predict(model1_without_ar_er,test_without_ar_er,train1_without_ar_er)
# #pred_theta1_without_ar_er_obs <- predict(model1_without_ar_er,s_obs,train1_without_ar_er)
# MSE1_without_ar_er <- mean((theta.test[,1] - pred.theta1_without_ar_er$expectation)^2)
# MSE1_without_ar_er
# 
# model1_without_ar_er_Q <- regAbcrf(g~.,data=train1_without_ar_er_Q,ntree=500)
# pred.theta1_without_ar_er_Q <- predict(model1_without_ar_er_Q,test_without_ar_er_Q,train1_without_ar_er_Q)
# #pred_theta1_without_ar_er_obs <- predict(model1_without_ar_er,s_obs,train1_without_ar_er)
# MSE1_without_ar_er_Q <- mean((theta.test[,1] - pred.theta1_without_ar_er_Q$expectation)^2)
# MSE1_without_ar_er_Q
# 
# 
# model2 <- regAbcrf(m~.,data=train2,ntree=500)
# pred.theta2 <- predict(model2,test,train2)
# 
# plot(theta.test[,2],pred.theta2$med)
# abline(a=0,b=1,col="red")
# 
# MSE2 <- mean((theta.test[,2] - pred.theta2$expectation)^2)
# MSE2 
# 
# model2_without_Q <- regAbcrf(m~.,data=train2_without_Q,ntree=500)
# pred.theta2_without_Q <- predict(model2_without_Q,test_without_Q,train2_without_QE)
# #pred_theta2_without_Q_obs <- predict(model2_without_Q,s_obs,train2_without_Q)
# MSE2_without_Q <- mean((theta.test[,2] - pred.theta2_without_Q$expectation)^2)
# MSE2_without_Q
# 
# model2_without_ar_er <- regAbcrf(m~.,data=train2_without_ar_er,ntree=500)
# pred.theta2_without_ar_er <- predict(model2_without_ar_er,test_without_ar_er,train2_without_ar_er)
# #pred_theta2_without_ar_er_obs <- predict(model2_without_ar_er,s_obs,train2_without_ar_er)
# MSE2_without_ar_er <- mean((theta.test[,2] - pred.theta2_without_ar_er$expectation)^2)
# MSE2_without_ar_er
# 
# model2_without_ar_er_Q <- regAbcrf(m~.,data=train2_without_ar_er_Q,ntree=500)
# pred.theta2_without_ar_er_Q <- predict(model2_without_ar_er_Q,test_without_ar_er_Q,train2_without_ar_er_Q)
# #pred_theta2_without_ar_er_obs <- predict(model2_without_ar_er,s_obs,train2_without_ar_er)
# MSE2_without_ar_er_Q <- mean((theta.test[,2] - pred.theta2_without_ar_er_Q$expectation)^2)
# MSE2_without_ar_er_Q
# 
# 
# for (i in 1:2) {
#   name <- paste("NMAE",i,sep="")
#   tr_n <- paste("train",i, sep="")
#   tr <- eval(parse(text=tr_n))
#   pr_n <- paste("pred.theta",i,"$expectation", sep="")
#   pr <- eval(parse(text=pr_n))
#   nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
#   assign(name, nmae)
#   
#   name <- paste("NMAE",i,"_without_Q",sep="")
#   tr_n <- paste("train",i,"_without_Q", sep="")
#   tr <- eval(parse(text=tr_n))
#   pr_n <- paste("pred.theta",i,"_without_Q$expectation", sep="")
#   pr <- eval(parse(text=pr_n))
#   nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
#   assign(name, nmae)
#   
#   name <- paste("NMAE",i,"_without_ar_er",sep="")
#   tr_n <- paste("train",i,"_without_ar_er", sep="")
#   tr <- eval(parse(text=tr_n))
#   pr_n <- paste("pred.theta",i,"_without_ar_er$expectation", sep="")
#   pr <- eval(parse(text=pr_n))
#   nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
#   assign(name, nmae)
#   
#   name <- paste("NMAE",i,"_without_ar_er_Q",sep="")
#   tr_n <- paste("train",i,"_without_ar_er_Q", sep="")
#   tr <- eval(parse(text=tr_n))
#   pr_n <- paste("pred.theta",i,"_without_ar_er_Q$expectation", sep="")
#   pr <- eval(parse(text=pr_n))
#   nmae <- mean(abs((tr[,1]-pr)/tr[,1]))
#   assign(name, nmae)
# }
# # NMAE <- mean(abs((train3_without_Q[,1]-pred.theta3_without_Q$expectation)/train3_without_Q[,1]))
# 
# oob1 <- predictOOB(model1, train1)
# oob1_without_Q <- predictOOB(model1_without_Q, train1_without_Q)
# oob1_without_ar_er <- predictOOB(model1_without_ar_er, train1_without_ar_er)
# oob1_without_ar_er_Q <- predictOOB(model1_without_ar_er_Q, train1_without_ar_er_Q)
# 
# OOB_results1 <- matrix(c(oob1$MSE,oob1_without_Q$MSE, oob1_without_ar_er$MSE, oob1_without_ar_er_Q$MSE,
#                          oob1$NMAE,oob1_without_Q$NMAE, oob1_without_ar_er$NMAE, oob1_without_ar_er_Q$NMAE,
#                          oob1$med_MSE,oob1_without_Q$med_MSE, oob1_without_ar_er$med_MSE, oob1_without_ar_er_Q$med_MSE,
#                          oob1$med_NMAE,oob1_without_Q$med_NMAE, oob1_without_ar_er$med_NMAE, oob1_without_ar_er_Q$med_NMAE),4,4, byrow=TRUE)
# colnames(OOB_results1) <- c("All","Without Qij", "Without ar,er", "Without ar,er,Qij")
# rownames(OOB_results1) <- c("MSE", "NMAE", "med MSE", "med NMAE")
# 
# 
# oob2 <- predictOOB(model2, train2)
# oob2_without_Q <- predictOOB(model2_without_Q, train2_without_Q)
# oob2_without_ar_er <- predictOOB(model2_without_ar_er, train2_without_ar_er)
# oob2_without_ar_er_Q <- predictOOB(model2_without_ar_er_Q, train2_without_ar_er_Q)
# 
# OOB_results2 <- matrix(c(oob2$MSE,oob2_without_Q$MSE, oob2_without_ar_er$MSE, oob2_without_ar_er_Q$MSE,
#                          oob2$NMAE,oob2_without_Q$NMAE, oob2_without_ar_er$NMAE, oob2_without_ar_er_Q$NMAE,
#                          oob2$med_MSE,oob2_without_Q$med_MSE, oob2_without_ar_er$med_MSE, oob2_without_ar_er_Q$med_MSE,
#                          oob2$med_NMAE,oob2_without_Q$med_NMAE, oob2_without_ar_er$med_NMAE, oob2_without_ar_er_Q$med_NMAE),4,4, byrow=TRUE)
# colnames(OOB_results2) <- c("All","Without Qij", "Without ar,er", "Without ar,er,Qij")
# rownames(OOB_results2) <- c("MSE", "NMAE", "med MSE", "med NMAE")
# 
# 
# # densityPlot(model1,test[1:100,],train1)
# densityPlot(model1,test[10,],train1)
# plot(model1)
# plot(model2)
# err.regAbcrf(model1, train1, paral=F)
# 
# op <- par(no.readonly = TRUE) 
# par(mfrow=c(2, 1)) 
# plot(model1, main="Estimation du paramètre g")
# plot(model2, main="Estimation du paramètre m")
# par(mfrow=c(1, 1))
# par(op)
# plot(model1_without_Q, main="Estimation du paramètre g")
# 
# biais1 <- mean(pred.theta1$expectation) - theta1_obs
# biais1_without_Q <- mean(pred.theta1_without_Q$expectation) - theta1_obs
# biais1_without_ar_er <- mean(pred.theta1_without_ar_er$expectation) - theta1_obs
# biais1_without_ar_er_Q <- mean(pred.theta1_without_ar_er_Q$expectation) - theta1_obs
# 
# biais2 <- mean(pred.theta2$expectation) - theta2_obs
# biais2_without_Q <- mean(pred.theta2_without_Q$expectation) - theta2_obs
# biais2_without_ar_er <- mean(pred.theta2_without_ar_er$expectation) - theta2_obs
# biais2_without_ar_er_Q <- mean(pred.theta2_without_ar_er_Q$expectation) - theta2_obs
# 
# 
# results1 <- rbind(c(mean(pred.theta1$expectation), mean(pred.theta1_without_Q$expectation), mean(pred.theta1_without_ar_er$expectation), mean(pred.theta1_without_ar_er_Q$expectation)),
#                   c(biais1, biais1_without_Q, biais1_without_ar_er, biais1_without_ar_er_Q),
#                   c(var(pred.theta1$expectation), var(pred.theta1_without_Q$expectation), var(pred.theta1_without_ar_er$expectation), var(pred.theta1_without_ar_er_Q$expectation)),
#                   c(MSE1, MSE1_without_Q, MSE1_without_ar_er,MSE1_without_ar_er_Q),                  
#                   c(NMAE1, NMAE1_without_Q, NMAE1_without_ar_er,NMAE1_without_ar_er_Q),
#                   OOB_results1)
# rownames(results1) <- c("Valeur moyenne","Biais","Variance", "MSE","NMAE", "OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")
# 
# results2 <- rbind(c(mean(pred.theta2$expectation), mean(pred.theta2_without_Q$expectation), mean(pred.theta2_without_ar_er$expectation), mean(pred.theta2_without_ar_er_Q$expectation)),
#                   c(biais2, biais2_without_Q, biais2_without_ar_er, biais2_without_ar_er_Q),
#                   c(var(pred.theta2$expectation), var(pred.theta2_without_Q$expectation), var(pred.theta2_without_ar_er$expectation), var(pred.theta2_without_ar_er_Q$expectation)),
#                   c(MSE2, MSE2_without_Q, MSE2_without_ar_er,MSE2_without_ar_er_Q),
#                   c(NMAE2, NMAE2_without_Q, NMAE2_without_ar_er,NMAE2_without_ar_er_Q),
#                   OOB_results2)
# rownames(results2) <- c("Valeur moyenne","biais", "Variance","MSE", "NMAE","OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")
# 
# 
# write.csv(results1, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results_gm_1.csv", sep=""))
# write.csv(results2, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results_gm_2.csv", sep=""))
# 
# save.image(file=paste("/home/vernierc/Documents/Logiciels/ABCRF/Rdata/",gsub(" ", "_", deb),"gm.Rdata", sep=""))
# print(deb)
# fin <- Sys.time()
# duree <- fin - deb