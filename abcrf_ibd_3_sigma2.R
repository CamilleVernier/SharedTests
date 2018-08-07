
rm(list = ls())
deb <- Sys.time()
deb <- gsub(" ", "_", deb)
name_dir_sys <- paste("./",deb,sep="")
dir.create(name_dir_sys)
para = TRUE
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
  data_file_name <- paste("Data_File",Seed, sep="")
  
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
                    
                    
                    Data_File_Name=",data_file_name,"
                    .txt_extension=true
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
            
      
      if(datafile==TRUE){file.copy(from=paste("./",data_file_name,".txt", sep=""), to=paste(".",name_dir_sys,sep=""))}
      setwd("../")
      if(tries<3) unlink(dir1, recursive =TRUE) 
      
      if(dim(data2)[1]>1) {data2<-t(as.matrix(data2))} else data2<-unlist(data2) # required to be read by Infusion
      
      break
      
    }
  }
  return(data2)
}


#setwd(dir="/home/vernierc/Documents/GitCamille/SharedTests/")

###########################    DONNEES OBSERVEES   ###########################    
#theta1_obs <- 0.575
theta1_obs <- 0.558615
theta2_obs <- 0.25
theta3_obs <- 70
theta4_obs <- (theta2_obs*(1+theta1_obs))/(1-theta1_obs)^2
theta4_obs
s_obs <- IBDSim_wrapper_IBD(g=theta1_obs,m=theta2_obs,habitatSize=theta3_obs)
s_obs <- as.data.frame(t(s_obs))
# sobs.test <- matrix(s_obs,12,1)


###########################    TABLE REF PARALLELISATION   ###########################    

N <- 2000
nb_sum_stats <- length(s_obs)
theta <- matrix(0,N,4) # matrice des parametres
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
  theta[i,4] <- (theta[i,2]*(1+theta[i,1]))/(1-theta[i,1])^2
}

library(doParallel)

cl <- makeCluster(nbcores) 
registerDoParallel(cl)   

train <- foreach(i=1:N, .combine="rbind") %dopar% 
{
  s[i,] <- IBDSim_wrapper_IBD(g=theta[i,1],m=theta[i,2],habitatSize=theta[i,3])
}
stopCluster(cl)    


train <- cbind(theta,train)
colnames(train)[1:4] <- c("g","m","habitatsize","sigma2")
head(train) 
train <- as.data.frame(train)
train1 <- train[,-c(2:4)]
train2 <- train[,-c(1,3:4)]
train3 <- train[,-c(1:2,4)]
train4 <- train[,-c(1:3)]


len <- length(train)

train_without_Q <- train[,-c(17:len)]
train1_without_Q <- train_without_Q[,-c(2:4)]
train2_without_Q <- train_without_Q[,-c(1,3:4)]
train3_without_Q <- train_without_Q[,-c(1:2,4)]
train4_without_Q <- train_without_Q[,-c(1:3)]


train_without_ar_er <- train[,-c(13:16)]
train1_without_ar_er <- train_without_ar_er[,-c(2:4)]
train2_without_ar_er <- train_without_ar_er[,-c(1,3:4)]
train3_without_ar_er <- train_without_ar_er[,-c(1:2,4)]
train4_without_ar_er <- train_without_ar_er[,-c(1:3)]

train_without_ar_er_Q <- train[,-c(13:len)]
train1_without_ar_er_Q <- train_without_ar_er_Q[,-c(2:4)]
train2_without_ar_er_Q <- train_without_ar_er_Q[,-c(1,3:4)]
train3_without_ar_er_Q <- train_without_ar_er_Q[,-c(1:2,4)]
train4_without_ar_er_Q <- train_without_ar_er_Q[,-c(1:3)]


###########################    TABLE TEST PARALLELISATION   ###########################    

Ntest <- 100

theta.test <- matrix(c(theta1_obs,theta2_obs, theta3_obs, theta4_obs),Ntest,4, byrow=TRUE)
s.test <- matrix(0,Ntest,nb_sum_stats)

library(doParallel)

cl <- makeCluster(nbcores) 
registerDoParallel(cl)   
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

#### Forêts ####

# Pour le paramètre g
model1 <- regAbcrf(g~.,data=train1,ntree=500, paral=para)
model1_without_Q <- regAbcrf(g~.,data=train1_without_Q,ntree=500, paral=para)
model1_without_ar_er <- regAbcrf(g~.,data=train1_without_ar_er,ntree=500, paral=para)
model1_without_ar_er_Q <- regAbcrf(g~.,data=train1_without_ar_er_Q,ntree=500, paral=para)

# Pour le paramètre m
model2 <- regAbcrf(m~.,data=train2,ntree=500, paral=para)
model2_without_Q <- regAbcrf(m~.,data=train2_without_Q,ntree=500, paral=para)
model2_without_ar_er <- regAbcrf(m~.,data=train2_without_ar_er,ntree=500, paral=para)
model2_without_ar_er_Q <- regAbcrf(m~.,data=train2_without_ar_er_Q,ntree=500, paral=para)

# Pour habitatsize
model3 <- regAbcrf(habitatsize~.,data=train3,ntree=500, paral=para)
model3_without_Q <- regAbcrf(habitatsize~.,data=train3_without_Q,ntree=500, paral=para)
model3_without_ar_er <- regAbcrf(habitatsize~.,data=train3_without_ar_er,ntree=500, paral=para)
model3_without_ar_er_Q <- regAbcrf(habitatsize~.,data=train3_without_ar_er_Q,ntree=500, paral=para)

# Pour sigma2
model4 <- regAbcrf(sigma2~.,data=train4,ntree=500, paral=para)
model4_without_Q <- regAbcrf(sigma2~.,data=train4_without_Q,ntree=500, paral=para)
model4_without_ar_er <- regAbcrf(sigma2~.,data=train4_without_ar_er,ntree=500, paral=para)
model4_without_ar_er_Q <- regAbcrf(sigma2~.,data=train4_without_ar_er_Q,ntree=500, paral=para)

plot(model1, main="Estimation du paramètre g")
plot(model2, main="Estimation du paramètre m")
plot(model3, main="Estimation du paramètre habitatsize")
plot(model4, main="Estimation du paramètre sigma2")

#### Predict ####

pred.theta1 <- predict(model1,test,train1)
pred.theta1_without_Q <- predict(model1_without_Q,test_without_Q,train1_without_Q)
pred.theta1_without_ar_er <- predict(model1_without_ar_er,test_without_ar_er,train1_without_ar_er)
pred.theta1_without_ar_er_Q <- predict(model1_without_ar_er_Q,test_without_ar_er_Q,train1_without_ar_er_Q)

pred.theta2 <- predict(model2,test,train2)
pred.theta2_without_Q <- predict(model2_without_Q,test_without_Q,train2_without_Q)
pred.theta2_without_ar_er <- predict(model2_without_ar_er,test_without_ar_er,train2_without_ar_er)
pred.theta2_without_ar_er_Q <- predict(model2_without_ar_er_Q,test_without_ar_er_Q,train2_without_ar_er_Q)

pred.theta3 <- predict(model3,test,train3)
pred.theta3_without_Q <- predict(model3_without_Q,test_without_Q,train3_without_Q)
pred.theta3_without_ar_er <- predict(model3_without_ar_er,test_without_ar_er,train3_without_ar_er)
pred.theta3_without_ar_er_Q <- predict(model3_without_ar_er_Q,test_without_ar_er_Q,train3_without_ar_er_Q)

pred.theta4 <- predict(model4,test[],train4)
pred.theta4_without_Q <- predict(model4_without_Q,test_without_Q,train4_without_Q)
pred.theta4_without_ar_er <- predict(model4_without_ar_er,test_without_ar_er,train4_without_ar_er)
pred.theta4_without_ar_er_Q <- predict(model4_without_ar_er_Q,test_without_ar_er_Q,train4_without_ar_er_Q)


#### Plot Predict ####

plot(theta.test[,1],pred.theta1$med)
abline(a=0,b=1,col="red")

plot(theta.test[,2],pred.theta2$med)
abline(a=0,b=1,col="red")

plot(theta.test[,3],pred.theta3$med)
abline(a=0,b=1,col="red")

plot(theta.test[,4],pred.theta4$med)
abline(a=0,b=1,col="red")

#### Plot densités a posteriori ####
densityPlot(model1,test,train1)
abline(v=theta1_obs,col="red", lwd=2)

densityPlot(model2,test,train2)
abline(v=theta2_obs,col="red", lwd=2)

densityPlot(model3,test,train3)
abline(v=theta3_obs,col="red", lwd=2)

densityPlot(model4,test,train4, xlim=c(0,50000))
abline(v=theta4_obs,col="red", lwd=2)


#### MSE ####

MSE1 <- mean((pred.theta1$expectation - theta.test[,1])^2)
MSE1_without_Q <- mean((theta.test[,1] - pred.theta1_without_Q$expectation)^2)
MSE1_without_Q
MSE1_without_ar_er <- mean((theta.test[,1] - pred.theta1_without_ar_er$expectation)^2)
MSE1_without_ar_er
MSE1_without_ar_er_Q <- mean((theta.test[,1] - pred.theta1_without_ar_er_Q$expectation)^2)
MSE1_without_ar_er_Q

MSE1_med <- mean((pred.theta1$med - theta.test[,1])^2)
MSE1_med_without_Q <- mean((theta.test[,1] - pred.theta1_without_Q$med)^2)
MSE1_med_without_Q
MSE1_med_without_ar_er <- mean((theta.test[,1] - pred.theta1_without_ar_er$med)^2)
MSE1_med_without_ar_er
MSE1_med_without_ar_er_Q <- mean((theta.test[,1] - pred.theta1_without_ar_er_Q$med)^2)
MSE1_med_without_ar_er_Q




MSE2 <- mean((theta.test[,2] - pred.theta2$expectation)^2)
MSE2 
MSE2_without_Q <- mean((theta.test[,2] - pred.theta2_without_Q$expectation)^2)
MSE2_without_Q
MSE2_without_ar_er <- mean((theta.test[,2] - pred.theta2_without_ar_er$expectation)^2)
MSE2_without_ar_er
MSE2_without_ar_er_Q <- mean((theta.test[,2] - pred.theta2_without_ar_er_Q$expectation)^2)
MSE2_without_ar_er_Q

MSE2_med <- mean((theta.test[,2] - pred.theta2$med)^2)
MSE2_med 
MSE2_med_without_Q <- mean((theta.test[,2] - pred.theta2_without_Q$med)^2)
MSE2_med_without_Q
MSE2_med_without_ar_er <- mean((theta.test[,2] - pred.theta2_without_ar_er$med)^2)
MSE2_med_without_ar_er
MSE2_med_without_ar_er_Q <- mean((theta.test[,2] - pred.theta2_without_ar_er_Q$med)^2)
MSE2_med_without_ar_er_Q

MSE3 <- mean((theta.test[,3] - pred.theta3$expectation)^2)
MSE3 
MSE3_without_Q <- mean((theta.test[,3] - pred.theta3_without_Q$expectation)^2)
MSE3_without_Q
MSE3_without_ar_er <- mean((theta.test[,3] - pred.theta3_without_ar_er$expectation)^2)
MSE3_without_ar_er
MSE3_without_ar_er_Q <- mean((theta.test[,3] - pred.theta3_without_ar_er_Q$expectation)^2)
MSE3_without_ar_er_Q

MSE3_med <- mean((theta.test[,3] - pred.theta3$med)^2)
MSE3_med 
MSE3_med_without_Q <- mean((theta.test[,3] - pred.theta3_without_Q$med)^2)
MSE3_med_without_Q
MSE3_med_without_ar_er <- mean((theta.test[,3] - pred.theta3_without_ar_er$med)^2)
MSE3_med_without_ar_er
MSE3_med_without_ar_er_Q <- mean((theta.test[,3] - pred.theta3_without_ar_er_Q$med)^2)
MSE3_med_without_ar_er_Q

MSE3_relatif_med <- mean((theta.test[,3] - pred.theta3$med)^2)/theta3_obs^2
MSE3_relatif_med 
MSE3_relatif_med_without_Q <- mean((theta.test[,3] - pred.theta3_without_Q$med)^2)/theta3_obs^2
MSE3_relatif_med_without_Q
MSE3_relatif_med_without_ar_er <- mean((theta.test[,3] - pred.theta3_without_ar_er$med)^2)/theta3_obs^2
MSE3_relatif_med_without_ar_er
MSE3_relatif_med_without_ar_er_Q <- mean((theta.test[,3] - pred.theta3_without_ar_er_Q$med)^2)/theta3_obs^2
MSE3_relatif_med_without_ar_er_Q

MSE4 <- mean((theta.test[,4] - pred.theta4$expectation)^2)
MSE4 
MSE4_without_Q <- mean((theta.test[,4] - pred.theta4_without_Q$expectation)^2)
MSE4_without_Q
MSE4_without_ar_er <- mean((theta.test[,4] - pred.theta4_without_ar_er$expectation)^2)
MSE4_without_ar_er
MSE4_without_ar_er_Q <- mean((theta.test[,4] - pred.theta4_without_ar_er_Q$expectation)^2)
MSE4_without_ar_er_Q

MSE4_med <- mean((theta.test[,4] - pred.theta4$med)^2)
MSE4_med 
MSE4_med_without_Q <- mean((theta.test[,4] - pred.theta4_without_Q$med)^2)
MSE4_med_without_Q
MSE4_med_without_ar_er <- mean((theta.test[,4] - pred.theta4_without_ar_er$med)^2)
MSE4_med_without_ar_er
MSE4_med_without_ar_er_Q <- mean((theta.test[,4] - pred.theta4_without_ar_er_Q$med)^2)
MSE4_med_without_ar_er_Q

MSE4_med_relatif <- mean((theta.test[,4] - pred.theta4$med)^2)/theta4_obs^2
MSE4_med_relatif 
MSE4_med_relatif_without_Q <- mean((theta.test[,4] - pred.theta4_without_Q$med)^2)/theta4_obs^2
MSE4_med_relatif_without_Q
MSE4_med_relatif_without_ar_er <- mean((theta.test[,4] - pred.theta4_without_ar_er$med)^2)/theta4_obs^2
MSE4_med_relatif_without_ar_er
MSE4_med_relatif_without_ar_er_Q <- mean((theta.test[,4] - pred.theta4_without_ar_er_Q$med)^2)/theta4_obs^2
MSE4_med_relatif_without_ar_er_Q

########################    Biais    ########################
### attention changement: biais relatif
biais1<- mean(pred.theta1$expectation) - theta1_obs
biais1_without_Q <- mean(pred.theta1_without_Q$expectation) - theta1_obs
biais1_without_ar_er <- mean(pred.theta1_without_ar_er$expectation) - theta1_obs
biais1_without_ar_er_Q <- mean(pred.theta1_without_ar_er_Q$expectation) - theta1_obs

biais1_relatif <- (mean(pred.theta1$expectation) - theta1_obs)/theta1_obs
biais1_relatif_without_Q <- (mean(pred.theta1_without_Q$expectation) - theta1_obs)/theta1_obs
biais1_relatif_without_ar_er <- (mean(pred.theta1_without_ar_er$expectation) - theta1_obs)/theta1_obs
biais1_relatif_without_ar_er_Q <- (mean(pred.theta1_without_ar_er_Q$expectation) - theta1_obs)/theta1_obs

biais1_med <- mean(pred.theta1$med) - theta1_obs
biais1_med_without_Q <- mean(pred.theta1_without_Q$med) - theta1_obs
biais1_med_without_ar_er <- mean(pred.theta1_without_ar_er$med) - theta1_obs
biais1_med_without_ar_er_Q <- mean(pred.theta1_without_ar_er_Q$med) - theta1_obs

biais1_med_relatif <- (mean(pred.theta1$med) - theta1_obs)/theta1_obs
biais1_med_relatif_without_Q <- (mean(pred.theta1_without_Q$med) - theta1_obs)/theta1_obs
biais1_med_relatif_without_ar_er <- (mean(pred.theta1_without_ar_er$med) - theta1_obs)/theta1_obs
biais1_med_relatif_without_ar_er_Q <- (mean(pred.theta1_without_ar_er_Q$med) - theta1_obs)/theta1_obs




biais2 <- mean(pred.theta2$expectation) - theta2_obs
biais2_without_Q <- mean(pred.theta2_without_Q$expectation) - theta2_obs
biais2_without_ar_er <- mean(pred.theta2_without_ar_er$expectation) - theta2_obs
biais2_without_ar_er_Q <- mean(pred.theta2_without_ar_er_Q$expectation) - theta2_obs

biais2_med <- mean(pred.theta2$med) - theta2_obs
biais2_med_without_Q <- mean(pred.theta2_without_Q$med) - theta2_obs
biais2_med_without_ar_er <- mean(pred.theta2_without_ar_er$med) - theta2_obs
biais2_med_without_ar_er_Q <- mean(pred.theta2_without_ar_er_Q$med) - theta2_obs


biais3 <- mean(pred.theta3$expectation) - theta3_obs
biais3_without_Q <- mean(pred.theta3_without_Q$expectation) - theta3_obs
biais3_without_ar_er <- mean(pred.theta3_without_ar_er$expectation) - theta3_obs
biais3_without_ar_er_Q <- mean(pred.theta3_without_ar_er_Q$expectation) - theta3_obs

biais3_med <- mean(pred.theta3$med) - theta3_obs
biais3_med_without_Q <- mean(pred.theta3_without_Q$med) - theta3_obs
biais3_med_without_ar_er <- mean(pred.theta3_without_ar_er$med) - theta3_obs
biais3_med_without_ar_er_Q <- mean(pred.theta3_without_ar_er_Q$med) - theta3_obs

biais4 <- mean(pred.theta4$expectation) - theta4_obs
biais4_without_Q <- mean(pred.theta4_without_Q$expectation) - theta4_obs
biais4_without_ar_er <- mean(pred.theta4_without_ar_er$expectation) - theta4_obs
biais4_without_ar_er_Q <- mean(pred.theta4_without_ar_er_Q$expectation) - theta4_obs

biais4_med <- mean(pred.theta4$med) - theta4_obs
biais4_med_without_Q <- mean(pred.theta4_without_Q$med) - theta4_obs
biais4_med_without_ar_er <- mean(pred.theta4_without_ar_er$med) - theta4_obs
biais4_med_without_ar_er_Q <- mean(pred.theta4_without_ar_er_Q$med) - theta4_obs

########################    NMAE    ########################
for (i in 1:4) {
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

NMAE1_med <- mean(abs((train1[1:100,1]-pred.theta1$med)/train1[1:100,1]))


NMAE1 <- mean(abs((theta1_obs-pred.theta1$expectation)/theta1_obs))
NMAE1_med <- mean(abs((theta1_obs-pred.theta1$med)/theta1_obs))


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
                  c(mean(pred.theta1$med), mean(pred.theta1_without_Q$med), mean(pred.theta1_without_ar_er$med), mean(pred.theta1_without_ar_er_Q$med)),
                  c(biais1, biais1_without_Q, biais1_without_ar_er, biais1_without_ar_er_Q),
                  c(var(pred.theta1$expectation), var(pred.theta1_without_Q$expectation), var(pred.theta1_without_ar_er$expectation), var(pred.theta1_without_ar_er_Q$expectation)),
                  c(MSE1, MSE1_without_Q, MSE1_without_ar_er,MSE1_without_ar_er_Q),                  
                  c(NMAE1, NMAE1_without_Q, NMAE1_without_ar_er,NMAE1_without_ar_er_Q),
                  OOB_results1)
rownames(results1) <- c("Valeur moyenne", "Valeur médiane","Biais","Variance", "MSE","NMAE", "OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")

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

results4 <- rbind(c(mean(pred.theta4$expectation), mean(pred.theta4_without_Q$expectation), mean(pred.theta4_without_ar_er$expectation), mean(pred.theta4_without_ar_er_Q$expectation)),
                  c(mean(pred.theta4$med), mean(pred.theta4_without_Q$med), mean(pred.theta4_without_ar_er$med), mean(pred.theta4_without_ar_er_Q$med)),
                  c(biais4, biais4_without_Q, biais4_without_ar_er, biais4_without_ar_er_Q),
                  c(var(pred.theta4$expectation), var(pred.theta4_without_Q$expectation), var(pred.theta4_without_ar_er$expectation), var(pred.theta4_without_ar_er_Q$expectation)),
                  c(MSE4, MSE4_without_Q, MSE4_without_ar_er,MSE4_without_ar_er_Q),                  
                  c(NMAE4, NMAE4_without_Q, NMAE4_without_ar_er,NMAE4_without_ar_er_Q),
                  OOB_results4)
rownames(results4) <- c("Valeur moyenne","Valeur médiane","Biais","Variance", "MSE","NMAE", "OOB MSE", "OOB NMAE", "OOB med MSE", "OOB med NMAE")


######################  Sauvegardes Rdata  ###################### 
dir_name_1 <- paste(gsub(" ", "_", deb),sep="")
dir_cr<- dir.create(dir_name_1)
write.csv(results1, file=paste(dir_name_1,"/Results1_100_meme_para_gmh.csv", sep=""))
write.csv(results2, file=paste(dir_name_1,"/Results2_100_meme_para_gmh.csv", sep=""))
write.csv(results3, file=paste(dir_name_1,"/Results3_100_meme_para_gmh.csv", sep=""))

save.image(file=paste("Rdata/",dir_name_1,".Rdata", sep=""))


# write.csv(results1, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results1_100_meme_para_gmh.csv", sep=""))
# write.csv(results2, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results2_100_meme_para_gmh.csv", sep=""))
# write.csv(results3, file=paste("/home/vernierc/Documents/Logiciels/ABCRF/",gsub(" ", "_", deb),"Results3_100_meme_para_gmh.csv", sep=""))
# 
# save.image(file=paste("/home/vernierc/Documents/Logiciels/ABCRF/Rdata/",gsub(" ", "_", deb),".Rdata", sep=""))
print(deb)
fin <- Sys.time()
duree <- fin - deb
#calcul de l'OOB error




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