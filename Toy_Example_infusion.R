########################################
#             INNFUSION
########################################

library(Infusion)


# Test estimation paramètres gaussienne
myrnorm <-function(mu, s2, sample.size) 
  {
  s <- rnorm(n=sample.size, mean=mu, sd=sqrt(s2))
  return(c(mean=mean(s), var=var(s)))
}

myrnorm(0,1,5)
myrnorm(0,1,10)
myrnorm(0,1,100)


S <- myrnorm(0,1,1000)
# Population variance = Sample_variance*(n-1)/n
S[2]*999/1000

set.seed(123) ## initialize the random generator
mu <- 4
s2 <- 1
n <- 40
ssize <- n

Sobs <- myrnorm(mu=mu, s2=s2, sample.size=n)
Sobs

# Population variance = Sample_variance*(n-1)/n
s2_ech <- Sobs[2]*(n-1)/n


# Confidence interval for mu (student)
#######################################################
error <- qt(0.975,df=n-1)*sqrt(Sobs[2])/sqrt(n)  #sample mean et sample sd
left <- Sobs[1] - error
right <- Sobs[1] + error
CI <- c(left, right)
CI <- matrix(CI, ncol=2,byrow=TRUE)
colnames(CI) = c("5%","95%")
CI 
# :)

logLnorm <- function(mu,s2,Sobs) { dnorm(Sobs[1],mean=mu,sd=sqrt(s2/ssize),log=TRUE)+
    dchisq((ssize-1)*Sobs[2]/s2,df = (ssize-1),log=TRUE)+log((ssize-1)/s2) }


# Confidence interval mu loi normale
#######################################################

error2 <- qnorm(0.975)*sqrt(Sobs[2])/sqrt(n)  #sample mean et sample sd
left2 <- Sobs[1] - error2
right2 <- Sobs[1] + error2
CI2 <- c(left2, right2)
CI2 <- matrix(CI2, ncol=2,byrow=TRUE)
colnames(CI2) = c("5%","95%")
CI2


# IC de la variance
####################################################
set.seed(123)
s <- rnorm(n=ssize, mean=mu, sd=sqrt(s2))
inf <-((n-1)*Sobs[2])/qchisq(p=0.975, df=n-1)
inf
sup <-((n-1)*Sobs[2])/qchisq(p=0.025, df=n-1)
sup
CI_var <- c(inf,sup)
CI_var


#####################################################
# init_grid pour générer un premier ensemble de points de paramètres (grille irrégulière)
# generating an irregular grid of parameter values, with duplicate values for some parameters.


options(digits=3)
library(Infusion)
parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                   upper=c(mu=5.2,s2=3,sample.size=40))
parsp
# Each row defines a list of arguments of vector of the function simulating the summary statistics.
# première colonne?


#  build a list of simulated distributions from this set of parameter values:
simuls <- add_simulation(NULL,Simulate="myrnorm",par.grid=parsp,verbose=FALSE)
# For each parameter point, add_simulation has simulated an empirical distri-
# bution (by default, of 1000 realizations). Here the simulation function myrnorm
# is available in the R session to produce samples from the distribution of the
# summary statistics. In more involved applications the simulation code may
# not be callable from R. This case is also handled by add_simulation, using its
# newsimuls argument.

# We estimate the probability density (“likelihood”) of the observed summary
# statistics for each simulated distribution:
densv <- infer_logLs(simuls,stat.obs=Sobs,verbose=FALSE)
# infer_logLs infers a probability density by smoothing the empirical distribu-
#   tion of the simulated statistics for each parameter value

# From all estimated likelihoods, we estimate a summary likelihood surface
# by smoothing the density estimates:
slik <- infer_surface(densv) ## infer a log-likelihood surface
slik## using REML to infer the S-likelihood surface...
slik <-MSL(slik)
confint(slik,"mu")
plot(slik,filled=TRUE)

maxit <- 2 ## maximum number of 'refine' iterations for all examples
slik2 <- refine(slik,filled=TRUE,verbose=FALSE,maxit=maxit)
slik2 <- MSL(slik2)
slik2
plot(slik2,filled=TRUE)

plot1Dprof(slik)
plot2Dprof(slik)

####################################################
# ?
blurred <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  s <- exp(s/4)
  return(c(mean=mean(s),var=var(s)))
}

set.seed(123)
dSobs <- blurred(mu=4,s2=1,sample.size=40) ## stands for the actual data to be analyzed
dsimuls <- add_simulation(NULL,Simulate="blurred",par.grid=parsp,verbose=FALSE)


###################################################
# Krigreage

# projecteurs de mu et sigma
mufit <- project("mu",stats=c("mean","var"),data=dsimuls)
s2fit <- project("s2",stats=c("mean","var"),data=dsimuls)

# représentation graphique 
mapMM(mufit,map.asp=1,plot.title=title(main="prediction of mean parameter",
                                       xlab="blurred mean",ylab="blurred var"))
mapMM(s2fit,map.asp=1,plot.title=title(main="prediction of variance parameter",
                                       xlab="blurred mean",ylab="blurred var"))

# valeurs des projections des stats résumées (projSobs) et de la table de simulation (projSimuls)
projSobs <- project(dSobs,projectors=list(MEAN=mufit,VAR=s2fit))
projSimuls <- project(dsimuls,projectors=list(MEAN=mufit,VAR=s2fit)) #?
cdensv <- infer_logLs(projSimuls,stat.obs=projSobs,verbose=FALSE)
cslik <- infer_surface(cdensv)
cslik <- MSL(cslik)
cslik <- refine(cslik,maxit=maxit,verbose=FALSE)
cslik



# Neural network
mufit <- project("mu",stats=c("mean","var"),data=dsimuls,method="neuralNet")
s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,method="neuralNet")
projSobs <- project(dSobs,projectors=list("MEAN"=mufit,"VAR"=s2fit))
projSimuls <- project(dsimuls,projectors=list("MEAN"=mufit,"VAR"=s2fit))
cdensv <- infer_logLs(projSimuls,stat.obs=projSobs,verbose=FALSE)
nnslik <- infer_surface(cdensv)
nnslik <- MSL(nnslik,verbose=FALSE)
nnslik <- refine(nnslik,maxit=maxit,verbose=FALSE)

