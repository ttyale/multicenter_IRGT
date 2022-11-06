#####################################################
###   This file reproduces Table 4 ( STPD trial ) ###
#####################################################

#function calculating number of centers
sample_size<-function(eff,t1,t2,n1,n2,rho1,rho2,alpha1,sigma1,sigma2){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha1: between-therapist ICC for treatment arm
  #sigma1: total variance for treatment arm
  #sigma2: total variance for control arm
  
  alpha=0.05        #alpha: type I error rate       
  beta=0.1          #beta: type II error rate  
  
  alpha2=alpha1*sigma1/sigma2            #alpha2: between-therapist ICC for control arm
  vb <- sigma1/n1/t1*((1+(n1-1)*rho1)+n1*(t1-1)*alpha1)+sigma2/n2/t2*((1+(n2-1)*rho2)+n2*(t2-1)*alpha2)-2*sqrt(sigma1*sigma2)*sqrt(alpha1*alpha2)  
  I <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))                   #I: number of centers
  percentage <- alpha1/rho1                                                       #percentage: ratio of between-therapist ICC versus within-therapist ICC
  return(c(rho1,rho2,alpha1,alpha2,t1,n1,t2,n2,percentage,I))
}

#function calculating realization of objective loss function
obf <- function(J,t1,t2,n1,n2,rho1,rho2,alpha1,sigma1,sigma2){
  #J: fixed total number of therapists for two arms
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha1: between-therapist ICC for treatment arm
  #sigma1: total variance for treatment arm
  #sigma2: total variance for control arm
  
  alpha2=alpha1*sigma1/sigma2            #alpha2: between-therapist ICC for control arm
  f <- sigma1*(1+(n1-1)*rho1+n1*(t1-1)*alpha1)/t1/n1+sigma2*(1+(n2-1)*rho2+n2*(J-t1-1)*alpha2)/(J-t1)/n2  #f: realization of loss function
  return(f)
}

#function calculating optimal number of therapists
opthera <- function(eff,J,n1,n2,rho1,rho2,alpha1,sigma1,sigma2){
  #eff: effect size 
  #J: fixed total number of therapists for two arms
  #n1: cluster size in treatment arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha1: between-therapist ICC for treatment arm
  #sigma1: total variance for treatment arm
  #sigma2: total variance for control arm
  
  alpha2=alpha1*sigma1/sigma2            #alpha2: between-therapist ICC for control arm
  T1op <- J*sqrt(sigma1)*sqrt(n2*(1-rho1+n1*(rho1-alpha1)))/                                           #T1op: non-integer optimal number of therapist for treatment arm
    ((sqrt(sigma1)*sqrt(n2*(1-rho1+n1*(rho1-alpha1))))+(sqrt(sigma2*n1*(1-rho2+n2*(rho2-alpha2))))) 
  T1op1 <- ceiling(T1op)                                                                               #T1op1: 1st candidate of optimal number of therapist for treatment arm after rounding up
  T1op2 <- T1op1 - 1                                                                                   #T1op2: 2nd candidate of optimal number of therapist for treatment arm after rounding down
  obf1 <- obf(J,T1op1,(J-T1op1),n1,n2,rho1,rho2,alpha1,sigma1,sigma2)                                  #obf1: realization of loss function by inserting 1st candidate
  obf2 <- obf(J,T1op2,(J-T1op2),n1,n2,rho1,rho2,alpha1,sigma1,sigma2)                                  #obf2: realization of loss function by inserting 2nd candidate
  if (obf1 >= obf2){
    T1op <- T1op2
  } else {T1op <- T1op1}
  T2op <- J-T1op                                                                                       #T2op: optimal number of therapist for the control arm
  opI <- sample_size(eff=eff,t1=T1op,t2=T2op,n1=n1,n2=n2,                                              #opI: optimal number of sites
                     rho1,rho2,alpha1,sigma1,sigma2)[10]
  return(c(T1op,T2op,opI))
}

#ICCs configuration for STPD trial
alpha1 <- 0.007                                                  #alpha1: between-therapist ICC for treatment arm
rho1<- c(rep(0.01,4),rep(0.05,4),rep(0.10,4),rep(0.20,4))        #rho1: within-therapist ICC for treatment arm
rho2 <- rep(c(0.01,0.05,0.10,0.20),4)                            #rho2: within-therapist ICC for control arm
T1 <- 66                        #T1: number of therapists in treatment arm
T2 <- 79                        #T2: number of therapists in control arm
N1 <- 2.2                       #N1: cluster size in treatment arm
N2 <- 1.7                       #N2: cluster size in control arm
sigma1 <- 169                   #sigma1: total variance for treatment arm
sigma2 <- c(11^2,13^2,14^2)     #sigma2: total variance for control arm
eff <- 1.3                      #eff: effect size

fidata <- c()                   #fidata: storing simulation results
for (k in 1:length(sigma2)){
  for (i in 1:length(rho1)){
    fidata <- rbind(fidata,sample_size(eff=eff,t1=T1,t2=T2,n1=N1,n2=N2,
                                       rho1=rho1[i],rho2=rho2[i],alpha1=alpha1,sigma1=sigma1,sigma2=sigma2[k]))
  }
}

J <- T1+T2                      #J: total number of therapists for two arms 
opdata <- c()                   #opdata: storing simulation results for optimal number of therapists allocation and the corresponding number of centers
for (k in 1:length(sigma2)){
  for (i in 1:length(rho1)){
      opdata <- rbind(opdata,opthera(eff=eff,J=J,n1=N1,n2=N2,
                                     rho1=rho1[i],rho2=rho2[i],alpha1=alpha1,sigma1=sigma1,sigma2=sigma2[k]))
  }
}

#name the columns of resulting data sets
colnames(opdata) <- c("optimalT1","optimalT2","optimal n.site")
opdata
colnames(fidata) <- c("rho1","rho2","alpha1","alpha2","T1","N1","T2","N2","Percentage","I")
fidata

