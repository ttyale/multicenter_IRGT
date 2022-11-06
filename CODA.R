#####################################################
###   This file reproduces Table 5 ( CODA trial ) ###
#####################################################

#function calculating number of centers
sample_size<-function(eff,t1,t2,n1,n2,rho1,rho2,sigma1,sigma2){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #t2: number of therapists in control arm
  #n1: cluster size in treatment arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #sigma1: total variance for treatment arm
  #sigma2: total variance for control arm
  
  alpha=0.05         #alpha: type I error rate             
  beta=0.1           #beta: type II error rate  
  alpha1 <- rho2*sigma2/sigma1       #alpha1: between-therapist ICC for treatment arm
  alpha2=rho2                        #alpha1: between-therapist ICC for control arm
  vb<-sigma1/n1/t1*((1+(n1-1)*rho1)+n1*(t1-1)*alpha1)+sigma2/n2/t2*((1+(n2-1)*rho2)+n2*(t2-1)*alpha2)-2*sqrt(sigma1*sigma2)*sqrt(alpha1*alpha2)
  I <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))                    #I: number of centers   
  percentage <- alpha1/rho1                                                        #percentage: ratio of between-therapist ICC versus within-therapist ICC
  return(c(rho1,rho2,alpha1,alpha2,t1,n1,t2,n2,percentage,I))
}

#ICCs configuration for CODA trial
rho2<- c(rep(0.01,4),rep(0.02,4),rep(0.04,3),rep(0.08,2))                  #rho2: within-therapist ICC in control arm
rho1 <- c(0.02,0.05,0.1,0.2,0.04,0.05,0.10,0.20,0.08,0.10,0.20,0.16,0.20)  #rho1: within-therapist ICC in treatment arm
T1 <- 48                            #T1: number of therapists in treatment arm
T2 <- 1                             #T2: number of therapists in control arm
N1 <- 16.16667                      #N1: cluster size in treatment arm
N2 <- 776                           #N2: cluster size in control arm
sigma1 <- 0.0169                    #sigma1: total variance for treatment arm
sigma2 <- c(0.10^2,0.11^2,0.12^2,0.13^2,0.14^2)     #sigma2: total variance for control arm
eff <- 0.008                                        #eff: effect size

fidata <- c()                      #fidata: storing simulation results
for (k in 1:length(sigma2)){
  for (i in 1:length(rho1)){
    fidata <- rbind(fidata,sample_size(eff=eff,t1=T1,t2=T2,n1=N1,n2=N2,
                                       rho1=rho1[i],rho2=rho2[i],sigma1=sigma1,sigma2=sigma2[k]))
  }
}

colnames(fidata) <- c("rho1","rho2","alpha1","alpha2","T1","N1","T2","N2","Percentage","I")
fidata #Table 5 output


