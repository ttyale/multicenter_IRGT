###############################################################################
###   This file reproduces simulation results in Table 2 and Web Table 1-5 ###
###############################################################################

require("nlme")
#Sample size calculation based on between-therapist ICCs and within-therapist ICCs
sample_size<-function(eff,t1,t2,n1,n2,rho1,alpha1,rho2,alpha2,sigma1){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #alpha1: between-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha2: between-therapist ICC for control arm
  #sigma1: total variance for treatment arm
  
  alpha=0.05        #alpha: type I error rate         
  beta=0.2          #beta: type II error rate  
  if (alpha2==0){sigma2=1}              #fix total variance for control arm if alpha2=0
  else{sigma2=sigma1*alpha1/alpha2}     #sigma2 is determined by sigma1, alpha1, alpha2
  vb <- sigma1/n1/t1*((1+(n1-1)*rho1)+n1*(t1-1)*alpha1)+sigma2/n2/t2*((1+(n2-1)*rho2)+n2*(t2-1)*alpha2)-2*sqrt(sigma1*sigma2)*sqrt(alpha1*alpha2)  
  I <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))  #I: number of centers
  return(I)
}

#predicted power calculation by formula
power<-function(I,eff,t1,t2,n1,n2,rho1,alpha1,rho2,alpha2,sigma1){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #alpha1: between-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha2: between-therapist ICC for control arm
  #sigma1: total variance for treatment arm
  
  alpha=0.05                                 #type I error rate 
  if (alpha2==0){sigma2=1}                   #fix total variance for control arm if alpha2=0
  else{sigma2=sigma1*alpha1/alpha2}          #sigma2 is determined by sigma1, alpha1, alpha2
  vb<-sigma1/n1/t1/I*((1+(n1-1)*rho1)+n1*(t1-1)*alpha1)+sigma2/n2/t2/I*((1+(n2-1)*rho2)+n2*(t2-1)*alpha2)-2*sqrt(sigma1*sigma2)*sqrt(alpha1*alpha2)/I
  pw<-pnorm(eff/sqrt(vb)-qnorm(1-alpha/2))   #pw: predicted power
  return(pw)
}

#############################################################################
#simulation for empirical power

#data generation process for a given scenario with parameter configuration
data_gen<-function(I=5,eff=0.2,t1=2,t2=3,n1=10,n2=13,rho1=0.05,alpha1=0.03,rho2=0.05,alpha2=0.02,sigma1=1){
  #I: number of centers
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #alpha1: between-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha2: between-therapist ICC for control arm
  #sigma1: total variance for treatment arm
  
  if (alpha2==0){sigma2=1}                  #fix total variance for control arm if alpha2=0
  else{sigma2=sigma1*alpha1/alpha2}         #sigma2 is determined by sigma1, alpha1, alpha2
  intercept=0.5                             #intercept: beta0
  N=I*(t1*n1+t2*n2)                         #N: total number of patients for two arms
  patid <- 1:N                              #patid: patient id 
  trt<-rep(c(rep(1,t1*n1),rep(0,t2*n2)),I)  #trt: treatment indicator X_ijk
  sigma_lambda<-sigma1*alpha1               #sigma_lambda: center-level random effect variance
  sigma_theta1<-rho1*sigma1-sigma_lambda    #sigma_theta1: therapist random effect variance in treatment arm
  sigma_theta2<-rho2*sigma2-sigma_lambda    #sigma_theta2: therapist random effect variance in control arm
  
  sigma_epsilon1<-sigma1-sigma_lambda-sigma_theta1  #sigma_epsilon1: individual residual variance in treatment arm
  sigma_epsilon2<-sigma2-sigma_lambda-sigma_theta2  #sigma_epsilon2: individual residual variance in control arm
  
  strata<-rep(1:I,each=(t1*n1+t2*n2))                         #strata: center id
  tid_s<-c(rep(1:t1,each=n1),rep((t1+1):(t1+t2),each=n2))     #tid_s: therapist id for each center
  tid<-c()                                                    #tid: overall therapist id 
  for(i in 1:I){
    tid<-c(tid,((i-1)*(t1+t2)+tid_s))
  }
  
  I_effect<-rep(rnorm(I,0,sqrt(sigma_lambda)),each=(t1*n1+t2*n2))  #I_effect: center random effect
  
  therapist_effect1<-rnorm(I*(t1+t2),0,sqrt(sigma_theta1))         #therapist_effect1: therapist random effect for treatment arm
  therapist_effect2<-rnorm(I*(t1+t2),0,sqrt(sigma_theta2))         #therapist_effect2: therapist random effect for control arm
  
  therapist_effect<-rep(therapist_effect1,as.numeric(table(tid)))*trt+rep(therapist_effect2,as.numeric(table(tid)))*(1-trt)  #therapist_effect: combined therapist random effect
  
  epsilon<-c()        #epsilon: combined measurement error for two arms
  for (i in 1:I){
    epsilon_s <-c(rnorm(t1*n1,0,sqrt(sigma_epsilon1)),rnorm(t2*n2,0,sqrt(sigma_epsilon2)))
    epsilon<-c(epsilon,epsilon_s)
  }
  
  y=intercept-eff*trt+I_effect+therapist_effect+epsilon  #y: outcome under alternative
  y0=intercept+I_effect+therapist_effect+epsilon         #y0: outcome under null
  
  df<-data.frame(y=y,y0=y0,trt=trt,tid=tid,sid=strata)
  return(df)
}

#simulation function for nsim iterations
sim_func<-function(eff=0.2,t1=2,t2=2,n1=10,n2=10,rho1=0.05,alpha1=0.03,rho2=0.05,alpha2=0.02,sigma1=1,nsim=1000,seed=1112){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #rho1: within-therapist ICC for treatment arm
  #alpha1: between-therapist ICC for treatment arm
  #rho2: within-therapist ICC for control arm
  #alpha2: between-therapist ICC for control arm
  #sigma1: total variance for treatment arm
  
  #nsim: number of simulations 
  #seed: seed number for reproducibility
  
  #model fitting
  set.seed(seed)
  pval<-rep(NA,nsim)         #pval: simulated power
  pval0<-rep(NA,nsim)        #pval0: simulated type I error rate
  I<-ceiling(sample_size(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,rho1=rho1,alpha1=alpha1,rho2=rho2,alpha2=alpha2,sigma1=sigma1))   #I: required number of centers from formula
  pw<-power(I=I,eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,rho1=rho1,alpha1=alpha1,rho2=rho2,alpha2=alpha2,sigma1=sigma1)             #pw: predicted power from formula
   
  count.singular <- 0    #count how many singular fitting cases
  
  for (s in 1:nsim){
    exit <- FALSE        #exit: non-singular fitting indicator 
    while (exit==FALSE){
      df<-data_gen(eff=eff,I=I,t1=t1,t2=t2,n1=n1,n2=n2,rho1=rho1,alpha1=alpha1,rho2=rho2,alpha2=alpha2,sigma1=sigma1)       #df: simulated data set
      fit0 = try(lme(y0 ~ trt, random = list(sid = ~ 1, tid = pdDiag(~trt)),weights = varIdent(form= ~ 1 | trt),            #fit0: model fitting using outcome generated under null
                     data=df),silent=T)
      fit1 = try(lme(y ~ trt, random = list(sid = ~ 1, tid = pdDiag(~trt)),weights = varIdent(form= ~ 1 | trt),             #fit1: model fitting using outcome generated under alternative
                     data=df),silent=T)
      count.singular <- count.singular + 1
      if(class(fit0)!="try-error"&class(fit1)!="try-error"){
        exit <- TRUE
      }
    }
    
    #z-test
    testan <- coef(summary(fit0))[2,1]/coef(summary(fit0))[2,2]       #testan: z-test statistics for outcome generated under null
    pval0[s] <- (min((1-pnorm(testan)),pnorm(testan)))*2              #calculate p-value under null 
    
    testaa <- coef(summary(fit1))[2,1]/coef(summary(fit1))[2,2]       #testaa: z-test statistics for outcome generated under alternative
    pval[s] <- (min((1-pnorm(testaa)),pnorm(testaa)))*2               #calculate p-value under alternative 
  }
  count.singular <- count.singular - nsim    
  tp1<-round(mean(pval0<0.05,na.rm=T) , 3)                            #tp1: empirical type I error rate     
  empower<-round(mean(pval<0.05,na.rm=T) , 3)                         #empower: empirical power
  pw <- round(pw,3)
  
  return(c(rho1,rho2,alpha1,alpha2,I,tp1,pw,empower,count.singular))
}

#simulation set 1 ICCs configuration
rho1 <- c(0.01,0.05,0.10,0.20)                #rho1: within-therapist ICC for treatment arm         
rho2 <- c(0.01,0.05,0.10,0.20)                #rho2: within-therapist ICC for control arm

#Table 2
#generate table 2 column 1
T1 <- 4            #T1: number of therapists in treatment arm
T2 <- 6            #T2: number of therapists in control arm
N1 <- 15           #N1: cluster size in treatment arm
N2 <- 10           #N2: cluster size in control arm
gamma <- 0.2       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs                

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate table 2 column 2
T1 <- 4            #T1: number of therapists in treatment arm
T2 <- 6            #T2: number of therapists in control arm
N1 <- 15           #N1: cluster size in treatment arm
N2 <- 10           #N2: cluster size in control arm
gamma <- 0.5       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs   

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate table 2 column 3
T1 <- 4            #T1: number of therapists in treatment arm
T2 <- 6            #T2: number of therapists in control arm
N1 <- 15           #N1: cluster size in treatment arm
N2 <- 10           #N2: cluster size in control arm
gamma <- 0.8       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#Web Table 1
#generate WEB table 1 column 1
T1 <- 4            #T1: number of therapists in treatment arm
T2 <- 6            #T2: number of therapists in control arm
N1 <- 15           #N1: cluster size in treatment arm
N2 <- 10           #N2: cluster size in control arm
gamma <- 0         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 1 column 2
T1 <- 4            #T1: number of therapists in treatment arm
T2 <- 6            #T2: number of therapists in control arm
N1 <- 15           #N1: cluster size in treatment arm
N2 <- 10           #N2: cluster size in control arm
gamma <- 1         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#Web Table 2
#generate WEB table 2 column 1
T1 <- 2            #T1: number of therapists in treatment arm
T2 <- 4            #T2: number of therapists in control arm
N1 <- 30           #N1: cluster size in treatment arm
N2 <- 15           #N2: cluster size in control arm
gamma <- 0.2       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 2 column 2
T1 <- 2            #T1: number of therapists in treatment arm
T2 <- 4            #T2: number of therapists in control arm
N1 <- 30           #N1: cluster size in treatment arm
N2 <- 15           #N2: cluster size in control arm
gamma <- 0.5       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 2 column 3
T1 <- 2            #T1: number of therapists in treatment arm
T2 <- 4            #T2: number of therapists in control arm
N1 <- 30           #N1: cluster size in treatment arm
N2 <- 15           #N2: cluster size in control arm
gamma <- 0.8       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#Web Table 3
#generate WEB table 3 column 1
T1 <- 2            #T1: number of therapists in treatment arm
T2 <- 4            #T2: number of therapists in control arm
N1 <- 30           #N1: cluster size in treatment arm
N2 <- 15           #N2: cluster size in control arm
gamma <- 0         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 3 column 2
T1 <- 2            #T1: number of therapists in treatment arm
T2 <- 4            #T2: number of therapists in control arm
N1 <- 30           #N1: cluster size in treatment arm
N2 <- 15           #N2: cluster size in control arm
gamma <- 1         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#Web Table 4
#generate WEB table 4 column 1
T1 <- 5            #T1: number of therapists in treatment arm
T2 <- 10           #T2: number of therapists in control arm
N1 <- 12           #N1: cluster size in treatment arm
N2 <- 6            #N2: cluster size in control arm
gamma <- 0.2       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 4 column 2
T1 <- 5            #T1: number of therapists in treatment arm
T2 <- 10           #T2: number of therapists in control arm
N1 <- 12           #N1: cluster size in treatment arm
N2 <- 6            #N2: cluster size in control arm
gamma <- 0.5       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 4 column 3
T1 <- 5            #T1: number of therapists in treatment arm
T2 <- 10           #T2: number of therapists in control arm
N1 <- 12           #N1: cluster size in treatment arm
N2 <- 6            #N2: cluster size in control arm
gamma <- 0.8       #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#Web Table 5
#generate WEB table 5 column 1
T1 <- 5            #T1: number of therapists in treatment arm
T2 <- 10           #T2: number of therapists in control arm
N1 <- 12           #N1: cluster size in treatment arm
N2 <- 6            #N2: cluster size in control arm
gamma <- 0         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

#generate WEB table 5 column 2
T1 <- 5            #T1: number of therapists in treatment arm
T2 <- 10           #T2: number of therapists in control arm
N1 <- 12           #N1: cluster size in treatment arm
N2 <- 6            #N2: cluster size in control arm
gamma <- 1         #gamma: ratio of between-therapist ICCs versus within-therapist ICCs

fidata <- c()      #fidata: storing simulation results
for (i in 1:length(rho1)){
  for (j in 1:length(rho2)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1,t2=T2,n1=N1,n2=N2,
                                    rho1=rho1[i],alpha1=gamma*rho1[i],rho2=rho2[j],alpha2=gamma*rho2[j],
                                    sigma1=1,nsim=5000,seed=920784642))
  }
}

 