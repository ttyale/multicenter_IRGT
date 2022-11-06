#########################################
###   This file reproduces Figure 3   ###
#########################################
require("ggplot2")
require("gridExtra")
require("nlme")
require("latex2exp")

#Sample size calculation based on random effects variance and residual variance
sample_size<-function(eff,t1,t2,n1,n2,sigma_lambda,sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #sigma_lambda: center-level random effect variance
  #sigma_theta1: therapist random effect variance in treatment arm
  #sigma_theta2: therapist random effect variance in control arm
  #sigma_epsilon1: individual residual variance in treatment arm
  #sigma_epsilon2: individual residual variance in control arm
  
  alpha=0.05 # type I error rate
  beta=0.2   # type II error rate
  
  sigma1 <- sigma_lambda+sigma_theta1+sigma_epsilon1 #sigma1: total variance for treatment arm
  sigma2 <- sigma_lambda+sigma_theta2+sigma_epsilon2 #sigma2: total variance for control arm
  rho1 <- sigma_lambda/sigma1                        #rho1: therapist ICC in treatment arm
  rho2 <- (sigma_lambda+sigma_theta1)/sigma1         #rho2: therapist ICC in control arm
  alpha1 <- sigma_lambda/sigma2                      #alpha1: center ICC in treatmeent arm
  alpha2 <- (sigma_lambda+sigma_theta2)/sigma2       #alpha2: center ICC in control arm
  
  vb<-sigma1/n1/t1*((1+(n1-1)*rho2)+n1*(t1-1)*rho1)+sigma2/n2/t2*((1+(n2-1)*alpha2)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)
  I <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))
  return(I)
}

#power calculation based on random effects variance and residual variance
power<-function(I,eff,t1,t2,n1,n2,sigma_lambda,sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2){
  #I: number of centers
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #sigma_lambda: center-level random effect variance
  #sigma_theta1: therapist random effect variance in treatment arm
  #sigma_theta2: therapist random effect variance in control arm
  #sigma_epsilon1: individual residual variance in treatment arm
  #sigma_epsilon2: individual residual variance in control arm
  
  alpha=0.05  # type I error rate
  
  sigma1 <- sigma_lambda+sigma_theta1+sigma_epsilon1 #sigma1: total variance for treatment arm
  sigma2 <- sigma_lambda+sigma_theta2+sigma_epsilon2 #sigma2: total variance for control arm
  rho1 <- sigma_lambda/sigma1                        #rho1: therapist ICC in treatment arm
  rho2 <- (sigma_lambda+sigma_theta1)/sigma1         #rho2: therapist ICC in control arm
  alpha1 <- sigma_lambda/sigma2                      #alpha1: center ICC in treatmeent arm
  alpha2 <- (sigma_lambda+sigma_theta2)/sigma2       #alpha2: center ICC in control arm
  
  vb<-sigma1/n1/t1/I*((1+(n1-1)*rho2)+n1*(t1-1)*rho1)+sigma2/n2/t2/I*((1+(n2-1)*alpha2)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)/I
  pw<-pnorm(eff/sqrt(vb)-qnorm(1-alpha/2))
  
  return(pw)
}

#data generation based on random effects variance and residual variance
data_gen<-function(I=5,eff=0.2,t1=2,t2=3,n1=10,n2=13,sigma_lambda,sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2){
  #I: number of centers
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #sigma_lambda: center-level random effect variance
  #sigma_theta1: therapist random effect variance in treatment arm
  #sigma_theta2: therapist random effect variance in control arm
  #sigma_epsilon1: individual residual variance in treatment arm
  #sigma_epsilon2: individual residual variance in control arm
  
  sigma1 <- sigma_lambda+sigma_theta1+sigma_epsilon1 #sigma1: total variance for treatment arm
  sigma2 <- sigma_lambda+sigma_theta2+sigma_epsilon2 #sigma2: total variance for control arm
  rho1 <- sigma_lambda/sigma1                        #rho1: therapist ICC in treatment arm
  rho2 <- (sigma_lambda+sigma_theta1)/sigma1         #rho2: therapist ICC in control arm
  alpha1 <- sigma_lambda/sigma2                      #alpha1: center ICC in treatment arm
  alpha2 <- (sigma_lambda+sigma_theta2)/sigma2       #alpha2: center ICC in control arm
  
  intercept=0.5                                      #intercept: beta0
  N=I*(t1*n1+t2*n2)                                  #N: total number of patients for two arms
  patid <- 1:N                                       #patid: patient id 
  trt<-rep(c(rep(1,t1*n1),rep(0,t2*n2)),I)           #trt: treatment indicator X_ijk
  
  strata<-rep(1:I,each=(t1*n1+t2*n2))                      #strata: center id
  tid_s<-c(rep(1:t1,each=n1),rep((t1+1):(t1+t2),each=n2))  #tid_s: therapist id for each center
  tid<-c()                                                 #tid: overall therapist id 
  for(i in 1:I){
    tid<-c(tid,((i-1)*(t1+t2)+tid_s))
  }
  
  I_effect<-rep(rnorm(I,0,sqrt(sigma_lambda)),each=(t1*n1+t2*n2))   #I_effect: center random effect
  
  therapist_effect1<-rnorm(I*(t1+t2),0,sqrt(sigma_theta1))          #therapist_effect1: therapist random effect for treatment arm
  therapist_effect2<-rnorm(I*(t1+t2),0,sqrt(sigma_theta2))          #therapist_effect2: therapist random effect for control arm
  
  therapist_effect<-rep(therapist_effect1,as.numeric(table(tid)))*trt+rep(therapist_effect2,as.numeric(table(tid)))*(1-trt) #therapist_effect: combined therapist random effect
  
  epsilon<-c()       #epsilon: combined measurement error for two arms                                           
  for (i in 1:I){
    epsilon_s <-c(rnorm(t1*n1,0,sqrt(sigma_epsilon1)),rnorm(t2*n2,0,sqrt(sigma_epsilon2)))
    epsilon<-c(epsilon,epsilon_s)
  }
  
  y=intercept-eff*trt+I_effect+therapist_effect+epsilon #y: outcome under alternative
  y0=intercept+I_effect+therapist_effect+epsilon        #y0: outcome under null
  
  df<-data.frame(y=y,y0=y0,trt=trt,tid=tid,sid=strata)
  return(df)
}

#simulation function for nsim iterations
sim_func<-function(eff=0.2,t1=2,t2=2,n1=10,n2=10,sigma_lambda,sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2,nsim=1000){
  #eff: effect size 
  #t1: number of therapists in treatment arm
  #n1: cluster size in treatment arm
  #t2: number of therapists in control arm
  #n2: cluster size in control arm
  
  #sigma_lambda: center-level random effect variance
  #sigma_theta1: therapist random effect variance in treatment arm
  #sigma_theta2: therapist random effect variance in control arm
  #sigma_epsilon1: individual residual variance in treatment arm
  #sigma_epsilon2: individual residual variance in control arm
  
  #nsim: number of simulations 
  
  #model fitting
  pval<-rep(NA,nsim)      #pval: simulated power
  pval0<-rep(NA,nsim)     #pval0: simulated type I error rate
  I<-ceiling(sample_size(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda,sigma_theta1=sigma_theta1,         #I: required number of centers from formula
                            sigma_theta2=sigma_theta2,sigma_epsilon1=sigma_epsilon1,sigma_epsilon2=sigma_epsilon2))  
  pw<-power(I=I,eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda,sigma_theta1=sigma_theta1,                  #pw: predicted power from formula
            sigma_theta2=sigma_theta2,sigma_epsilon1=sigma_epsilon1,sigma_epsilon2=sigma_epsilon2)
  
  count.singular <- 0  #count how many singular fitting cases
  
  for (s in 1:nsim){
    exit <- FALSE     #exit: non-singular fitting indicator                      
    while (exit==FALSE){
      df<-data_gen(eff=eff,I=I,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda,sigma_theta1=sigma_theta1,        #df: simulated data set
                   sigma_theta2=sigma_theta2,sigma_epsilon1=sigma_epsilon1,sigma_epsilon2=sigma_epsilon2)
      fit0 = try(lme(y0 ~ trt, random = list(sid = ~ 1, tid = pdDiag(~trt)),weights = varIdent(form= ~ 1 | trt),   #fit0: model fitting using outcome generated under null
                     control=list(msMaxIter = 1000, msMaxEval = 1000), data=df),silent=T)
      fit1 = try(lme(y ~ trt, random = list(sid = ~ 1, tid = pdDiag(~trt)),weights = varIdent(form= ~ 1 | trt),    #fit1: model fitting using outcome generated under alternative
                     control=list(msMaxIter = 1000, msMaxEval = 1000), data=df),silent=T)
      count.singular <- count.singular + 1                                                  
      if(class(fit0)!="try-error"&class(fit1)!="try-error"){
        exit <- TRUE
      }
    }
    
    #z-test
    testan <- coef(summary(fit0))[2,1]/coef(summary(fit0))[2,2]      #testan: z-test statistics for outcome generated under null
    pval0[s] <- (min((1-pnorm(testan)),pnorm(testan)))*2             #calculate p-value under null 
    
    testaa <- coef(summary(fit1))[2,1]/coef(summary(fit1))[2,2]      #testaa: z-test statistics for outcome generated under alternative
    pval[s] <- (min((1-pnorm(testaa)),pnorm(testaa)))*2              #calculate p-value under alternative 
    
  }
  
  count.singular <- count.singular - nsim                     
  tp1<-round(mean(pval0<0.05,na.rm=T) , 3)                      #tp1: empirical type I error rate     
  empower<-round(mean(pval<0.05,na.rm=T) , 3)                   #empower: empirical power
  pw <- round(pw,3)                   
  
  return(c(I,tp1,pw,empower,count.singular))
}


#set simulation parameters for panel (A) and (B)
#(A),(B) correspond to the first and second column of parameters sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2
sigma_lambda <- seq(0.02,0.16,by=0.02)
sigma_theta1 <- c(0.01,0.05) 
sigma_theta2 <- c(0.02,0.10)
sigma_epsilon1 <- c(0.97,0.9)
sigma_epsilon2 <- c(0.96,0.9)
eff <- 0.2
t1 <- 4
t2 <-6
n1 <- 15
n2 <- 10

#set seed
set.seed(920784642)

#data set storing the simulation results for panel (A)
fidata1 <- c()
#panel (A)
i <- 1         #(A) corresponds to the first column    
#iterates across sigma_lambda
for (j in 1:length(sigma_lambda)){         
  fidata1 <- rbind(fidata1,sim_func(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda[j],sigma_theta1=sigma_theta1[i],
                                  sigma_theta2=sigma_theta2[i],sigma_epsilon1=sigma_epsilon1[i],sigma_epsilon2=sigma_epsilon2[i],nsim=5000))
}
colnames(fidata1) <- c("I","empirical.type.I","predicted.power","empirical.power","n.singular") #name the simulation results

#panel (B)
#data set storing the simulation results for panel (B)
fidata2 <- c()
i <- 2          #(B) corresponds to the second column  
#iterates across sigma_lambda
for (j in 1:length(sigma_lambda)){
  fidata2 <- rbind(fidata2,sim_func(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda[j],sigma_theta1=sigma_theta1[i],
                                   sigma_theta2=sigma_theta2[i],sigma_epsilon1=sigma_epsilon1[i],sigma_epsilon2=sigma_epsilon2[i],nsim=5000))
}
colnames(fidata2) <- c("I","empirical.type.I","predicted.power","empirical.power","n.singular") #name the simulation results

#set simulation parameters for panel (C) and (D)
#(C),(D) correspond to the first and second column of parameters sigma_theta1,sigma_theta2,sigma_epsilon1,sigma_epsilon2
sigma_lambda <- seq(0.02,0.16,by=0.02)
sigma_theta1 <- c(0.02,0.04)
sigma_theta2 <- c(0,0)
sigma_epsilon1 <- c(0.97,0.9)
sigma_epsilon2 <- c(0.96,0.4)

#panel (C)
#data set storing the simulation results for panel (C)
fidata3 <- c() 
i <- 1       #(C) corresponds to the first column 
#iterates across sigma_lambda
for (j in 1:length(sigma_lambda)){
  fidata3 <- rbind(fidata3,sim_func(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda[j],sigma_theta1=sigma_theta1[i],
                                  sigma_theta2=sigma_theta2[i],sigma_epsilon1=sigma_epsilon1[i],sigma_epsilon2=sigma_epsilon2[i],nsim=5000))
}
colnames(fidata3) <- c("I","empirical.type.I","predicted.power","empirical.power","n.singular")  #name the simulation results

#panel (D)
#data set storing the simulation results for panel (D)
fidata4 <- c()
i <- 2       #(D) corresponds to the second column 
#iterates across sigma_lambda
for (j in 1:length(sigma_lambda)){
  fidata4 <- rbind(fidata4,sim_func(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,sigma_lambda=sigma_lambda[j],sigma_theta1=sigma_theta1[i],
                                  sigma_theta2=sigma_theta2[i],sigma_epsilon1=sigma_epsilon1[i],sigma_epsilon2=sigma_epsilon2[i],nsim=5000))
}
colnames(fidata4) <- c("I","empirical.type.I","predicted.power","empirical.power","n.singular") #name the simulation results

#generate plot
as.data.frame(fidata1)
as.data.frame(fidata2)
as.data.frame(fidata3)
as.data.frame(fidata4)
#add column for sigma_lambda
fidata1$sigmalambda <- sigma_lambda
fidata2$sigmalambda <- sigma_lambda
fidata3$sigmalambda <- sigma_lambda
fidata4$sigmalambda <- sigma_lambda

#function generating figure 3 
plotfig3 <- function(z,xlab,ylab,title,pw){
  #z: data set containing all simulation results for each panel
  #xlab: labels for x-axis of each panel
  #ylab: labels for y-axis of each panel
  #title: title for each panel
  #pw: predicted power cauclated from formula
  
  p1 <- ggplot(z, aes(x=sigmalambda, y=empirical.power)) + 
    geom_point(size=2.5) + 
    geom_line() +
    geom_hline(yintercept=pw+0.012, lty=2, alpha=0.5) +  #dashed line for upper bound of 95% CI around predicted power
    geom_hline(yintercept=pw-0.012, lty=2, alpha=0.5) +  #dashed line for lower bound of 95% CI around predicted power
    xlab(TeX(xlab)) +
    ylab(TeX(ylab)) +
    ggtitle(TeX(title)) +
    xlim(0,0.18) +
    ylim(0.75,0.9) +
    theme_bw(base_size = 12)
  return(p1)
}

#generate panel (A)
p1 <- plotfig3(fidata1,"$\\sigma_\\lambda^2$","$Empirical Power$",
           "$(A): I=8,\\sigma_{\\theta_1}^2=0.01,\\sigma_{\\theta_2}^2=0.02,\\sigma_{\\epsilon_1}^2=0.97,\\sigma_{\\epsilon_2}^2=0.96$",fidata1$predicted.power)
#generate panel (B)
p2 <- plotfig3(fidata2,"$\\sigma_\\lambda^2$","$Empirical Power$",
           "$(B): I=12,\\sigma_{\\theta_1}^2=0.05,\\sigma_{\\theta_2}^2=0.10,\\sigma_{\\epsilon_1}^2=0.9,\\sigma_{\\epsilon_2}^2=0.9$",fidata2$predicted.power)
#generate panel (C)
p3 <- plotfig3(fidata3,"$\\sigma_\\lambda^2$","$Empirical Power$",
           "$(C): I=8,\\sigma_{\\theta_1}^2=0.02,\\sigma_{\\theta_2}^2=0,\\sigma_{\\epsilon_1}^2=0.97,\\sigma_{\\epsilon_2}^2=0.96$",fidata3$predicted.power)
#generate panel (D)
p4 <- plotfig3(fidata4,"$\\sigma_\\lambda^2$","$Empirical Power$",
           "$(D): I=7,\\sigma_{\\theta_1}^2=0.04,\\sigma_{\\theta_2}^2=0,\\sigma_{\\epsilon_1}^2=0.9,\\sigma_{\\epsilon_2}^2=0.4$",fidata4$predicted.power)
#combine 4 panels to produce figure 3
grid.arrange(p1,p2,p3,p4)

