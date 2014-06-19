
model {

#Intercept Priors on N and p
mu.lam ~ dnorm(0,0.001)  
mu.p	~ dunif(0,0.001)
		
sigma.lam ~ dunif(0,10)  
tau.lam <- pow(sigma.lam,-2)    
sigma.p ~ dunif(0,10)
tau.p <- pow(sigma.p,-2)
    
#Detection Covariate Priors    
b.temp~dunif(0,1)
mu.temp <- log(b.temp/(1-b.temp))
s.temp~dunif(0,0.5)
tau.temp <- pow(s.temp,-2)
    
b.jd~dunif(0,1)
mu.jd <- log(b.jd/(1-b.jd))
s.jd~dunif(0,10)
tau.jd <- pow(s.jd,-2)
    
b.eff~dunif(0,1)
mu.eff <- log(b.eff/(1-b.eff))
s.eff~dunif(0,7)
tau.eff <- pow(s.eff,-2)
    
#Treatment effect on detection priors
    
pb.ccl~dunif(0,1)
pmu.ccl <- log(pb.ccl/(1-pb.ccl))
ps.ccl~dunif(0,7)
ptau.ccl <- pow(ps.ccl,-2)

pb.csh~dunif(0,1)
pmu.csh <- log(pb.csh/(1-pb.ccl))
ps.csh~dunif(0,7)
ptau.csh <- pow(ps.csh,-2)

pb.04~dunif(0,1)
pmu.04 <- log(pb.04/(1-pb.04))
ps.04~dunif(0,7)
ptau.04 <- pow(ps.04,-2)

pb.2~dunif(0,1)
pmu.2 <- log(pb.2/(1-pb.2))
ps.2~dunif(0,7)
ptau.2 <- pow(ps.2,-2)
    
#Abundance Priors
b.aspect~dunif(0,1)
mu.aspect <- log(b.aspect/(1-b.aspect))
s.aspect~dunif(0,7)
tau.aspect <- pow(s.aspect,-2)
    
b.ccl~dunif(0,1)
mu.ccl <- log(b.ccl/(1-b.ccl))
s.ccl~dunif(0,7)
tau.ccl <- pow(s.ccl,-2)

b.csh~dunif(0,1)
mu.csh <- log(b.csh/(1-b.ccl))
s.csh~dunif(0,7)
tau.csh <- pow(s.csh,-2)

b.04~dunif(0,1)
mu.04 <- log(b.04/(1-b.04))
s.04~dunif(0,7)
tau.04 <- pow(s.04,-2)

b.2~dunif(0,1)
mu.2 <- log(b.2/(1-b.2))
s.2~dunif(0,7)
tau.2 <- pow(s.2,-2)    
    
b.mast~dunif(0,1)
mu.mast <- log(b.mast/(1-b.mast))
s.mast~dunif(0,7)
tau.mast <- pow(s.mast,-2)

   
#Model Specification
for (i in 1:nspecies){
  #Generate random slopes and intercepts for each species
  alpha.lam[i] ~ dnorm(mu.lam,tau.lam)
  alpha.p[i] ~ dnorm(mu.p, tau.p)
  beta.aspect[i] ~ dnorm(mu.aspect,tau.aspect)
  beta.c.cl[i] ~ dnorm(mu.ccl,tau.ccl)
  beta.c.sh[i] ~ dnorm(mu.csh,tau.csh)
  beta.temp[i] ~ dnorm(mu.temp,tau.temp)
  beta.jd[i] ~ dnorm(mu.jd,tau.jd)
  beta.eff[i] ~ dnorm(mu.eff,tau.eff)
  beta.mast[i] ~ dnorm(mu.mast,tau.mast)
  beta.04[i] ~ dnorm(mu.04,tau.04)
  beta.2[i] ~ dnorm(mu.2,tau.2)
  pbeta.c.cl[i] ~ dnorm(pmu.ccl,ptau.ccl)
  pbeta.c.sh[i] ~ dnorm(pmu.csh,ptau.csh)
  pbeta.04[i] ~ dnorm(pmu.04,ptau.04)
  pbeta.2[i] ~ dnorm(pmu.2,ptau.2)
  #Iterate through sites in first 3 years (all sites)
  for (j in 1:nsites){
    for (t in 1:3){
      N[j,i,t] ~ dpois(lambda[j,i,t])
      log(lambda[j,i,t]) <- alpha.lam[i] 
        + beta.aspect[i]*aspect[j] + beta.c.cl[i]*treat[j,1,t] 
        + beta.c.sh[i]*treat[j,2,t] + beta.mast[i]*mast[j,t]
        + beta.04[i]*treat[j,4,t] + beta.2[i]*treat[j,5,t]
      for (k in sindex[i]:5){
        lin.pred[j,k,i,t] <- alpha.p[i] + beta.temp[i]*temp[j,t] + beta.jd[i]*jd[j,t]
        + beta.eff[i]*eff[j,i,t] + pbeta.c.cl[i]*treat[j,1,t] 
        + pbeta.c.sh[i]*treat[j,2,t] + pbeta.04[i]*treat[j,4,t] + pbeta.2[i]*treat[j,5,t]
        p[j,k,i,t] <- 1/(1+exp(-lin.pred[j,k,i,t])) #logit() DOES NOT WORK
        counts[j,k,i,t] ~ dbin(p[j,k,i,t],N[j,i,t])
        
        
      }
    pmean[j,i,t] <- sum(p[j,2:5,i,t]) / 4

    }
    #Iterate through sites in years 4 and 5 (not all sampled)
    for (t in first[j]:5){
      N[j,i,t] ~ dpois(lambda[j,i,t])
      log(lambda[j,i,t]) <- alpha.lam[i]
        + beta.aspect[i]*aspect[j] + beta.c.cl[i]*treat[j,1,t] 
        + beta.c.sh[i]*treat[j,2,t] + beta.mast[i]*mast[j,t]
        + beta.04[i]*treat[j,4,t] + beta.2[i]*treat[j,5,t]
      for (k in sindex[i]:5){
        lin.pred[j,k,i,t] <- alpha.p[i] + beta.temp[i]*temp[j,t] + beta.jd[i]*jd[j,t]
        + beta.eff[i]*eff[j,i,t] + pbeta.c.cl[i]*treat[j,1,t] 
        + pbeta.c.sh[i]*treat[j,2,t] + pbeta.04[i]*treat[j,4,t] + pbeta.2[i]*treat[j,5,t]
        p[j,k,i,t] <- 1/(1+exp(-lin.pred[j,k,i,t]))
        counts[j,k,i,t] ~ dbin(p[j,k,i,t],N[j,i,t])
      }
      pmean[j,i,t] <- sum(p[j,2:5,i,t]) / 4
}
}}
    
#Derived parameters
    
for (i in 1:2){ #just mice and chipmunks??
    
    pre.control.mean[i] <-   (sum(N[1:2,i,1])+sum(N[7:8,i,1])+sum(N[13:16,i,1])
                                +sum(N[1:2,i,2])+sum(N[7:8,i,2])+sum(N[13:16,i,2])) / 16
    pre.clear.mean[i] <-  (sum(N[9,i,1:2])+sum(N[11,i,1:2])+sum(N[17,i,1:2])+sum(N[20,i,1:2])
                                +sum(N[31,i,1:2])+sum(N[32,i,1:2])) / 12
    pre.shelter.mean[i] <- (sum(N[10,i,1:2])+sum(N[12,i,1:2])+sum(N[18,i,1:2])+sum(N[19,i,1:2])
                                  +sum(N[29,i,1:2])+sum(N[30,i,1:2])) / 12
    pre.2.mean[i] <-   (sum(N[3:4,i,1])+sum(N[22,i,1])+sum(N[26,i,1])+sum(N[28,i,1])+
                        sum(N[3:4,i,2])+sum(N[22,i,2])+sum(N[26,i,2])+sum(N[28,i,2])) / 10   
    pre.04.mean[i] <-  (sum(N[5:6,i,1])+sum(N[21,i,1])+sum(N[23:25,i,1])+sum(N[27,i,1])+
                        sum(N[5:6,i,2])+sum(N[21,i,2])+sum(N[23:25,i,2])+sum(N[27,i,2])) / 14
 
   
    ##Post-harvest
   
    post.control.mean[i] <-   (sum(N[1:2,i,3])+sum(N[7:8,i,3])+sum(N[13:16,i,3])
                                +sum(N[7:8,i,4])+sum(N[13:16,i,4])
                                +sum(N[1:2,i,5])+sum(N[7:8,i,5])+sum(N[13:16,i,5])) / 22
    post.clear.mean[i] <-  (sum(N[9,i,3:5])+sum(N[11,i,3:5])+sum(N[17,i,3:5])+sum(N[20,i,3:5])
                                +sum(N[31,i,3:5])+sum(N[32,i,3:5])) / 18
    post.shelter.mean[i] <- (sum(N[10,i,3:5])+sum(N[12,i,3:5])+sum(N[18,i,3:5])+sum(N[19,i,3:5])
                                  +sum(N[29,i,3:5])+sum(N[30,i,3:5])) / 18  
    post.2.mean[i] <-   (sum(N[3:4,i,3])+sum(N[22,i,3])+sum(N[26,i,3])+sum(N[28,i,3])+
                        sum(N[3:4,i,5])+sum(N[22,i,5])+sum(N[26,i,5])+sum(N[28,i,5])) / 10   
    post.04.mean[i] <-  (sum(N[5:6,i,3])+sum(N[21,i,3])+sum(N[23:25,i,3])+sum(N[27,i,3])
                        +sum(N[24:25,i,4])+
                        sum(N[5:6,i,5])+sum(N[21,i,5])+sum(N[23:25,i,5])+sum(N[27,i,5])) / 16
}
############Mean p values by treatment

for (i in 1:4){ #just mice and chipmunks??

    pre.control.p[i] <-   (sum(pmean[1:2,i,1])+sum(pmean[7:8,i,1])+sum(pmean[13:16,i,1])
                                +sum(pmean[1:2,i,2])+sum(pmean[7:8,i,2])+sum(pmean[13:16,i,2])) / 16
    pre.clear.p[i] <-  (sum(pmean[9,i,1:2])+sum(pmean[11,i,1:2])+sum(pmean[17,i,1:2])+sum(pmean[20,i,1:2])
                                +sum(pmean[31,i,1:2])+sum(pmean[32,i,1:2])) / 12
    pre.shelter.p[i] <- (sum(pmean[10,i,1:2])+sum(pmean[12,i,1:2])+sum(pmean[18,i,1:2])+sum(pmean[19,i,1:2])
                                  +sum(pmean[29,i,1:2])+sum(pmean[30,i,1:2])) / 12
    pre.2.p[i] <-   (sum(pmean[3:4,i,1])+sum(pmean[22,i,1])+sum(pmean[26,i,1])+sum(pmean[28,i,1])+
                        sum(pmean[3:4,i,2])+sum(pmean[22,i,2])+sum(pmean[26,i,2])+sum(pmean[28,i,2])) / 10   
    pre.04.p[i] <-  (sum(pmean[5:6,i,1])+sum(pmean[21,i,1])+sum(pmean[23:25,i,1])+sum(pmean[27,i,1])+
                        sum(pmean[5:6,i,2])+sum(pmean[21,i,2])+sum(pmean[23:25,i,2])+sum(pmean[27,i,2])) / 14
 
   
    ##Post-harvest
   
    post.control.p[i] <-   (sum(pmean[1:2,i,3])+sum(pmean[7:8,i,3])+sum(pmean[13:16,i,3])
                                +sum(pmean[7:8,i,4])+sum(pmean[13:16,i,4])
                                +sum(pmean[1:2,i,5])+sum(pmean[7:8,i,5])+sum(pmean[13:16,i,5])) / 22
    post.clear.p[i] <-  (sum(pmean[9,i,3:5])+sum(pmean[11,i,3:5])+sum(pmean[17,i,3:5])+sum(pmean[20,i,3:5])
                                +sum(pmean[31,i,3:5])+sum(pmean[32,i,3:5])) / 18
    post.shelter.p[i] <- (sum(pmean[10,i,3:5])+sum(pmean[12,i,3:5])+sum(pmean[18,i,3:5])+sum(pmean[19,i,3:5])
                                  +sum(pmean[29,i,3:5])+sum(pmean[30,i,3:5])) / 18  
    post.2.p[i] <-   (sum(pmean[3:4,i,3])+sum(pmean[22,i,3])+sum(pmean[26,i,3])+sum(pmean[28,i,3])+
                        sum(pmean[3:4,i,5])+sum(pmean[22,i,5])+sum(pmean[26,i,5])+sum(pmean[28,i,5])) / 10   
    post.04.p[i] <-  (sum(pmean[5:6,i,3])+sum(pmean[21,i,3])+sum(pmean[23:25,i,3])+sum(pmean[27,i,3])
                        +sum(pmean[24:25,i,4])+
                        sum(pmean[5:6,i,5])+sum(pmean[21,i,5])+sum(pmean[23:25,i,5])+sum(pmean[27,i,5])) / 16
    
    p.all[i] <- (22*post.control.p[i]+18*post.clear.p[i]+18*post.shelter.p[i]+10*post.2.p[i]+16*post.04.p[i]+
                   16*pre.control.p[i]+12*pre.clear.p[i]+12*pre.shelter.p[i]+10*pre.2.p[i]+14*pre.04.p[i]) / 148
}   
#####################################################
       
for (i in 1:nspecies){
    for (t in 1:2){
      yearly.control[i,t] <- (sum(N[1:2,i,t])+sum(N[7:8,i,t])+sum(N[13:16,i,t])) / 8
                        
      yearly.clear[i,t] <-  (sum(N[9,i,t])+sum(N[11,i,t])+sum(N[17,i,t])+sum(N[20,i,t])
                                +sum(N[31,i,t])+sum(N[32,i,t])) / 6
      yearly.shelter[i,t] <- (sum(N[10,i,t])+sum(N[12,i,t])+sum(N[18,i,t])+sum(N[19,i,t])
                                  +sum(N[29,i,t])+sum(N[30,i,t])) / 6
      yearly.patch04[i,t] <-  (sum(N[5:6,i,t])+N[21,i,t]+sum(N[23:25,i,t]) + N[27,i,t]) / 7
    
      yearly.patch2[i,t] <- (sum(N[3:4,i,t])+N[22,i,t]+ N[26,i,t] + N[28,i,t]) / 5
    }
    for (t in 3:3){
      yearly.control[i,t] <- (sum(N[1:2,i,t])+sum(N[7:8,i,t])+sum(N[13:16,i,t])) / 8
                        
      yearly.clear[i,t] <-  (sum(N[9,i,t])+sum(N[11,i,t])+sum(N[17,i,t])+sum(N[20,i,t])
                                +sum(N[31,i,t])+sum(N[32,i,t])) / 6
      yearly.shelter[i,t] <- (sum(N[10,i,t])+sum(N[12,i,t])+sum(N[18,i,t])+sum(N[19,i,t])
                                  +sum(N[29,i,t])+sum(N[30,i,t])) / 6
      yearly.patch04[i,t] <-  (sum(N[5:6,i,t])+N[21,i,t]+sum(N[23:25,i,t]) + N[27,i,t]) / (7)
    
      yearly.patch2[i,t] <- (sum(N[3:4,i,t])+N[22,i,t]+ N[26,i,t] + N[28,i,t]) / 5
    }
    for (t in 4:4){
      yearly.control[i,t] <- (sum(N[7:8,i,t])+sum(N[13:16,i,t])) / 6
                        
      yearly.clear[i,t] <-  (sum(N[9,i,t])+sum(N[11,i,t])+sum(N[17,i,t])+sum(N[20,i,t])
                                +sum(N[31,i,t])+sum(N[32,i,t])) / 6
      yearly.shelter[i,t] <- (sum(N[10,i,t])+sum(N[12,i,t])+sum(N[18,i,t])+sum(N[19,i,t])
                                  +sum(N[29,i,t])+sum(N[30,i,t])) / 6
      yearly.patch2[i,t] <- 0
    
      yearly.patch04[i,t] <-  (sum(N[24:25,i,t])) / (2)
    
    }  
    for (t in 5:5){
      yearly.control[i,t] <- (sum(N[1:2,i,t])+sum(N[7:8,i,t])+sum(N[13:16,i,t])) / 8
                        
      yearly.clear[i,t] <-  (sum(N[9,i,t])+sum(N[11,i,t])+sum(N[17,i,t])+sum(N[20,i,t])
                                +sum(N[31,i,t])+sum(N[32,i,t])) / 6
      yearly.shelter[i,t] <- (sum(N[10,i,t])+sum(N[12,i,t])+sum(N[18,i,t])+sum(N[19,i,t])
                                  +sum(N[29,i,t])+sum(N[30,i,t])) / 6
      yearly.patch04[i,t] <-  (sum(N[5:6,i,t])+N[21,i,t]+sum(N[23:25,i,t]) + N[27,i,t]) / (7)
    
      yearly.patch2[i,t] <- (sum(N[3:4,i,t])+N[22,i,t]+ N[26,i,t] + N[28,i,t]) / 5
    }
}    
}

