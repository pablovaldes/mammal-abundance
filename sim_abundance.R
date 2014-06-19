################################################################
## Simulation code for comparison of binomial/poisson         ##
## mixture model of abundance with traditional capture models ##
################################################################

#Utility functions

# Simulate a site's individual encounter histories with this function;
# Every individual gets a history with 5 occasions
# The output from this can be put directly into Rcapture once NAs are removed
gen.data <- function(Nvec,p,treat,bef,sef,hvar,tef){
  
  nsites <- length(Nvec)
  
  #Create empty output array
  eh <- array(data=NA, dim=c(max(Nvec),5,nsites))
  
  #Iterate over sites
  for (j in 1:nsites){
    #set site/treatment effect
    site.effect <- 0
    if(treat[j]==1){site.effect=sef} 
    
    #Iterate over individual animals at site j
    for (i in 1:Nvec[j]){
      #set heterogeneity effect
      het.effect <- rnorm(1,mean=0,sd=hvar)  
      
      #Iterate over encounters for animal i and site j
      for (k in 1:5){
        #set time effect (increases or decreases linearly with time)
        time.effect <- tef*k
        #set behavioral response (changes if animal has been trapped at least once)
        behave.effect <- 0
        if(k>1){
          if(sum(eh[i,1:(k-1),j])>0){
            behave.effect <- bef
          }
        }
        
        intercept <- exp(p)/(1+exp(p))
        
        #Calculate linear predictor
        lin.pred <- intercept + het.effect + site.effect + behave.effect + time.effect
        #Transform to 0-1 scale
        p.eff <- log(lin.pred/(1-lin.pred))
        
        #Generate observed encounter history for animal i,j
        eh[i,k,j] <- rbinom(1,1,p.eff)
      }      
    }
  }
  return(eh)    
}

#Obtain count data for Binomial-Poisson mixture model with this function
#Generates counts of observed individuals at a site at each sampling occasion
gen.mix <- function(data){
  nsites <- dim(data)[3]
  counts <- array(data=NA, dim=c(nsites,5))
  for (i in 1:nsites){
    for (j in 1:5){
      counts[i,j] <- sum(na.omit(data[,j,i]))
    }}
  return(counts)
}

################################################################

#Simulation function Documentation:

# nsims: number of iterations to simulate
# Nvec: vector of abundance values
# p: base (mean) value for deteciton (p)
# bef: behavioral effect for animals captured at least once
# sef: site effect for sites with treatment=1 
# treat: vector of 0s and 1s designating site treatment
# hvar: standard deviation of individual heterogeneity effect
# bef/sef/hvar have linear relationship with p on logit scale
# tef: coefficient on effect of time (i.e., trapping occasion) on detection
# model.file: location of BUGS model to be used

################################################################

#Function code

fitsim <- function(nsims, Nvec, p, bef=0,sef=0,treat=rep(0,length(Nvec)),
                   hvar=0,tef=0,model.file='models/model_sim_Nmix.R'){

nsites <- length(Nvec)

#Create empty output container
output <- array(data=NA, dim=c(nsites,7,6,nsims))

#Load necessary libraries
require(jagsUI)
require(Rcapture)

#Run actual simulation

for (i in 1:nsims){
  #Generate encounter histories
  raw <- gen.data(Nvec,p,treat,bef,sef,hvar,tef)
  
  #Convert raw data to counts
  counts <- gen.mix(raw)
  
  #Initial values for N
  Nst <- array(data=NA,dim=c(nsites))
  for (k in 1:nsites){
    Nst[k] <- max(counts[k,1:5])+1  	
  }
  
  inits <- function(){  
    list(N=Nst      
    )}
  
  print(inits())
  
  #Bundle data
  data <- list("counts","nsites","treat")

  	
  # Parameters to save
  params <- c('mu.lam','p',	"N")
     
  #Send everything to JAGS with specified parameters  

  fit <- jags(data,inits,params,model.file, n.chains=3, 
             n.iter=30000, n.burnin=15000, n.thin=30, parallel=TRUE)

  #Begin generation of output file
  #Structure is as follows - nsites x generated values x model type x iteration
  #Values for site j: 
  #1. Actual N, 2. estimated N, 3. est. N se, 4. est. p, 5. est. p se,
  #6. %bias     7. Is actual N in 95% confidence interval for est? (indicator)
  
  #Generate output for traditional capture models at each site j
  #Currently using M0 (true, since p is constant), Mhc = Mh Chao, MthC = Mth Chao, MKNA, ncaps
  #More or different models can easily be added
  for (j in 1:nsites){ 
    #Clean up data (i.e., remove individuals that were never detected)
    site.data <- na.omit(raw[,,j])
    
    #Generate MKNA for site
    alive.test = vector(length=(dim(site.data)[1]))
    for (k in 1:(dim(site.data)[1])){
      if(sum(site.data[k,])>0){
        alive.test[k]=1}
    }
    mkna = sum(alive.test)
    
    #Generate number of captures for site
    ncaps = sum(site.data)
    
    #Plug into output file
    try(output[j,1,1:6,i] <- Nvec[j],silent=TRUE)
    try(output[j,2,5,i] <- mkna,silent=TRUE)
    try(output[j,2,6,i] <- ncaps,silent=TRUE)   
    
    #Run MRR analysis if animals were captured
    if(sum(site.data)>0){
      #Run in Rcapture
      cap.fit <- closedp.t(site.data)
      position <- c(NA,1,3,7)
      #Name of models to save results for
      mnames <- c(NA,"M0","MhC", "MthC")
        
    #Move results to output file for each MRR model
    for (m in 2:length(position)){
      #Estimated N
      try(nest <- cap.fit$results[position[m],1],silent=TRUE)
      #Standard error
      try(se <- cap.fit$results[position[m],2],silent=TRUE)
      #2.5 and 97.5% confidence interval
      try(lower <- nest-(1.96*cap.fit$results[position[m],2]),silent=TRUE)
      try(upper <- nest+(1.96*cap.fit$results[position[m],2]),silent=TRUE)
      #Estimate of p
      try(pest <- eval(parse(text=paste("cap.fit$parameters$",mnames[m],"[,2]", sep=""))),silent=TRUE)
      
      #Plut into output
      try(output[j,2,m,i] <- nest, silent=TRUE)
      try(output[j,3,m,i] <- se,silent=TRUE)
      try(output[j,4,m,i] <- pest, silent=TRUE)
    
      #Calculate bias
      try(c.bias <- 100*(nest-Nvec[j])/Nvec[j],silent=TRUE)
      try(output[j,6,m,i] <- c.bias, silent=TRUE)
     
      #Test if true N is in confidence interval
      if(Nvec[j]>=lower){
       if(Nvec[j]<=upper){output[j,7,m,i]=1}
        else{output[j,7,m,i]=0}}
      else{output[j,7,m,i]=0}
      
      #Correct unrealistic estimates of N to NA (100 is cutoff for now)
      if(output[j,2,m,i]>100){output[j,2:7,m,i] = NA}
    }}
    
    #Plug in 0s if no captures
    if(sum(site.data)==0){output[j,2:7,m,i] = 0}
    
    #Move BUGS results to output file
    bugs.out <- fit$sims.list$N[,j]
    bugs.p <- fit$sims.list$p[,j]
    nhat <- mean(bugs.out)
    sdev <- sd(bugs.out)
    phat <- mean(bugs.p)
    sdevp <- sd(bugs.p)
    interval <- quantile(bugs.out, c(0.025,0.975))
    output[j,2,1,i] <- nhat
    output[j,3,1,i] <- sdev
    output[j,4,1,i] <- phat
    output[j,5,1,i] <- sdevp    
    
    #Calculate bias for Nmixture model
    bias <- 100*(nhat-Nvec[j])/Nvec[j]
    output[j,6,1,i] = bias
    
    #Test if true N is in credible interval
    if(Nvec[j]>=interval[1]){
      if(Nvec[j]<=interval[2]){output[j,7,1,i]=1}
      else{output[j,7,1,i]=0}}
    else{output[j,7,1,i]=0}      
  }
    
}

return(output)
}