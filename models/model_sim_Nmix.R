model {
  
  #Priors
  mu.lam ~ dnorm(0,0.01)  # intercept
  mu.p ~ dnorm(0,0.01)
  
  sigma.lam ~ dunif(0,10)  
  tau.lam <- pow(sigma.lam,-2)  
  sigma.p ~ dunif(0,10)  
  tau.p <- pow(sigma.lam,-2)
  
  
  #Model Specification
  
  for (i in 1:nsites){
    alpha.p[i] ~ dnorm(mu.p,tau.p)#I(-5,5) 
    alpha.lam[i] ~ dnorm(mu.lam,tau.lam)#I(-5,5) #Can cause crashes when alpha.lam is small
    N[i] ~ dpois(lambda[i])
    logit(p[i]) <- alpha.p[i]
    log(lambda[i]) <- alpha.lam[i]
    for (j in 1:5){
      counts[i,j] ~ dbin(p[i],N[i])
      
    }}
  
}