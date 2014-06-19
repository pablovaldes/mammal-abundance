
model{
  
  #Constant priors
  
  p0 ~ dunif(0,1) #probability of capture with infinite exposure (baseline)
  sigma ~ dunif(0,10) #scale parameter for individual exposure (normal kernel)
  psi ~ dunif(0,1) #zero-inflation parameter for augmentation (sample inclusion probability)
  
  
  #Model kernel
  
  for (i in 1:(nind+nzero)){ #Iterate over individuals + augmented individuals
    
    #Individual priors
    
    Sx[i] ~ dunif(xlower,xupper) #x value of individual i activity center, parameters are constants
    Sy[i] ~ dunif(ylower,yupper) #y value 
    z[i] ~ dbern(psi) #Determine if individual i is included in sample (for augmentation)
    
    for (j in 1:ntraps){ #Iterate over traps in the grid 
      
      D2[i,j] <- pow(Sx[i] - x[j,1],2) + pow(Sy[i]-x[j,2],2) 
      #Squared distance between the location of the activity center of individual i (Sx,Sy),est. parameter
      #and trap j (x[,]),input data
      
      K[i,j] <- exp(-D2[i,j]/sigma)
      #Exposure of individual i to trap j (normal kernel)
      
      gamma[i,j] <- K[i,j]/E[i]
      #Exposure of individual i to trap j divided by the total exposure of individual i E[i]
      #E[i] is estimated below
      #That is, the conditional probability of capture in trap j given that individual i is captured; sum to 1
      
      dprob[i,j] <- Paug[i]*gamma[i,j] 
      #Unconditional probability of individual i in trap j (given that i was sampled; see below)
      
      cp[i,j+1] <- dprob[i,j]
      #Fill in final matrix of capture probs at each trap; cp[,1] left blank for instance of no capture
    }  
    
    E[i] <- sum(K[i,]) #total exposure of individual i to trapping grid (used above)
    
    p[i] <- p0*exp(1/-E[i]) #probability of capture of individual i (not trap-specific)
    
    Paug[i] <- p[i]*z[i] #Capture probability of individual i (0 by default if not sampled)
    
    nocap[i] <- 1-Paug[i] #Probability individual i is not captured (1 by default if not sampled)
    
    cp[i,1] <- nocap[i] #Fill in first slot of capture prob matrix with non-capture probability
    
    for (k in 1:noccs){ #Iterate over the number of trapping occasions
      
      h[i,k] ~ dcat(cp[i,1:(ntraps+1)])
      #h is the input data  matrix, containing the # of the trap individual i was captured at in time k,  
      #or a 1 if individual i was not captured at time k. Given that this is a Sherman trap
      #array, individuals can be captured at most a single time during occasion k.
      #Probability of capture in a given trap is modeled as a multinomial distribution (dcat in WinBUGS)
      
    }
    
  }
  
  #Derived parameters
  
  #Population size
  
  N <- sum(z[1:(nind+nzero)]) #Total number of individuals on the grid, observed and unobserved
  
  #Effective trap area
  
  for (i in 1:(ngridcells)){  #Iterate over an arbitrary number of grid cells of a certain size
    #(i.e., an arbitrary activity center s; should be bigger than actual grid)
    
    for (j in 1:ntraps){ #Iterate over each trap in the grid to calculate exposure of cell s
      
      D_grid2[i,j] <- pow(grid_xy[i,1] - x[j,1],2) + pow(grid_xy[i,2]-x[j,2],2) 
      #grid_xy are supplied coordinates of grid cell i; x[] are locations of actual traps
      #Distance is calculated as above
      
      K_grid[i,j] <- exp(-D_grid2[i,j]/sigma)
      
      #Exposure of grid cell i to trap j
      
    }
    
    E_grid[i] <- sum(K_grid[i,]) #Total exposure of grid cell i
    
    p_grid[i] <- p0*exp(1/(-E_grid[i])) #prob of exposure of grid cell i ###ADDED
    
    E_grid_temp[i] <- 1 - p_grid[i] #Temporary step to simplify calculation ##CHANGED
    
    E_grid_total[i] <- 1 - pow(E_grid_temp[i],noccs) #Probability of exposure for i in 'noccs' occasions
    
    grid_eff_area[i] <- E_grid_total[i]*cell_area #Probability multiplied by area of cell (supplied data)
  }
  
  A_eff <- sum(grid_eff_area[]) #All cells are summed to generate total area; in same units as supplied cell area
  
  #End model specification  
}
