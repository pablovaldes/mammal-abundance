######################################################
## N-mixture Analysis of HEE Small Mammal Abundance ##
######################################################

#Import capture data

#Read in trap.parse() function script
source('function_trap_parse.R')

#Arguments for trap.parse()

#Species included in analysis
species.list = c('EACH','WFMO','SHSH','PIVO')

#Missing grids in 2010
missing.grids = c(1001,1002,1105,1113,1121,1125,
                  1703,1705,1716,1820,1821,1830)

#Only looking at traps inside opening
traps = c('A1','A2','A3','A4','A5','A6',
          'B1','B2','B3','B4','B5','B6',
          'C1','C2','C3','C4','C5','C6')

#Generate occ.list dataframe
EACH = c('1P','2P','3P','4P','5P',NA,NA,NA,NA)
WFMO = c('2A','3A','4A','5A',NA,NA,NA,NA,NA)
SHSH = c('1P','2A','2P','3A','3P','4A','4P','5A','5P')
PIVO = c('1P','2A','2P','3A','3P','4A','4P','5A','5P')
#bind together
occ.list = cbind(EACH,WFMO,SHSH,PIVO)

nyears = 5
nspecies = length(species.list)
nsites = 32

#Create data array to fill
counts.raw = array(NA, dim=c(nsites,dim(occ.list)[1],nspecies,nyears))

#Run function for each year (this will take a while!)
counts.raw[,,,1] = trap.parse("data/2007_captures_scrub.csv", traps=traps,
                          checks=occ.list, species=species.list,output="NmixCounts")$NmixCounts
counts.raw[,,,2] = trap.parse("data/2008_captures_scrub.csv", traps=traps,
                          checks=occ.list, species=species.list,output="NmixCounts")$NmixCounts
counts.raw[,,,3] = trap.parse("data/2009_captures_scrub.csv", traps=traps,
                          checks=occ.list, species=species.list,output="NmixCounts")$NmixCounts
counts.raw[,,,4] = trap.parse("data/2010_captures_scrub.csv", missing.grids = missing.grids, traps=traps,
                          checks=occ.list, species=species.list,output="NmixCounts")$NmixCounts
counts.raw[,,,5] = trap.parse("data/2011_captures_scrub.csv", traps=traps,
                          checks=occ.list, species=species.list,output="NmixCounts")$NmixCounts

#Collapse everything into a 5-sample-occasion array

counts = array(NA,dim=c(nsites,5,nspecies,nyears))
#EACH - afternoon counts
counts[,,1,] = counts.raw[,1:5,1,]
#WFMO - morning counts
counts[,2:5,2,] = counts.raw[,1:4,2,]
#SHSH + PIVO - entire day counts (both morning and night checks)
counts[,1,3,] = counts.raw[,1,3,]
counts[,1,4,] = counts.raw[,1,4,]
index = 2
for (i in 2:5){
counts[,i,3,] = counts.raw[,index,3,] + counts.raw[,(index+1),3,]
counts[,i,4,] = counts.raw[,index,4,] + counts.raw[,(index+1),4,]
index = index+2
}

#Convert into density to account for unequal sample areas (only 0.4 ha patch cuts post harvest are scaled)
#Based on SCR output
#Values average for shrews/voles

counts[c(5,6,21,23,24,25),,1,3:5] = as.integer(counts[c(5,6,21,23,24,25),,1,3:5] / 0.6936)
counts[c(5,6,21,23,24,25),,2,3:5] = as.integer(counts[c(5,6,21,23,24,25),,2,3:5] / 0.6502)
counts[c(5,6,21,23,24,25),,3:4,3:5] = as.integer(counts[c(5,6,21,23,24,25),,3:4,3:5] / 0.6719)

##########################################################################

#Generate Indexes

#Was site sampled in 2010? 4 if yes, 5 if no
first = c(5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,4,4,5,5,5,4,4,4,4)
#What day was species i first sampled? 2 for mice, 1 for others
sindex = c(1,2,1,1)
#Initial values for N array: max reported count (required for JAGS)
Nst = array(data=NA,dim=c(nsites,nspecies,nyears))
for (i in 1:nspecies){
  for (j in 1:nsites){
    for (t in 1:3){
      Nst[j,i,t] = max(counts[j,1:5,i,t],na.rm=TRUE)+1
    }
    for (t in first[j]:5){
      Nst[j,i,t] = max(counts[j,1:5,i,t],na.rm=TRUE)+1
    }
}}

##########################################################################

#Input Detection Covariates

#Read in grid-level data
grid.data = read.csv('data/hee_grid_covs.csv',header=TRUE)

#Temperature and Julian day
temp.raw = as.matrix(grid.data[,1:5])
temp = (temp.raw - mean(temp.raw, na.rm=TRUE))/sd(temp.raw,na.rm=TRUE)
jd.raw = as.matrix(grid.data[,6:10])
jd = (jd.raw - mean(jd.raw, na.rm=TRUE))/sd(jd.raw,na.rm=TRUE)

#Effort using trap.parse()
eff.raw = array(data=NA, dim=c(32,4,5))
eff.raw[,,1] = trap.parse('data/2007_captures_scrub.csv', species=species.list, traps=traps,
                          checks=occ.list, output="TrapNights")$TrapNights
eff.raw[,,2] = trap.parse('data/2008_captures_scrub.csv', species=species.list, traps=traps,
                          checks=occ.list, output="TrapNights")$TrapNights
eff.raw[,,3] = trap.parse('data/2009_captures_scrub.csv', species=species.list, traps=traps,
                          checks=occ.list, output="TrapNights")$TrapNights
eff.raw[,,4] = trap.parse('data/2010_captures_scrub.csv', species=species.list, traps=traps,
                          checks=occ.list, missing.grids=missing.grids, output="TrapNights")$TrapNights
eff.raw[,,5] = trap.parse('data/2011_captures_scrub.csv', species=species.list, traps=traps,
                          checks=occ.list, output="TrapNights")$TrapNights
eff = (eff.raw - mean(eff.raw, na.rm=TRUE))/sd(eff.raw,na.rm=TRUE)

##########################################################################

#Input Abundance Covariates

#Treatment Index (grid x treatment: clear, shelter, patch(all), 0.4/1 ha patch, 2 ha patch)
treat = array(data=NA, dim=c(32,5,5))
treat[,,1:2] = as.matrix(grid.data[,17:21])
treat[,,3:5] = as.matrix(grid.data[,22:26])

#Aspect
aspect = grid.data$aspect

#Mast
mast.raw = as.matrix(grid.data[,12:16]) 
mast = (mast.raw - mean(as.numeric(mast.raw)))/sd(as.numeric(mast.raw))

###########################################################################

#JAGS arguments

#Bind all data together

data <- list("counts","sindex","first","nspecies","nsites","jd",
             "temp","eff","aspect","treat","mast"
)

# Parameters to save
params <- c(
  'mu.lam','mu.p', 'alpha.lam', 'alpha.p',	
  #Detection parameters
  'beta.temp','beta.jd','beta.eff',
  'pbeta.c.cl', 'pbeta.c.sh', 'pbeta.04', 'pbeta.2',
  #Abundance parameters
  "beta.aspect", "beta.c.cl","beta.c.sh",
  "beta.mast", "beta.04", "beta.2",
  #Derived parameters
  "pre.control.mean","pre.clear.mean","pre.shelter.mean","pre.04.mean", "pre.2.mean",
  "post.control.mean","post.clear.mean","post.shelter.mean","post.04.mean", "post.2.mean",
  "pre.control.p","pre.clear.p","pre.shelter.p","pre.04.p", "pre.2.p",
  "post.control.p","post.clear.p","post.shelter.p","post.04.p", "post.2.p", "p.all",
  "yearly.control","yearly.clear","yearly.shelter","yearly.patch04","yearly.patch2"
)

# Initial values function
inits <- function(){  
  list(N=Nst      
  )}

#Model file
modFile = "models/model_abundance_Nmix.R"  

###########################################################################

#Run analysis in JAGS
library(jagsUI)
abundance.out = jags(data, inits, params, model.file=modFile, 
                    n.chains=4, n.iter=15000, n.burnin=11000, n.thin=40, parallel=TRUE)

#Save output
save(abundance.out,file="output_hee_abundance_Nmix.Rdata")

###########################################################################

#Validation of Nmix method via simulation

#Load simulation functions
source('sim_abundance.R')

#Generate true abundance vector
nsites <- 150
lambda <- 10
Nvec <- rpois(n=nsites,lambda=lambda)

#Example simulation run with default values
results5 = fitsim(1,Nvec,p=0.5)

save(results5,file='output_sim_Nmix.Rdata')
