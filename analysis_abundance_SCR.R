################################################
## Spatially Explicit Capture-Recapture (SCR) ##
## Analysis of HEE MRR Data                   ##
################################################

#Input data

#Get trap.parse() function
source('function_trap_parse.R')

#Specify arguments for trap.parse()

species.list = c('EACH','WFMO')

#List of grids where animals were marked
#Grid 1202 was replaced with 1504 in 2011
grid.list08 = c('1202','1312','1502','1728','1803','1907')
grid.list11 = c('1504','1312','1502','1728','1803','1907')

#Generate occ.list dataframe
EACH = c('1P','2P','3P','4P','5P')
WFMO = c('2A','3A','4A','5A',NA)

#bind together
occ.list = cbind(EACH,WFMO)

#Indexes
nspecies = length(species.list)
noccs = c(5,4)
nsites = length(grid.list08)
nyears = 3

#Parse data
scr1 = trap.parse('data/2008_captures_scrub.csv',species=species.list,check=occ.list,
                 grids=grid.list08,output='SCR')$SCR
scr2 = trap.parse('data/2010_captures_scrub.csv',species=species.list,check=occ.list,
                 grids=grid.list08,output='SCR')$SCR
scr3 = trap.parse('data/2011_captures_scrub.csv',species=species.list,check=occ.list,
                 grids=grid.list11,output='SCR')$SCR

###############################################################################

#Break into separate species x grid size datasets

##########################EACH
each.08 = scr1[,,1,]
each.10 = scr2[,,1,]
each.11 = scr3[,,1,]

#Small grids (0.4 ha)
each.small = rbind( each.08[apply(each.08[,,4],1,function(x)any(!is.na(x))),,4],
                    each.08[apply(each.08[,,5],1,function(x)any(!is.na(x))),,5],
                    each.10[apply(each.10[,,4],1,function(x)any(!is.na(x))),,4],
                    each.10[apply(each.10[,,5],1,function(x)any(!is.na(x))),,5],
                    each.11[apply(each.11[,,4],1,function(x)any(!is.na(x))),,4],
                    each.11[apply(each.11[,,5],1,function(x)any(!is.na(x))),,5])

#Large grids (1 ha)
each.large = rbind( each.08[apply(each.08[,,1],1,function(x)any(!is.na(x))),,1],
                    each.08[apply(each.08[,,2],1,function(x)any(!is.na(x))),,2],
                    each.08[apply(each.08[,,3],1,function(x)any(!is.na(x))),,3],
                    each.08[apply(each.08[,,6],1,function(x)any(!is.na(x))),,6],
                    each.10[apply(each.10[,,1],1,function(x)any(!is.na(x))),,1],
                    each.10[apply(each.10[,,2],1,function(x)any(!is.na(x))),,2],
                    each.10[apply(each.10[,,3],1,function(x)any(!is.na(x))),,3],
                    each.10[apply(each.10[,,6],1,function(x)any(!is.na(x))),,6],
                    each.11[apply(each.11[,,1],1,function(x)any(!is.na(x))),,1],
                    each.11[apply(each.11[,,2],1,function(x)any(!is.na(x))),,2],
                    each.11[apply(each.11[,,3],1,function(x)any(!is.na(x))),,3],
                    each.11[apply(each.11[,,6],1,function(x)any(!is.na(x))),,6])

#############################WFMO
wfmo.08 = scr1[,,2,]
wfmo.10 = scr2[,,2,]
wfmo.11 = scr3[,,2,]

#Small grids (0.4 ha)
wfmo.small = rbind( wfmo.08[apply(wfmo.08[,,4],1,function(x)any(!is.na(x))),,4],
                    wfmo.08[apply(wfmo.08[,,5],1,function(x)any(!is.na(x))),,5],
                    wfmo.10[apply(wfmo.10[,,4],1,function(x)any(!is.na(x))),,4],
                    wfmo.10[apply(wfmo.10[,,5],1,function(x)any(!is.na(x))),,5],
                    wfmo.11[apply(wfmo.11[,,4],1,function(x)any(!is.na(x))),,4],
                    wfmo.11[apply(wfmo.11[,,5],1,function(x)any(!is.na(x))),,5])

#Large grids (1 ha)
wfmo.large = rbind( wfmo.08[apply(wfmo.08[,,1],1,function(x)any(!is.na(x))),,1],
                    wfmo.08[apply(wfmo.08[,,2],1,function(x)any(!is.na(x))),,2],
                    wfmo.08[apply(wfmo.08[,,3],1,function(x)any(!is.na(x))),,3],
                    wfmo.08[apply(wfmo.08[,,6],1,function(x)any(!is.na(x))),,6],
                    wfmo.10[apply(wfmo.10[,,1],1,function(x)any(!is.na(x))),,1],
                    wfmo.10[apply(wfmo.10[,,2],1,function(x)any(!is.na(x))),,2],
                    wfmo.10[apply(wfmo.10[,,3],1,function(x)any(!is.na(x))),,3],
                    wfmo.10[apply(wfmo.10[,,6],1,function(x)any(!is.na(x))),,6],
                    wfmo.11[apply(wfmo.11[,,1],1,function(x)any(!is.na(x))),,1],
                    wfmo.11[apply(wfmo.11[,,2],1,function(x)any(!is.na(x))),,2],
                    wfmo.11[apply(wfmo.11[,,3],1,function(x)any(!is.na(x))),,3],
                    wfmo.11[apply(wfmo.11[,,6],1,function(x)any(!is.na(x))),,6])

#####################################################################################

#Input trap location data (based on HEE grid sizes)

ntraps = 27

###########################

#Generate trap locations for large grids
x.temp = c(36,34,32,30,28,26,21,16,11,36,34,32,30,28,26,21,16,11,36,34,32,30,28,26,21,16,11)
y.temp = c(rep(11,9),rep(15,9),rep(19,9)) 

x.large = cbind(x.temp,y.temp)

#Generate points for large grids
gridx.large = matrix(NA,nrow=1426,ncol=2)
index=1
for (i in 1:46){
  for (j in 1:31){
    gridx.large[index,1]=i
    gridx.large[index,2]=j
    index = index+1
  }
}

####################################

#Generate trap locations for small grids
x.temp = c(rep(c(27,22,17,16,15,14,13,12,11),3))
y.temp = c(rep(11,9),rep(13,9),rep(15,9)) 

x.small = cbind(x.temp,y.temp)

#Points for small grids
gridx.small = matrix(NA,nrow=864,ncol=2)
index=1
for (i in 1:36){
  for (j in 1:24){
    gridx.small[index,1]=i
    gridx.small[index,2]=j
    index = index+1
  }
}

######################################

#Package separate datasets by species and grid size (4 total)


########WFMO LARGE##################
nind = dim(wfmo.large)[1]
nzero = 150-nind
h = rbind(wfmo.large,matrix(data=1,nrow=nzero,ncol=5))
x = x.large
noccs = 4
grid_xy = gridx.large
ngridcells = dim(grid_xy)[1] #(46 by 31; 100m buffers on all sides of grid)
cell_area = 100 #10 by 10 m cells
xlower = 1
xupper = 46
ylower = 1
yupper = 31

#For jags parallel

wfmo_large_data=c('h','nind','nzero','x','noccs','ngridcells','ntraps','grid_xy','cell_area',
                     'xlower','xupper','ylower','yupper')

#############EACH LARGE##########################
nind = dim(each.large)[1]
nzero = 150-nind
h = rbind(each.large,matrix(data=1,nrow=nzero,ncol=9))
x = x.large
noccs = 9
grid_xy = gridx.large
ngridcells = dim(grid_xy)[1] #(46 by 31; 100m buffers on all sides of grid)
cell_area = 100 #10 by 10 m cells
xlower = 1
xupper = 46
ylower = 1
yupper = 31

each_large_data=list('h','nind','nzero','x','noccs','ngridcells','ntraps','grid_xy','cell_area',
                     'xlower','xupper','ylower','yupper')

##############WFMO SMALL#########################
nind = dim(wfmo.small)[1]
nzero = 150-nind
h = rbind(wfmo.small,matrix(data=1,nrow=nzero,ncol=9))
x = x.small
noccs = 4
grid_xy = gridx.small
ngridcells = dim(grid_xy)[1] #(46 by 31; 100m buffers on all sides of grid)
cell_area = 100 #10 by 10 m cells
xlower = 1
xupper = 36
ylower = 1
yupper = 24

wfmo_small_data=list('h','nind','nzero','x','noccs','ngridcells','ntraps','grid_xy','cell_area',
                     'xlower','xupper','ylower','yupper')

##############EACH SMALL#########################
nind = dim(each.small)[1]
nzero = 150-nind
h = rbind(each.small,matrix(data=1,nrow=nzero,ncol=9))
x = x.small
noccs = 9
grid_xy = gridx.small
ngridcells = dim(grid_xy)[1] 
cell_area = 100 #10 by 10 m cells
xlower = 1
xupper = 36
ylower = 1
yupper = 24

each_small_data=list('h','nind','nzero','x','noccs','ngridcells','ntraps','grid_xy','cell_area',
                     'xlower','xupper','ylower','yupper')

############################################################################################

#Other JAGS arguments

#Inits for jags
inits <- function(){
  list(p0=0.1,
       psi=0.9)}

#Parameters to save

params = c('sigma', 'psi', 'p0', 'N', 'A_eff')

#Begin model specification
modFile = "models/model_abundance_SCR.R"

#Run JAGS analysis (in parallel)
library(jagsUI)
#change output name and dataset to correspond to each species/size combo
wfmo_large_output = jags(wfmo_large_data, inits, params, model.file=modFile, 
                               n.chains=4, n.iter=10000, n.burnin=1000, n.thin=20,parallel=TRUE)

#output
save(each_small_output,file="output/output_scr_each_small.Rdata")

#Actual estimated sample area sizes
#In this analysis only individuals that were detected at least once
#in the interior 18 traps were included.
mice.large = 3.3388
mice.small = 2.1709
mouse.ratio = 0.6502

each.large = 3.4903
each.small = 2.4210
each.ratio = 0.6936