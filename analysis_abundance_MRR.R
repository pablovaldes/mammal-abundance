##############################
## MRR analysis of HEE Data ##
##############################

#Load required library
library(Rcapture)

#Get trap.parse() function
source('function_trap_parse.R')

#Arguments for trap.parse()

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

#Read in individual encounter histories for years with marked animals
eh1 = trap.parse('data/2008_captures_scrub.csv',species=species.list,check=occ.list,
                  grids=grid.list08,output='IndHistories')$IndHistories
eh2 = trap.parse('data/2010_captures_scrub.csv',species=species.list,check=occ.list,
                   grids=grid.list08,output='IndHistories')$IndHistories
eh3 = trap.parse('data/2011_captures_scrub.csv',species=species.list,check=occ.list,
                   grids=grid.list11,output='IndHistories')$IndHistories

#Combine into single array
eh = array(data=NA,dim=c(19,noccs,nspecies,nsites,nyears))
eh[1:dim(eh1)[1],,,,1] = eh1
eh[1:dim(eh2)[1],,,,2] = eh2
eh[1:dim(eh3)[1],,,,3] = eh3

#Fix errors
eh[,,2,2,2] = NA
eh[,,2,3,2] = NA
eh[,,2,5,2] = NA

#Run Analyses site x species x year
for(i in 1:nsites){
  for (j in 1:nspecies){
    for (k in 1:nyears){
      #Separate dataset for i,j,k
      hold = na.omit(eh[,1:noccs[j],j,i,k])
      if(dim(hold)[1] == 0){print('no data')
      }else{print(c(j,i,k))
            print(closedp.t(hold))}
  
}}}
