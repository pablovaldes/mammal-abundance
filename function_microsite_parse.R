########################################################
## Function to parse HEE microsite data and summarize ##
########################################################

#Arguments
####################################################################

# file = Filepath to a CSV file with a single year of trapping data, 
#        with each row representing one "trapping event" 
#        (could be a capture, or a trap death, or a sprung trap, etc.). 
#
#     The file should contain at minimum the following column headings:
#     1. Grid – trapping grid number or name
#     2. Trap – trap numbers or locations or codes
#     3-6. Microhabitat variable names (herb, wood, cwd, litter)
#
# grids = If supplied, a vector of grid names to include in the parse. 
#         If not supplied, all grids included.
#
# var.list = A character vector of the microsite variables to include in the output; 
#            must match microhabitat variable column headings in input file.
#
# traps = Vector of trap numbers/locations/codes. If not supplied, all are included.
#
# missing.grids = A vector of grid names not sampled in a given year. 
#                 They will be added into the dataset as ordered anyway, 
#                 and filled with NAs (to preserve structure for multiple years of data). 
#                 Should not be specified at the same time as grids option.
#

#Output
#########################################################################

# Output is a list with two elements: 
# 
# micro – A formatted array of microsite values, 
#         with dimensions trap locations x microhabitat variables x grids
# 
# micro.mean - Yearly means and standard deviations for each microhabitat variable.

#Begin function
microsite.parse <- function(file, grids=FALSE, traps=FALSE, var.list, missing.grids=FALSE){
  
  #Read in data
  data <- read.csv(file, header=TRUE)
  
  #List of grids trapped and traps to include
  
  if(grids[1]==FALSE){grid.list <- sort(unique(data$grid))}else{grid.list=grids}
  
  if(missing.grids[1]!=FALSE){
    #If grids that actually have data specified as missing, slice out from list...
    add = missing.grids[!missing.grids%in%grid.list]
    #and print warning
    if(length(add)!=length(missing.grids)){print('Data found for some grids specified missing')}
    grid.list = sort(c(grid.list,add))}
  
  if(traps[1]==FALSE){trap.list = sort(as.character(unique(data$trap)))}else{trap.list=traps}
  
  ngrids <- length(grid.list)
  ntraplocs <- length(trap.list)
  nvars <- length(var.list)
  
  #Ordinal to numeric conversion vector (for woody/herbaceous veg)
  conv <- c(0.5,3,15,37.5,62.5,85,97.5)
  
  #Create empty output objects
  micro.mean <- matrix(data=NA, nrow=ngrids, ncol=nvars)
  micro <- array(data=NA,dim=c(ntraplocs,nvars,ngrids))
  
  #Iterate through grids
  for (i in 1:ngrids){
    
    #Work only in grids that are not missing
    if(!grid.list[i]%in%missing.grids){   
    
    #Iterate through variables and trap locations
    for (j in 1:nvars){
      for (k in 1:ntraplocs){    
        
        #Select data
        hold <- filter(data,grid==grid.list[i],trap==trap.list[k])[,var.list[j]]
        
        #Convert from ordinal to numeric if herbaceous or woody vegetation variable
        if(var.list[j] == 'herb' || var.list[j] == 'wood'){
          ord <- as.integer(ceiling(as.numeric(hold)))
          hold <- conv[ord]
        }
        
        #Fill in output object
        micro[k,j,i] <- hold
      }
      
      #Fill in means summary output object
      micro.mean[i,j] <- mean(micro[,j,i],na.rm=TRUE)
    }
  }}
  
  #Return both output objects
  return(list(micro=micro, micro.mean=micro.mean))
}