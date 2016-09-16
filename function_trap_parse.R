############################################################
## Function to convert HEE trap data into various formats ##
############################################################

#Arguments
####################################################################

# file = Filepath to a CSV file with a single year of trapping data, 
#        with each row representing one "trapping event" 
#        (could be a capture, or a trap death, or a sprung trap, etc.). 
#
#     The file should contain at minimum the following column headings:
#     1. Grid – trapping grid number or name
#     2. Trap – trap numbers or locations or codes
#     3. Species – names or codes of species
#     4. Checks – trap check identifier (e.g. 5A for day 5 AM)
#     5. ID – individual tag ID numbers if animals were tagged
#     6. Fate – trap fate (if trap-night calculation is required). 
#     Disturbed and sprung traps should be given Fate values of 7 and 8, respectively.
#
# grids = If supplied, a vector of grid names to include in the parse. 
#         If not supplied, all grids included.
#
# missing.grids = A vector of grid names not sampled in a given year. 
#                 They will be added into the dataset as ordered anyway, 
#                 and filled with NAs (to preserve structure for multiple years of data). 
#                 Should not be specified at the same time as grids option.
#
# traps = Vector of trap numbers/locations/codes. If not supplied, all are included.
#
# species = Vector of species names or codes to include. If not supplied, all are included.
#
# checks = A vector of which trap check occasions to include in the parse. 
#          A matrix (where each column is a vector of checks to include for a given species, 
#          with column header matching the species name in the species vector) can also be supplied. 
#          The matrix must be populated with NAs at the end of each column,
#          if not every species is included in the same number of checks.

#Output
#########################################################################

# Character vector of options for which results to calculate and return. 
# More than one can be calculated and returned at once. Possibilities are:
# 
# 1. TrapOccupancy – trap-level encounter histories, dimension trap*occasion*species*grid
# 2. TrapCounts – count of species caught in a particular trap (binomial, N=# of occasions)
# 3. NmixCounts – capture counts by grid/trapping occasion/species, grid*occasion*species
# 4. SiteOccupancy – site-level encounter histories, grid*occasion*species
# 5. IndHistories – individual animal-level encounter histories (for tagged animals only), 
#    individual*occasion*species*grid
# 6. SCR – spatial capture-recapture data: instead of encounter history with 0 for individual not captured
#    and 1 if it was, returns 1 if not captured and the trap # caught at if it was.
# 7. TrapNights – effective trap nights (removing disturbed traps), grid*species

#Begin function
trap.parse = function(file, grids=FALSE, missing.grids=FALSE, traps=FALSE, 
                      species=FALSE, checks=FALSE,output){
  
  #Read in datafile. Should be CSV format only.
  data = as.data.frame(read.csv(file,header=TRUE))
  
  #############################
  #If defaults are selected, make lists of grids, traps, etc.
  
  if(grids[1]==FALSE){grid.list = sort(unique(data$Grid))}else{grid.list=grids}
  
  if(missing.grids[1]!=FALSE){
    #If grids that actually have data specified as missing, slice out from list...
    add = missing.grids[!missing.grids%in%grid.list]
    #and print warning
    if(length(add)!=length(missing.grids)){print('Data found for some grids specified missing')}
    grid.list = sort(c(grid.list,add))}
  
  #Get all values for traps, species, checks if not provided
  if(traps[1]==FALSE){trap.list = sort(as.character(unique(data$Trap)))}else{trap.list=traps}
  
  if(species[1]==FALSE){species.list = as.character(unique(data$Species))}else{species.list=species}
  
  if(checks[1]==FALSE){occ.list = as.data.frame(rep(sort(unique(data$Check)),length(species.list))
                                         ,ncol=length(species.list))
                       names(occ.list) = species.list
  }else{occ.list=checks} 
  
  ##############################
  #Set additional index variables based on above
  
  ngrids = length(grid.list)
  ntraplocs = length(trap.list)
  nspecies = length(species.list)
  noccs = dim(occ.list)[1]
  
  ##############################
  
  #Run this code if either trap-level occupancy data or Nmixture model data are requested.
  
  if('TrapOccupancy'%in%output || 'NmixCounts'%in%output || 'SiteOccupancy'%in%output || 'TrapCounts'%in%output){    
  
  #Set up empty trap-level encounter history (non-tagged data)
  tr.encounter = array(NA,dim=c(ntraplocs,noccs,nspecies,ngrids))
  
  #Begin base-level encounter history generation (trap-level)
  #Iterate through all combinations of grid/trap/species/trapping occasion
  for (i in 1:ngrids){
    for (j in 1:ntraplocs){
      for (k in 1:nspecies){
        for (l in 1:length(na.omit(as.character(occ.list[,species.list[k]])))){     
          
          #Don't do anything if the grid is in the missing grids list (leave NA)
          if(!grid.list[i]%in%missing.grids){
          
          #Create new dataframe of only rows in given combination of grid/trap/species/occ
          hold = filter(data,Grid==grid.list[i]
                        ,Trap==trap.list[j]
                        ,Species==species.list[k]
                        ,Check==occ.list[l,species.list[k]])
          
          #Check to see if the trapping occasion is in the new dataframe hold[]
          if (occ.list[l,species.list[k]]%in%hold$Check){
            #if it is, change the NA or 0 to a 1 for present
            if(is.na(tr.encounter[j,l,k,i])){tr.encounter[j,l,k,i] = 1}
            else if(tr.encounter[j,l,k,i]==0){tr.encounter[j,l,k,i] = 1}
          }else{
            #if not, change to a 0 for not present
            tr.encounter[j,l,k,i] = 0}
          
          }
    }}}}
  }

  ##################################################
  
  #Run this code if trap-level counts (binomial) are requested
  
  #Sum over each traploc x species x grid combination
  
  if('TrapCounts'%in%output){
    tr.counts = array(data=NA,dim=c(ntraplocs,nspecies,ngrids))
    for (i in 1:ntraplocs){
      for (j in 1:nspecies){
        for (k in 1:ngrids){
          #Don't do anything if the grid is in the missing grids list (leave NA)
          if(!grid.list[k]%in%missing.grids){
          tr.counts[i,j,k] = sum(tr.encounter[i,,j,k], na.rm=TRUE)
          }
    }}}
       
  }
  
  ##################################################
  
  #Run this code if individual encounter histories (tagged animals)
  #are requested (basic or SCR)
  
  if('IndHistories'%in%output || 'SCR'%in%output){
  
  #Figure out how many rows the array should have (max # of individuals caught in a 
  #single one of the grids sampled)
    
  maxcount = vector(length=(ngrids*nspecies))
  index=1
  #Iterate through grids/species and put the number of individuals in maxcount vector
  for (i in 1:ngrids){                               
    for (j in 1:nspecies){
      hold = filter(data,Grid==grid.list[i]
                    ,Species==species.list[j])  
      maxcount[index] = length(unique(hold$ID[!is.na(hold$ID)]))
      index=index+1
    }}
  
  #Set the number of rows in the array to the maximum observed species/grid combo
  #To insure all other grid/species combos will fit in the arraya
  nid = max(maxcount,na.rm=TRUE)
  
  #Set up empty individual encounter array. Individual x occasion x species x grid
  #Separate encounter array for SCR.
  if('IndHistories'%in%output){ind.encounter = array(data=NA, dim=c(nid,noccs,nspecies,ngrids))}
  if('SCR'%in%output){scr.encounter = array(data=NA,dim=c(nid,noccs,nspecies,ngrids))}
  
  #Iterate through grids and species and select data
  for (i in 1:ngrids){                               
    for (j in 1:nspecies){
      hold = filter(data,Grid==grid.list[i]
                    ,Species==species.list[j])
      
      #Generate list of unique tag IDs in grid x species combo
      id.list = as.vector(na.omit(as.numeric(as.vector(unique(hold$ID)))))
      
      #Iterate through ID list to generate encounter history for each individual
      for (k in 1:length(id.list)){
        
        hold.id = filter(hold,ID==id.list[k])
        
        for (l in 1:length(na.omit(as.character(occ.list[,species.list[j]])))){
          
          #If individual was captured in a particular occasion l
          if (occ.list[l,species.list[j]]%in%hold.id$Check){ 
            
            #Do this if individual encounter history needed
            if('IndHistories'%in%output){
            
              #Set encounter as 1
              if(is.na(ind.encounter[k,l,j,i])){ind.encounter[k,l,j,i] = 1}
              else if(ind.encounter[k,l,j,i] == 0){ind.encounter[k,l,j,i] = 1}
            }
            
            #Do this if SCR encounter history needed
            if('SCR'%in%output){
              if(is.na(scr.encounter[k,l,j,i])){
                
                scr.encounter[k,l,j,i]= (which(trap.list==hold.id$Trap[
                  which(hold.id$Check==occ.list[l,species.list[j]])])+1)
              
              }
              
              
              
            }
              #Otherwise set to 'not captured' code (0 for basic, 1 for SCR)
              }else {
                
               if('IndHistories'%in%output){ind.encounter[k,l,j,i] = 0} 
               if('SCR'%in%output){scr.encounter[k,l,j,i] = 1}  
                
              }
          
  }}}}
  }
  
  ##################################################
  
  #Run this code if site (grid) - level occupancy is requested
  
  if('SiteOccupancy'%in%output){
    #Grid-level occupancy - grid x occasion x species
    
    #Set up empty array
    siteocc = array(NA,dim=c(length(grid.list),noccs,length(species.list)))
    
    #Begin iterating over grids, species, and occasions
    for (i in 1:ngrids){
      for (j in 1:nspecies){
        for (k in 1:length(na.omit(occ.list[,species.list[j]]))){  
          #Set site occupancy to 1 if at least one animal is detected
          if(sum(tr.encounter[,k,j,i]>0)){siteocc[i,k,j] = 1
          #Otherwise set to 0
          }else{siteocc[i,k,j] = 0}
        
    }}}  
  }
  
  ##################################################
  
  #Run this code if Nmixture output is requested.
  
  if('NmixCounts'%in%output){
    #Grid-level abundance - grid x occasion x species
    
    #Set up empty array
    counts = array(NA,dim=c(length(grid.list),dim(occ.list)[1],length(species.list)))         
    
    #Begin iterating
    #Sum over all traps in a given grid/occasion/species combination
    for (i in 1:ngrids){                                                  
      for (j in 1:nspecies){
        for (k in 1:length(na.omit(occ.list[,species.list[j]]))){ 
          #Sum all individuals of species j in grid i on occasion k
          counts[i,k,j] = sum(tr.encounter[,k,j,i])
        }}} 
  }
  
  ##################################################
  
  #Run this code if Trap nights output is required
  
  if('TrapNights'%in%output){
    
    #Set up empty trap nights array (grids x species)
    eff = array(NA, dim=c(length(grid.list),nspecies))
    
    #Iterate through all selected grids
    for (i in 1:ngrids){ 
      #Ignore if grid is missing
      if(!grid.list[i]%in%missing.grids){
      
      #Select data entries for grid i
      hold1 = filter(data,Grid==grid.list[i])
      
      #Iterate over species for grid i
      for (j in 1:nspecies){
        
        #Select data entries in sample occasions for species j
        hold2 = hold1[hold1$Check%in%na.omit(occ.list[,species.list[j]]),]  
        
        #Sum all traps in these data entries that are disturbed or sprung
        d = sum(hold2$Fate==7)+sum(hold2$Fate==8)
        #If NA for some reason, set d to 0
        if(is.na(d)){d=0}
        
        #Calculate effort: # of traps * number of sample occasions for species j - number disturbed/sprung
        eff[i,j] <- ntraplocs*length(na.omit(occ.list[,species.list[j]])) - d
        
    }}}  
  }
  
  ##################################################
  
  #Return output
  result = list("TrapOccupancy"=NULL,"TrapCounts"=NULL,"SiteOccupancy"=NULL,
                "NmixCounts"=NULL,"IndHistories"=NULL,"SCR"=NULL,"TrapNights"=NULL)
  
  if("TrapOccupancy"%in%output){result$TrapOccupancy = tr.encounter}
  if("TrapCounts"%in%output){result$TrapCounts = tr.counts}
  if('SiteOccupancy'%in%output){result$SiteOccupancy = siteocc}
  if('NmixCounts'%in%output){result$NmixCounts = counts}
  if('IndHistories'%in%output){result$IndHistories = ind.encounter}
  if('SCR'%in%output){result$SCR = scr.encounter}
  if('TrapNights'%in%output){result$TrapNights = eff}
  
  return(result)
  
}
  
##End script