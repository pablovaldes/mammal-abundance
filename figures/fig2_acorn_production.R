##########################################
## Figure 2: Acorn production over time ##
##########################################

#Read in data
count = read.csv('data/hee_production.csv',header=FALSE)
treedata = read.csv('data/hee_treedata.csv',header=TRUE)

#Generate plot variables
nspecies = 2
nyears = 5

#Create empty matrix
graph.data = matrix(NA,nrow=nspecies,ncol=nyears)

#Fill matrix with sum by species/year (0=black oak, 1=white oak)
for (i in 1:nspecies){
  for (j in 1:nyears){
    
    graph.data[i,j] = sum(count[treedata$species==(i-1),j],na.rm=TRUE)
    
  }}

#Plot data
barplot(graph.data, beside=TRUE, names=c("2006","2007","2008","2009","2010"), 
        xlab="Year", ylab="Collected Acorns (Oct-Dec)", 
        col=c('black','gray'), ylim=c(0,600))

legend(1,600, c('Black Oak', 'White Oak'), fill=c('black','gray'), bty='n')  