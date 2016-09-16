####################################################
## Non-parametric Analysis of HEE Microsite Data, ##
## Between treatments and years                   ##
####################################################

#Read in data (summarized CSV)
micro.data = read.csv('data/hee_microsite_means.csv',header=TRUE)
micro.data$grid = as.factor(micro.data$grid)
micro.data$year = as.factor(micro.data$year)

####################################################################

#Overall differences between years, by treatment
#Kruskal-Wallis test (analogue to ANOVA)

#Summarize data by treatment
micro.data.clear = micro.data[micro.data$treat==2,]
micro.data.patch04 = micro.data[micro.data$treat==3,]
micro.data.patch2 = micro.data[micro.data$treat==4,]
micro.data.shelter = micro.data[micro.data$treat==5,]
micro.data.control = micro.data[micro.data$treat==1,]

#Change data in this line to test different treatments
kruskal.test(herb ~ year, data=micro.data.clear)


####################################################################

#Post-hoc Multiple comparisons between years

#Non-parametric multiple comparisons
#Siegel and Castellan 1988 Nonparametric tests for the behavioral sciences

#| RBari - RBarj | >= Z*Sqrt[ (N*(N+1)/12)*(1/ni + 1/nj) ]

#where RBari, RBarj, ni, nj are the mean of the ranks and the sample sizes 
#associated with the i-th and j-th groups. N is the total sample size and Z 
#is the critical value from the standard normal curve.

library(pgirmess)
kruskalmc(micro.data.clear$litter, micro.data.clear$year)
kruskalmc(micro.data.clear$herb, micro.data.clear$year)
kruskalmc(micro.data.clear$wood, micro.data.clear$year)
kruskalmc(micro.data.clear$cwd, micro.data.clear$year)

#Patch
kruskalmc(micro.data.patch04$litter, micro.data.patch04$year)
kruskalmc(micro.data.patch04$herb, micro.data.patch04$year)
kruskalmc(micro.data.patch04$wood, micro.data.patch04$year)
kruskalmc(micro.data.patch04$cwd, micro.data.patch04$year)

kruskalmc(micro.data.patch2$litter, micro.data.patch2$year)
kruskalmc(micro.data.patch2$herb, micro.data.patch2$year)
kruskalmc(micro.data.patch2$wood, micro.data.patch2$year)
kruskalmc(micro.data.patch2$cwd, micro.data.patch2$year)

#Shelter
kruskalmc(micro.data.shelter$litter, micro.data.shelter$year)
kruskalmc(micro.data.shelter$herb, micro.data.shelter$year)
kruskalmc(micro.data.shelter$wood, micro.data.shelter$year)
kruskalmc(micro.data.shelter$cwd, micro.data.shelter$year)

#Control
kruskalmc(micro.data.control$litter, micro.data.control$year)
kruskalmc(micro.data.control$herb, micro.data.control$year)
kruskalmc(micro.data.control$wood, micro.data.control$year)
kruskalmc(micro.data.control$cwd, micro.data.control$year)
