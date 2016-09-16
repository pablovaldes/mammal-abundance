##############################################
## Figure: Bias comparisons from simulation ##
##############################################

#Compare correlation and absolute bias for different abundance models
#In master's thesis but not final paper in JWM

#Load simulation workspace
load('output/output_sim_Nmix.Rdata')
data = array(data=NA, dim=c(150,7,6,30,5))

data[,,,,1] = results1
data[,,,,2] = results3
data[,,,,3] = results5
data[,,,,4] = results7
data[,,,,5] = results9


#fix placement of ncaps and mknr
for (i in 5:6){
  for (j in 2:2){
    for (k in 1:30){
      data[,2,i,k,j] = data[,1,i,k,j]
      data[,1,i,k,j] = data[,1,1,k,j]
    }
  }
}

#Correlations
meancor = matrix(data=NA,nrow=6,ncol=2)
corhold = array(data=NA, dim=c(6,2,30))

for (i in 1:6){
  for (j in 1:2){
    for (k in 1:30){
      corhold[i,j,k] = cor((data[,1,i,k,j])[!is.na(data[,2,i,k,j])],na.omit(data[,2,i,k,j]))
    }
    meancor[i,j] = mean(corhold[i,j,])
  }
}

#Plot correlations
barplot(meancor, beside=TRUE, col=gray(seq(0.1,0.9,length=6)),ylim=c(0,1.3),
        ylab="Mean Correlation with True Abundance",xlab=c('Detection Probability on a given Trapping Occasion'),names=c('HBT','Constant 0.5'),
        main = "Abundance Model Comparison")
legend(1,1.3,c('NMixture',expression(M[0]),expression(Chao[h]),expression(Chao[th]),"MKNR",'Ncaps'),
       fill=gray(seq(0.1,0.9,length=6)))

#Absolute bias
meanbias = matrix(data=NA,nrow=4,ncol=5)
biashold = array(data=NA, dim=c(4,5,30,150))

for (i in 1:4){
  for (j in 1:5){
    for (k in 1:30){
      for (l in 1:150){
      biashold[i,j,k,l] = abs(data[l,1,i,k,j]-data[l,2,i,k,j])
    }}
    meanbias[i,j] = mean(biashold[i,j,,],na.rm=TRUE)
  }
}

#Plot absolute bias
barplot(meanbias, beside=TRUE, col=gray(seq(0.1,0.9,length=4)),
        ylab="Mean Absolute Bias",xlab=c('Detection Probability on a given Trapping Occasion'),names=c(0.1,0.3,0.5,0.7,0.9),
        main = "Abundance Model Comparison")
legend(3,10,c('NMixture',expression(M[0]),expression(Chao[h]),expression(Chao[th])),
       fill=gray(seq(0.1,0.9,length=4)))