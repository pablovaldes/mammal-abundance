##########################################################
## Figure 5: Comparison of abundance estimation methods ##
##########################################################

ngrids = 32
nspecies = 4
noccs = 5

#Read in abundance estimates

#Total Captures
#Requires 'counts' array from 'analysis_abundance_Nmix.R'
totalcount = array(data=NA, dim=c(32,4,5))
for (i in 1:ngrids){
  for (j in 1:nspecies){
    for (k in 1:noccs){
      totalcount[i,j,k] = sum(counts[i,,j,k],na.rm=TRUE)
}}}

#Effort (Trap Nights)
#Requires 'eff.raw' array from 'analysis_abundance_Nmix.R'

caps100 = (totalcount/eff.raw)*100

#Extract only sites with marked animal data
final.count = matrix(data=NA, nrow=2,ncol=18)
for (j in 1:2){
k=1
for (t in c(2,4)){
  for (i in c(7,9,15,24,25,29)){
    final.count[j,k] = caps100[i,j,t]
    k = k+1
  }}
for (i in c(16,9,15,24,25,29)){
  final.count[j,k] = caps100[i,j,5]
  k = k+1
  }}

####################################

#Chipmunk MRR and Nmix data for graph
caps.raw = final.count[1,]
mhc.raw = c(8.8,NA,NA,2.0,1.1,NA,NA,5.6,NA,1.1,1.1,NA,NA,9.4,NA,8.9,NA,NA)
mhc.se.raw = c(3.2,NA,NA,0,0.5,NA,NA,1.0,NA,0.5,0.5,NA,NA,1.6,NA,7.1,NA,NA)
recaps = c(2,0,0,1,1,0,0,3,0,1,1,0,0,3,0,1,0,0)
naive.raw = c(6,7,1,2,1,0,1,5,1,1,1,1,0,8,0,4,3,0)
mix.raw = c(7.03,6.535,1.660,2.112,1.227,0.283,1.312,5.520,1.550,
            2.073,3.002,1.225,0.303,12.230,0.285,3.967,3.550,0.340)
mix.se.raw = c(1.67,1.168,0.799,0.336,0.461,0.532,0.573,1.857,
               0.780,1.072,1.423,0.524,0.573,2.865,0.530,1.426,1.423,0.644)

####################################

#Mouse MRR and Nmix data for graph
caps.m.raw = final.count[2,]
mhc.m.raw = c(43.3,NA,5.3,13.2,30,6.8,2,NA,NA,2.8,NA,NA,3.7,29.4,6.6,14.1,2,6.6)
mhc.m.se.raw = c(37.1,NA,3.7,2.2,16.9,2.6,0.1,NA,NA,1.7,NA,NA,1.3,17.1,1.0,3.8,0.1,1.0)
recaps.m = c(1,0,1,4,2,2,2,0,0,1,0,0,2,2,5,4,2,4)
naive.m.raw = c(12,9,3,11,12,5,2,0,0,2,0,2,3,11,6,10,2,6)
mix.m.raw = c(11.945,10.625,4.898,8.253,7.197,9.403,7.097,0.993,3.928,
              4.173,1.597,1.777,6.680,8.885,10.682,12.740,6.277,7.707)
mix.se.m.raw = c(1.982,2.334,1.774,1.099,1.123,2.051,2.139,1.111,2.099,
                 1.474,1.978,0.905,1.942,1.852,2.080,3.066,2.460,1.697)


##############

#Bind data for both species together
each = cbind(mhc.raw,mhc.se.raw,mix.raw,mix.se.raw,naive.raw,caps.raw)
wfmo = cbind(mhc.m.raw,mhc.m.se.raw,mix.m.raw,mix.se.m.raw,naive.m.raw,caps.m.raw)
both = rbind(each,wfmo)

#Count # of recounts at each site
recaps = c(5,0,1,2,1,1,0,4,2,1,1,2,2,6,1,1,2,1)
recaps.m = c(3,0,1,8,3,5,3,0,0,1,0,0,3,2,6,5,2,4)
totrecap = c(recaps,recaps.m)

#Include only sites with 3 or more recaps
select = cbind(both,totrecap)
final = as.data.frame(select[select[,7]>=3,])

#Remove 2 final sites
final = final[c(-4,-6),]

##################################################

#Figure code

#MRR estimates vs. Nmix
plot(x=final$mhc.raw,y=final$mix.raw,xlab="N estimate from MRR (AIC best model)", 
     ylab="(Relative) abundance estimate",xlim=c(0,16),ylim=c(0,40),pch=19,cex=1.5)
#MRR estimates vs. caps/100 trap nights
points(x=final$mhc.raw,y=final$caps.raw,pch=17,cex=1.5)
#MRR estimates vs. MKNA
points(x=final$mhc.raw, y=final$naive.raw, pch=21, cex=1.5)

#Draw trendlines
abline(lm(final$mix.raw~final$mhc.raw -1),lwd=3,lty=1)
abline(lm(final$caps.raw~final$mhc.raw -1),lwd=3,lty=2)
abline(lm(final$naive.raw~final$mhc.raw -1),lwd=3,lty=3)

#Add trendline formulas to figure
summary(lm(test5$mix.raw~test5$mhc.raw -1))
text(12,15.5,labels="y=1.00x, r=0.95")
summary(lm(test5$naive.raw~test5$mhc.raw -1))
text(12,6.5,labels="y=0.81x, r=0.99")
summary(lm(test5$caps.raw~test5$mhc.raw -1))
text(9,27.5,labels="y=2.41x, r=0.98")

#Add legend
legend(0,40,rev(c("N-mixture model","MKNA","Caps/100 trap nights")),
       pch=rev(c(19,21,17)),bty='n',lty=c(3,2,1))

