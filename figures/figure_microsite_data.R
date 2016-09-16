######################################
## Figure: Microsite data summaries ##
######################################

#Summaries of microsite information, used in a poster

#Read in microsite.parse() function
source('function_microsite_parse.R')

#Arguments for microsite.parse()
ngrids = 32
nyears = 5
ntraplocs = 27

missing.grids2007 = c(1113,1125,1312,1317,1401,1404,1803,1820,1821,1830)

missing.grids2010 = c(1001,1002,1105,1113,1121,1125,1703,1705,1716,1820,1821,1830)

traps.list = c('A1','A2','A3','A4','A5','A6',
             'B1','B2','B3','B4','B5','B6',
             'C1','C2','C3','C4','C5','C6')

var.list = c('herb','wood','cwd','litter')

#Parse microsite data and get mean values for grids
micro = array(data=NA, dim=c(ngrids,length(var.list),nyears))
micro[,,1] = microsite.parse(file='data/2007_microsite_scrub.csv',var.list=var.list,
                              traps=traps.list,missing.grids=missing.grids2007)$micro.mean
micro[,,2] = microsite.parse(file='data/2008_microsite_scrub.csv',var.list=var.list,traps=traps.list)$micro.mean
micro[,,3] = microsite.parse(file='data/2009_microsite_scrub.csv',var.list=var.list,traps=traps.list)$micro.mean
micro[,,4] = microsite.parse(file='data/2010_microsite_scrub.csv',var.list=var.list,
                              traps=traps.list,missing.grids=missing.grids2010)$micro.mean
micro[,,5] = microsite.parse(file='data/2011_microsite_scrub.csv',var.list=var.list,traps=traps.list)$micro.mean

#Set grid treatment indexes
control = c(1,2,7,8,13,14,15,16)
clearcut = c(9,11,17,20,31,32)
patch = c(3,4,5,6,21,22,23,24,25,26,27,28)
shelter = c(10,12,18,19,29,30)

#Create treatment x year microsite means
clearmean = matrix(NA,nrow=4,ncol=5)
patchmean = matrix(NA,nrow=4,ncol=5)
sheltermean = matrix(NA,nrow=4,ncol=5)
conmean = matrix(NA,nrow=4,ncol=5)

for (i in 1:4){
  conmean[i,1:5] = colMeans(micro[control,i,1:5],na.rm=TRUE)
  clearmean[i,1:5] = colMeans(micro[clearcut,i,1:5],na.rm=TRUE)
  patchmean[i,1:5] = colMeans(micro[patch,i,1:5],na.rm=TRUE)
  sheltermean[i,1:5] = colMeans(micro[shelter,i,1:5],na.rm=TRUE)
}

############################################################################

#Begin figure code

#Set up charg structure 5 x 4
par(mfrow=c(4,5), oma = c(4,1,1,1), mar=c(1,4.5,3,1.6))

#Clearcut Figures

#Acorn production data hardcoded in
barplot(c(2.68,7.53,4.37,0,24.88),ylab="Clearcut", 
        main="Acorns / sq m", col=rgb(red=244,green=125,blue=66, maxColorValue=255), ylim=c(0,25))
abline(v=2.52, lty=2)

barplot(clearmean[1,1:5],
        main="% Herb Cover", col=rgb(red=244,green=125,blue=66, maxColorValue=255), ylim=c(0,70), cex.lab=1.5)
abline(v=2.52, lty=2)

barplot(clearmean[2,1:5], main="% Woody Cover", col=rgb(red=244,green=125,blue=66, maxColorValue=255),ylim=c(0,70))
abline(v=2.52, lty=2)

barplot(clearmean[3,1:5], main="CWD (m)", col=rgb(red=244,green=125,blue=66, maxColorValue=255),ylim=c(0,3))
abline(v=2.5, lty=2)

barplot(clearmean[4,1:5], main="Litter Depth (cm)", col=rgb(red=244,green=125,blue=66, maxColorValue=255),ylim=c(0,4))
abline(v=2.52, lty=2)

#Patch Cut Figures

barplot(c(14.42,20.62,0,0.31,20.78),ylab="Patch Cut", 
        col=rgb(red=147,green=146,blue=20, maxColorValue=255), ylim=c(0,25))
abline(v=2.52, lty=2)

barplot(patchmean[1,1:5], col=rgb(red=147,green=146,blue=20, maxColorValue=255), ylim=c(0,70), 
        cex.lab=1.5)
abline(v=2.52, lty=2)

barplot(patchmean[2,1:5], col=rgb(red=147,green=146,blue=20, maxColorValue=255),ylim=c(0,70))
abline(v=2.52, lty=2)

barplot(patchmean[3,1:5], col=rgb(red=147,green=146,blue=20, maxColorValue=255),ylim=c(0,3))
abline(v=2.5, lty=2)

barplot(patchmean[4,1:5], col=rgb(red=147,green=146,blue=20, maxColorValue=255),ylim=c(0,4))
abline(v=2.52, lty=2)

#Shelterwood

barplot(c(5.46,17.36,3.04,2.06,15.54),ylab="Shelterwood", 
        col=rgb(red=141,green=213,blue=18, maxColorValue=255), ylim=c(0,25))
abline(v=2.52, lty=2)

barplot(sheltermean[1,1:5], col=rgb(red=141,green=213,blue=18, maxColorValue=255), ylim=c(0,70), 
        cex.lab=1.5)
abline(v=2.52, lty=2)

barplot(sheltermean[2,1:5], col=rgb(red=141,green=213,blue=18, maxColorValue=255),ylim=c(0,70))
abline(v=2.52, lty=2)

barplot(sheltermean[3,1:5], col=rgb(red=141,green=213,blue=18, maxColorValue=255),ylim=c(0,3))
abline(v=2.5, lty=2)

barplot(sheltermean[4,1:5], col=rgb(red=141,green=213,blue=18, maxColorValue=255),ylim=c(0,4))
abline(v=2.52, lty=2)

#Control

barplot(c(8.98,16.54,5.75,0.44,19.30),ylab="Control", 
        col=rgb(red=75,green=142,blue=26, maxColorValue=255), ylim=c(0,25))
abline(v=2.52, lty=2)

barplot(conmean[1,1:5], col=rgb(red=75,green=142,blue=26, maxColorValue=255), ylim=c(0,70), 
        xlab='Year',cex.lab=1.5)
abline(v=2.52, lty=2)

barplot(conmean[2,1:5], col=rgb(red=75,green=142,blue=26, maxColorValue=255),ylim=c(0,70),xlab='Year')
abline(v=2.52, lty=2)

barplot(conmean[3,1:5], col=rgb(red=75,green=142,blue=26, maxColorValue=255),ylim=c(0,3),xlab='Year')
abline(v=2.5, lty=2)

barplot(conmean[4,1:5], col=rgb(red=75,green=142,blue=26, maxColorValue=255),ylim=c(0,4),xlab='Year')
abline(v=2.52, lty=2)

