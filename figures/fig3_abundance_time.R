#######################################################
## Figure 3: Estimated abundance by species and year ##
#######################################################

#Read in data
load('output/output_hee_abundance_Nmix.Rdata')

#Grab output from JAGS
control = abundance.out$sims.list$yearly.control
shelter = abundance.out$sims.list$yearly.shelter
patch04 = abundance.out$sims.list$yearly.patch04
patch2 = abundance.out$sims.list$yearly.patch2
clear = abundance.out$sims.list$yearly.clear

#Bundle into data array
data = array(data=NA,dim=c(abundance.out$mcmc.info$n.samples,4,5,5))
data[,,,1] = control
data[,,,2] = shelter
data[,,,3] = patch04
data[,,,4] = patch2
data[,,4,4] = NA
data[,,,5] = clear

#Graph function
abun.graph = function(species,ymax,title,ylabel,xlabel,leg){

#Check if ylabel should be included
if(ylabel==TRUE){
  plot(1,type='n', xaxt="n",xlab="",ylab="Relative abundance", main=title, 
  xlim=c(0.5,5.5), ylim=c(0,ymax))}
if(ylabel==FALSE){
  plot(1,type='n', xaxt="n",xlab="",ylab="", main=title, 
  xlim=c(0.5,5.5), ylim=c(0,ymax))}

#Check if xlabel should be included
if(xlabel==TRUE){
axis(side=1,at=c(1,2,3,4,5), labels=c('2007', '2008', '2009','2010','2011'),
 tck=0)}

#Set colors/species index
colors = gray(seq(0.4,0.1,length=5))
s=species

#Generate coordinates for lines
#Plot error bars
xcoord = matrix(nrow=5,ncol=5)
ycoord = matrix(nrow=5,ncol=5)
structure = c(1,2,3,4,5)
offset = c(-0.2,-0.1,0,0.1,0.2)
points = c(17,18,21,22,24)
for (i in 1:5){
  for (j in 1:5){
    ycoord[i,j] = mean(data[,s,i,j])
    xcoord[i,j] = i + offset[j]
    lims <- quantile(data[,s,i,j],c(0.1587,0.8413),na.rm=T)
    points(xcoord[i,j],ycoord[i,j],pch=points[j],col=colors[j],cex=1.5)
    segments(x0=xcoord[i,j],y0=lims[1],x1=xcoord[i,j],y1=lims[2], col=colors)
  }
}

#Draw pre/post harvest separation line
abline(v=2.5)

#Draw lines connecting years (accounting for missing 2 ha data in 2010)
for (k in 1:3){
  lines(x=xcoord[,k],y=ycoord[,k],type="l",pch=points[k],cex=1.5,col=colors[k])
}
for (k in 4:4){
  lines(x=xcoord[1:3,k],y=ycoord[1:3,k],type="l",pch=points[k],cex=1.5,col=colors[k])
}
for (k in 5:5){
  lines(x=xcoord[,k],y=ycoord[,k],type="l",pch=points[k],cex=1.5,col=colors[k])
}

#Plot legend
if(leg==TRUE){
legend('topleft',rev(c("Control","Shelter","4 ha cut","2 ha cut", "0.4 ha cut")),cex=1,
pch=c(21,22,24,18,17),lwd=1,col=rev(colors),bg="white",bty="n")}
}


#Run abundance graph function for each species
par(fig=c(0,0.53,0.43,1),new=FALSE)
abun.graph(1,14,"Eastern chipmunk",ylabel=TRUE,xlabel=FALSE,leg=TRUE)
par(fig=c(0.47,1,0.43,1),new=TRUE)
abun.graph(2,20,"White-footed mouse",ylabel=FALSE,xlabel=FALSE,leg=FALSE)
par(fig=c(0,0.53,0,0.57),new=TRUE)
abun.graph(3,5,"Short-tailed shrew",ylabel=TRUE,xlabel=TRUE,leg=FALSE)
par(fig=c(0.47,1,0,0.57),new=TRUE)
abun.graph(4,2.5,"Pine vole",ylabel=FALSE,xlabel=TRUE,leg=FALSE)
