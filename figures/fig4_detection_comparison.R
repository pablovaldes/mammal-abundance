##########################################################
## Figure 4: Comparing estimated p by harvest treatment ##
##########################################################

#Read in data p values are obtained from
load('output/output_hee_abundance_Nmix.Rdata')

#Set up figure structure
par(mfrow = c(1,2))

#Pre/Post Harvest p estimates - Chipmunk
#a = 0.4ha, b = 2ha, c = 4ha (clearcut), d = shelter, e = control
a = c(0.486,0.399)
b = c(0.443,0.384)
c = c(0.416,0.366)
d = c(0.400,0.500)
e = c(0.423,0.396)

data = cbind(a,b,c,d,e)

#Plot
barplot(data, beside=TRUE, col=gray(seq(0.1,0.9,length=2)),ylim=c(0,0.6),
        ylab="Mean detection probability",xlab=c('Treatment'),names=c('0.4 ha','2 ha','4 ha','Shelt','Control'),
        main="Chipmunk")
legend(1,0.6,c("Before harvest", "After harvest"),
       fill=gray(seq(0.1,0.9,length=2)))

####################################################################

#Pre/Post Harvest p estimates - Mouse
a = c(0.539,0.442)
b = c(0.401,.337)
c = c(0.370,0.508)
d = c(0.338,0.500)
e = c(0.383,0.251)

data = cbind(a,b,c,d,e)

#Plot
barplot(data, beside=TRUE, col=gray(seq(0.1,0.9,length=2)),ylim=c(0,0.6),
        ylab="Mean detection probability",xlab=c('Treatment'),names=c('0.4 ha','2 ha','4 ha','Shelt','Control'),
main="Mouse")
