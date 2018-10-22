## Plotting codes
## plopulation dynamics
d<- read.csv("out.txt",sep="\t")
dev.new()
par(mar=c(5,4,2,2))
plot(d[,3],d[,2],type="l",lwd=2,col="skyblue",xlab="Time",ylab="Population size",
     ylim=c(0,1500))
points(d[,3],d[,1],type="l",lwd=2,col="orange")

## environmental time series
raw<-read.csv("env_log.txt",sep="\t")
dev.new()
par(mar=c(5,4,4,2))
plot(raw[,1],raw[,2],type="l",lwd=2,col="light grey",xlab="Time",ylab="Environmental conditions",
     ylim=c(35,115))
abline(h=75,lty=2)
