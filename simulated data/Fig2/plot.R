## Plotting codes
raw<-read.csv("out.txt",sep="\t")
dev.new()
par(mar=c(5,4,2,2))
plot(raw[,3],raw[,2],type="l",lwd=2,col="skyblue",xlab="Time",ylab="Population size",
     ylim=c(0,1500))
points(raw[,3],raw[,1],type="l",lwd=2,col="orange")

dev.new()
par(mar=c(5,4,4,2))
plot(raw[,3],raw[,4],type="l",lwd=2,col="light grey",xlab="Time",ylab="Encironmental conditions",
     ylim=c(35,115))
abline(h=75,lty=2)
