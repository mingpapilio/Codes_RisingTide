## Population dynamics
raw<-read.csv("pop_dmc.txt",sep="\t")
dev.new()
plot(raw[,1],raw[,3],type="l",ylim=c(0,1050),xlab="Time",ylab="Population",col="skyblue")
points(raw[,1],raw[,4],type="l",col="orange")
legend("topleft",legend=c("rising-tide","bet-hedging"),col=c("orange","skyblue"),lty=1)

## Total population versus opportunity for selection
d<- read.csv("pop_dmc.txt",sep="\t")
span= 200
dev.new()
with(d, plot(time, Cooperative, type="l", col="black", 
             ylab="Population size",ylim=c(200,1250)))
par(new = T)
with(d, plot(time, Total,type="n", col="red2", axes=F, xlab=NA, ylab=NA, ylim=c(0.0,3.0)))
axis(side = 4)
mtext(side = 4, line = 3, 'Opportunity for selection')
lines(x=d$time-span/2,y=d$Total,col="red2")
legend("topright",
       legend=c("Total Population", "Opportunity for Selection"),
       lty=c(1,1), col=c("black", "red2"))


