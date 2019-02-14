## Populaiton dynamics (discrete)

raw<-read.csv("env_log.txt",sep="\t")
data<- read.csv("out.txt", sep="\t")
dev.new()
plot(data[1:14,3],data[1:14,1], type="l", col= "orange",ylim=c(0,1000),xlab="Time",ylab="Population size")
lines(data[1:14,3],data[1:14,2], type="l", col= "skyblue")

par(new = T)
with(raw, plot(Generation[1:260], env[1:260], type="l", col="grey", axes=F, xlab=NA, ylab=NA, ylim=c(0.0,100.0)))
axis(side = 4)
mtext(side = 4, line = 3, 'Environmental conditions')

## per capita growth rate (discrete)

raw<-read.csv("env_log.txt",sep="\t")
data<- read.csv("out.txt", sep="\t")
dev.new()
plot(data[1:14,3]+0.5,data[2:15,4], type="l", col= "orange",xlim=c(0,13),xlab="Time",ylab="per capita growth rate")
lines(data[1:14,3]+0.5,data[2:15,5], type="l", col= "skyblue")

par(new = T)
with(raw, plot(Generation[1:260], env[1:260], type="l", col="grey", axes=F, xlab=NA, ylab=NA, ylim=c(0.0,100.0)))
axis(side = 4)
mtext(side = 4, line = 3, 'Environmental conditions')

## Populaiton dynamics (continuous)

raw<-read.csv("summary.txt",sep="\t")
ending<- 1000
dev.new()
plot(raw[1:ending,3],raw[1:ending,1], type="l", col= "orange", ylim=c(0,1000),xlab="Time",ylab="Population size")
lines(raw[1:ending,3],raw[1:ending,2], type="l", col= "skyblue")

par(new = T)
with(raw, plot(time[1:ending], environment[1:ending], type="l", col="grey", axes=F, xlab=NA, ylab=NA, ylim=c(0.0,100.0)))
axis(side = 4)
mtext(side = 4, line = 3, 'Environmental conditions')

## per capita growth rate (continuous)

raw<-read.csv("summary.txt",sep="\t")
ending<- 1000
dev.new()
plot(raw[1:ending,3],raw[1:ending,5], type="l", col= "orange",ylim=c(-1.5, 1.5),xlab="Time",ylab="per capita growth rate")
lines(raw[1:ending,3],raw[1:ending,6], type="l", col= "skyblue")

par(new = T)
with(raw, plot(time[1:ending], environment[1:ending], type="l", col="grey", axes=F, xlab=NA, ylab=NA, ylim=c(0.0,100.0)))
axis(side = 4)
mtext(side = 4, line = 3, 'Environmental conditions')
