
fdata<- read.csv("summary.txt", sep="\t")

dev.new()
plot(fdata$Environment,fdata$abf_r,
     type="l",
     ## ylim=c(-0.5, 9),
     ## ylim=c(-0.05, 0.2),
     ylim=c(-0.5, 5),
     xlab="environmental condition",
     ylab="per capita growth rate",
     col="orange")
points(fdata$Environment,fdata$abf_b, type="l",col="skyblue")
polygon(fdata$Environment,fdata$abf_b,col=rgb(0.5,0.5,0.5,0.5),border=NA)

c(mean(fdata$abf_r), sd(fdata$abf_r))
c(mean(fdata$abf_b), sd(fdata$abf_b))