## FFT analysis
## This uses the same data as Fig. 2
data<-read.csv("out.txt",sep="\t",skip=2)
aa.spec<- spectrum(data[,4],log="no",span=5,plot=F)
dev.new()
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l",ylim=c(5e-5,2e+5))
