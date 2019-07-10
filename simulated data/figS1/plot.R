## Environmental time series
raw<-read.csv("summary.txt",sep="\t")

dev.new()
plot(raw[2:4000,3],raw[2:4000,4], type="l", col="black",
     xlab="Time",
     ylab="Environmental conditions"
     , ylim=c(0,100)
)
abline(h=50,lty=2)

## FFT analysis
data<-read.csv("summary.txt",sep="\t",skip=2)
aa.spec<- spectrum(data[,4],log="no",span=5,plot=F)
dev.new()
plot(aa.spec$spec~aa.spec$freq,xlab="frequency",ylab="spectral density",log="xy",type="l",ylim=c(5e-5,2e+5))

