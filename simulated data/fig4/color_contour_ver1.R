library(plotly)
library(viridis)

rr<- rep(NA, 11)
gg= bb= rr
data<-read.csv("summary.txt",sep="\t")
rep<- 100

rr[1:6]<- round(seq(86, 202, length.out= 6))
gg[1:6]<- round(seq(53, 122, length.out= 6))
bb[1:6]<- round(seq(46, 44, length.out= 6))

rr[6:11]<- round(seq(202, 250, length.out= 6))
gg[6:11]<- round(seq(122, 214, length.out= 6))
bb[6:11]<- round(seq(44, 137, length.out= 6))

trace1<- list(
  z= matrix(data$dom_1,nrow=5,ncol=5),
  type= "contour"
  ,colorscale= cbind(seq(0,1,by=1/10),c(
    rgb(rr[1]/255, gg[1]/255, bb[1]/255),
    rgb(rr[2]/255, gg[2]/255, bb[2]/255),
    rgb(rr[3]/255, gg[3]/255, bb[3]/255),
    rgb(rr[4]/255, gg[4]/255, bb[4]/255),
    rgb(rr[5]/255, gg[5]/255, bb[5]/255),
    rgb(rr[6]/255, gg[6]/255, bb[6]/255),
    rgb(rr[7]/255, gg[7]/255, bb[7]/255),
    rgb(rr[8]/255, gg[8]/255, bb[8]/255),
    rgb(rr[9]/255, gg[9]/255, bb[9]/255),
    rgb(rr[10]/255, gg[10]/255, bb[10]/255),
    rgb(rr[11]/255, gg[11]/255, bb[11]/255)
  ))
)
p <- plot_ly(
  contours = list(
    start= 0.1*rep,
    end= 0.9*rep,
    size= 0.1*rep
  ))
p<- add_trace(p, z=trace1$z, colorscale=trace1$colorscale, type=trace1$type)
p


#############
rr[1:6]<- round(seq(15, 0, length.out= 6))
gg[1:6]<- round(seq(37, 92, length.out= 6))
bb[1:6]<- round(seq(64, 175, length.out= 6))

rr[6:11]<- round(seq(0, 88, length.out= 6))
gg[6:11]<- round(seq(92, 178, length.out= 6))
bb[6:11]<- round(seq(175, 220, length.out= 6))

trace1<- list(
  # z= matrix(data$dom_1,nrow=5,ncol=5),
  z= matrix(data$dom_2, nrow=5, ncol=5),
  type= "contour"
  ,colorscale= cbind(seq(0,1,by=1/10),c(
    rgb(rr[1]/255, gg[1]/255, bb[1]/255),
    rgb(rr[2]/255, gg[2]/255, bb[2]/255),
    rgb(rr[3]/255, gg[3]/255, bb[3]/255),
    rgb(rr[4]/255, gg[4]/255, bb[4]/255),
    rgb(rr[5]/255, gg[5]/255, bb[5]/255),
    rgb(rr[6]/255, gg[6]/255, bb[6]/255),
    rgb(rr[7]/255, gg[7]/255, bb[7]/255),
    rgb(rr[8]/255, gg[8]/255, bb[8]/255),
    rgb(rr[9]/255, gg[9]/255, bb[9]/255),
    rgb(rr[10]/255, gg[10]/255, bb[10]/255),
    rgb(rr[11]/255, gg[11]/255, bb[11]/255)
  ))
)
p <- plot_ly(
  contours = list(
    start= 0.1*rep,
    end= 0.9*rep,
    size= 0.1*rep
  ))
p<- add_trace(p, z=trace1$z, colorscale=trace1$colorscale, type=trace1$type)
p

###########
rr[1:6]<- round(seq(54, 27, length.out= 6))
gg[1:6]<- round(seq(86, 129, length.out= 6))
bb[1:6]<- round(seq(60, 62, length.out= 6))

rr[6:11]<- round(seq(27, 168, length.out= 6))
gg[6:11]<- round(seq(129, 216, length.out= 6))
bb[6:11]<- round(seq(62, 185, length.out= 6))

trace1<- list(
  # z= matrix(data$dom_1,nrow=5,ncol=5),
  # z= matrix(data$dom_2, nrow=5, ncol=5),
  z= matrix(data$coexist,nrow=5, ncol=5),
  type= "contour"
  ,colorscale= cbind(seq(0,1,by=1/10),c(
    rgb(rr[1]/255, gg[1]/255, bb[1]/255),
    rgb(rr[2]/255, gg[2]/255, bb[2]/255),
    rgb(rr[3]/255, gg[3]/255, bb[3]/255),
    rgb(rr[4]/255, gg[4]/255, bb[4]/255),
    rgb(rr[5]/255, gg[5]/255, bb[5]/255),
    rgb(rr[6]/255, gg[6]/255, bb[6]/255),
    rgb(rr[7]/255, gg[7]/255, bb[7]/255),
    rgb(rr[8]/255, gg[8]/255, bb[8]/255),
    rgb(rr[9]/255, gg[9]/255, bb[9]/255),
    rgb(rr[10]/255, gg[10]/255, bb[10]/255),
    rgb(rr[11]/255, gg[11]/255, bb[11]/255)
  ))
)
p <- plot_ly(
  contours = list(
    start= 0.1*rep,
    end= 0.9*rep,
    size= 0.1*rep
  ))
p<- add_trace(p, z=trace1$z, colorscale=trace1$colorscale, type=trace1$type)
p