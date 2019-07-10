
result<- matrix(NA,3,2)
rownames(result)<- c("0.3","0.5","0.7")
colnames(result)<- c("rit","bet")

rep<- 1000
aa<- 0
bb<- 0
ibr_alt<- read.csv("summary_0.3.txt",sep="\t",header=T)
for(i in 1:rep){
  if(ibr_alt[i,2]==1){ aa= aa+1}  ## bet-hedging
  if(ibr_alt[i,2]==-1){ bb= bb-1} ## rising tide
}
result[1,2]=aa/rep
result[1,1]=-bb/rep
##
aa<- 0
bb<- 0
ibr_alt<- read.csv("summary_0.5.txt",sep="\t",header=T)
for(i in 1:rep){
  if(ibr_alt[i,2]==1){ aa= aa+1}
  if(ibr_alt[i,2]==-1){ bb= bb-1}
}
result[2,2]=aa/rep
result[2,1]=-bb/rep
##
aa<- 0
bb<- 0
ibr_alt<- read.csv("summary_0.7.txt",sep="\t",header=T)
for(i in 1:rep){
  if(ibr_alt[i,2]==1){ aa= aa+1}
  if(ibr_alt[i,2]==-1){ bb= bb-1}
}
result[3,2]=aa/rep
result[3,1]=-bb/rep
##
result