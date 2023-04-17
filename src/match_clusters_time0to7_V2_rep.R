################################################# 
#        match_clusters_time0to7_V2.R
################################################# 
#  In this version the match has already been done
#  with the Perl script  TraceSamples_fromPartTime0toPartTime7.pl,
#  and here we create a matrix plot with the number of communities in partitions
#  at time 0 going to partitions in time 7
# ##############################
#  A.P-G., ETH, 
#  Zürich, 31 July 2018
# ##############################

rm(list=ls())

library(reshape2)
library(ggplot2)
#library(vcd)


# START EDITING -----------
label="SJD" # "SJD" if all samples are considered,  "SJD_Time0restrictPar" for the restricted set 
#fileIn="Matrix_TracingSamplesTime0toTime7_SJD.out"
fileIn=paste("Matrix_TracingSamplesTime0toTime7_",label,".out",sep="")

# STOP editing here

# Load data ----------------
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)
setwd("../7.2_match")

if(label=="SJD"){ # 5 if all samples are considered, 3 for the restricted set 
  Npar=5
}else{
  Npar=3  
}
match.matrix=as.matrix(read.table(file=fileIn,sep="\t"))
match.matrix=match.matrix[1:Npar,]
norm=match.matrix[,1]+match.matrix[,2] # Number of communities for each partition at time0
norm2=colSums(match.matrix)
match.matrix.norm=round(sweep(match.matrix,1,norm,"/"),digits=2)
#match.matrix.norm.trans=t(match.matrix.norm)
dimnames(match.matrix)=list(Time0=rownames(match.matrix),Time7=colnames(match.matrix))
match.df=melt(match.matrix)
colnames(match.df)[3]="Num. Communities"

dimnames(match.matrix.norm)=list(Time0=rownames(match.matrix.norm),Time7=colnames(match.matrix.norm))
match.df.norm=melt(match.matrix.norm)
colnames(match.df.norm)[3]="% Communities"

fileOut=paste("TraceSamples_Time0toTime7_",label,".pdf",sep="")
pdf(file=fileOut,width=10,height=8)
ggplot(match.df,aes(Time7, Time0)) +
  scale_fill_gradient(low="white",high="red1")+
  geom_tile(aes(fill = `Num. Communities`), colour = "black")+
  geom_text(aes(label=`Num. Communities`),hjust=0, vjust=0)
dev.off()

fileOut=paste("TraceSamples_Time0toTime7_Percent_",label,".pdf",sep="")
pdf(file=fileOut,width=11,height=4)
ggplot(match.df.norm,aes(Time7, Time0,size=12)) +
  scale_fill_gradient(low="white",high="red1")+
  geom_tile(aes(fill = `% Communities`), size=1,colour = "black")+
  geom_text(aes(label=`% Communities`),size=8,hjust=0, vjust=0,nudge_x = -0.3,nudge_y = -0.2)+
  theme(legend.title = element_text(size=28),legend.text = element_text(size = 23)) # change size legend title
dev.off()

