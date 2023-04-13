##################################################
# compare_clusters.R
##################################################
# In this script I compare the mean beta diversity distance between
# samples belonging to different clusters in different experiments
# and I test the significance.
#
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list=ls())

#library(RDPutils) # Notused
library(ggplot2)
library(usedist)
library(gplots) # needed for heatmap.2

# START EDITING -----------
# --- Set the files needed
file.meta="metadata_Time0D-7D-4M_May2022_wJSDpart.csv" # metadata with partitions
fileDist="Dist_JSD_Time0D-7D-4M.RDS" # beta diversity distance
# --- Select category to classify the distances and prepare structures
category="replicate.partition" #"experiment.location" #"replicate.partition"
# STOP EDITING -----------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)
source("output_matrices.R")
source("functionArgsList.R")
source("heatmap.2.mod.R")

#  --- Load metadata
setwd("../7.1_classes")
sample_metadata <-read.table(file.meta,sep="\t",header=TRUE)
sample_metadata$replicate=as.factor(paste("Rep",sample_metadata$replicate,sep=""))
sample_metadata$partition=as.factor(paste("Par",sample_metadata$partition,sep=""))
sample_metadata$replicate.partition=as.factor(paste(sample_metadata$replicate,"."
                                                    ,sample_metadata$partition,sep=""))
sample_metadata$experiment.replicate.partition=as.factor(paste(sample_metadata$Experiment,"."
                                                    ,sample_metadata$replicate.partition,sep="."))
sample_metadata$experiment.location=as.factor(paste(sample_metadata$Experiment,"."
                                                               ,sample_metadata$Location,sep="."))

# --- Load beta diversity distance
dist=readRDS(fileDist) # finally used the one computed internally in ordinate, see below

# --- Purge undesired experiments
exclude="4M"
id.keep=as.character(sample_metadata$sampleid[which(sample_metadata$Experiment != exclude)])
#sample_metadata=sample_metadata[sample_metadata$Experiment != exclude,]
exclude="Rep0.Par6"
id.keep=c(id.keep,as.character(sample_metadata$sampleid[which(sample_metadata$replicate.partition != exclude)]))
id.keep=id.keep[duplicated(id.keep)] # we retain only those that fulfilled both conditions
matched=match(id.keep, sample_metadata$sampleid)
sample_metadata=sample_metadata[matched,]

# Start comparison -----

partitions=as.character(unique(sample_metadata[,category]))
cum.dist=matrix(data=0,nrow=length(partitions),ncol=length(partitions))
cum.pairs=matrix(data=0,nrow=length(partitions),ncol=length(partitions))
colnames(cum.dist)=partitions; rownames(cum.dist)=partitions
colnames(cum.pairs)=partitions; rownames(cum.pairs)=partitions

# ---Extract distances and accumulate values
Nsamples=dim(sample_metadata)[1]
for(i in 1:Nsamples){
  for(j in i:Nsamples){
    sample1=as.character(sample_metadata$sampleid[i])
    sample2=as.character(sample_metadata$sampleid[j])
    par1=as.character(sample_metadata[i,category])
    par2=as.character(sample_metadata[j,category])
    if(i == j){
      dist.tmp=0
    }else{
      dist.tmp=dist_get(dist,sample1,sample2)
    }
    cum.dist[par1,par2]=cum.dist[par1,par2]+dist.tmp
    cum.dist[par2,par1]=cum.dist[par2,par1]+dist.tmp
    cum.pairs[par1,par2]=cum.pairs[par1,par2]+1
    cum.pairs[par2,par1]=cum.pairs[par2,par1]+1
    
  }
}
mean.dist=cum.dist/cum.pairs
#mean.dist=mean.dist[1:13,1:13] # exclude class 6

# Plot --------
# **Warning** the following string was created manually, and it depends
#             on the order of aggregation. Check rownames(mean.dist)
# --- 4M and Class 6 excluded
#"Rep1.Par1" "Rep2.Par1" "Rep3.Par1" "Rep4.Par1" 
# "Rep2.Par2" "Rep3.Par2" "Rep4.Par2" "Rep1.Par2"
#"Rep0.Par1" "Rep0.Par2" "Rep0.Par3" "Rep0.Par4" "Rep0.Par5"
newNames=c("Final.Rep1.Class1","Final.Rep2.Class1","Final.Rep3.Class1","Final.Rep4.Class1", 
           "Final.Rep2.Class2","Final.Rep3.Class2", "Final.Rep4.Class2","Final.Rep1.Class2",
          "Starting.Class1", "Starting.Class2","Starting.Class3" ,"Starting.Class4", "Starting.Class5" )
rownames(mean.dist)=newNames
colnames(mean.dist)=newNames

fileDist=paste0("meanDist_",category)
output_matrices(mean.dist,name.mat =fileDist,plot.heatmap = TRUE,
                par.pdf = list(width=30,height=30), # for categories with many elements
                par.heatmap = list(margins=c(40,40),
                                   cexRow=5,cexCol=5,
                                   keysize=1,
                                   key.par = list(cex.main=0.2,cex.axis=4, # for categories with many elements
                                                  mar=c(10,5,5,5),
                                                  mgp=c(6,3,0),
                                                  cex.lab=6)))
