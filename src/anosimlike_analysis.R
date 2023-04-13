##################################################
# anosimlike_analysis.R
##################################################
# In this script I compare the betadiversity variability
# between samples sharing a commong factor, such as 
# belonging to the same class and replicate, same
# class independently of the replicate, or same parent
# community, among other possibilities.
#
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list = ls())
library(usedist)
library(vegan)

# START EDITING -----------
# --- Set the files needed
file.meta="metadata_Time0D-7D-4M_May2022_wJSDpart.csv" # metadata with partitions
fileDist="Dist_JSD_Time0D-7D-4M.RDS" # beta diversity distance

# --- Select category to classify the distances and prepare structures
category="parent" # "partition" #"time" #"experiment.location" #"replicate.partition"# "Experiment"

# --- Select a method to analyse the variance between groups
setMethod="anosim" # One of "anosim", "adonis2" or "mrpp"
Nperm0 = 999 # number of permutations to run the test

# --- Experiments to purge (optionally, 4M and class 6 will always be purged)
purge = "0D" # one of "0D", "7D" or "none"

# --- Replicate to keep (if you want to analyse a single replicate)
keep = "none" # one of "Rep1", "Rep2", "Rep3", "Rep4" or "none"
# STOP EDITING -----------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

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
sample_metadata$time = as.character(sample_metadata$Experiment)
id.7D = grep("7D",sample_metadata$time)
sample_metadata$time[id.7D] = "7D"


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

# --- Optional experiments to purge

if(purge != "none"){
  id.keep = as.character(sample_metadata$sampleid[which(sample_metadata$Experiment != purge)])
  matched=match(id.keep, sample_metadata$sampleid)
  sample_metadata=sample_metadata[matched,]
}


# --- Extract a single replicate
if(keep != "none"){
  id.keep = as.character(sample_metadata$sampleid[which(sample_metadata$replicate == keep)])
  matched=match(id.keep, sample_metadata$sampleid)
  sample_metadata=sample_metadata[matched,]
}


# --- Extract distance matrix with desired samples
Nsamples=dim(sample_metadata)[1]
dist.tmp = matrix(data = NA, nrow = Nsamples, ncol = Nsamples)
for(i in 1:Nsamples){
  for(j in i:Nsamples){
    sample1=as.character(sample_metadata$sampleid[i])
    sample2=as.character(sample_metadata$sampleid[j])
    par1=as.character(sample_metadata[i,category])
    par2=as.character(sample_metadata[j,category])
    if(i == j){
      dist.tmp[i, j] = 0
    }else{
      dist.tmp[i, j] = dist_get(dist,sample1,sample2)
      dist.tmp[j, i] = dist_get(dist,sample1,sample2)
      
    }
  }
}
dist.tmp = as.dist(dist.tmp)

# Start comparison -----

partition=sample_metadata[[category]] # this is a vector
if(setMethod=="anosim"){
  AA=anosim(dist.tmp, as.factor(partition), permutations = Nperm0)
  (anosimSumm=AA$statistic)
  (anosimSig=AA$signif)
}else if(setMethod=="mrpp"){
  AA=mrpp(dist.tmp, as.factor(partition), permutations = Nperm0)
  (anosimSumm=AA$delta)
  (anosimSig=AA$Pvalue)
}else{
  tmp.frame=as.factor(partition)
  colnames(tmp.frame)="X"
  AA=adonis2(dist.tmp ~ X , tmp.frame, permutations = Nperm0)
  (anosimSumm=AA$F[1])
  (anosimSig=AA$`Pr(>F)`[1])
}
