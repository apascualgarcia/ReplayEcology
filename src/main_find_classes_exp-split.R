 ##########################################
#  main_find_classes_exp-split.R
##########################################
# Here I extend the pipeline suggested to analyse enterotypes 
# (http://enterotype.embl.de/enterotypes.html#sim) to analyse
# the existence of classes. I expand the more basic script
# "main_find_classes.R" which is intended to generate an all-against-all
# distance matrix across all experiments and to perform a clustering
# analysis, to repeat this analysis for each specific experiment
# and replicate separately.
# 
# author = apascualgarcia.github.io
# date = May 11th, 2022. ETH-Zürich
#
library(fpc) # cluster.stats, pamk

rm(list=ls())
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)
source("pam.clustering.R")
source("dist.JSD.R")
source("otu_table_to_classes.R")
source("plot_calinski.R")
source("clustering_to_partition.R")
source("indexG1.check.R")
source("print_pamk_to_partition.R")
source("add_metadata.R")
source("dist_to_pcoa.R")


#source("")

# Read files ------------
# --- Read ASVs table
setwd("../6_finalfiles")
fileOTU="seqtable_readyforanalysis.csv"
ASV.table=read.table(fileOTU,sep="\t")
dim(ASV.table)
ASV.table=as.matrix(ASV.table)
colnames(ASV.table)[1:5]
rownames(ASV.table)[1:5]
ASV.list=list()

# ..... read metadata
setwd("../4_dada2/")
fileMD="metadata_Time0D-7D-4M_May2022.csv"
sample_md<-read.table(fileMD,sep="\t",header=TRUE)

# ..... read distance matrix (This distance is computed in main_find_classes.R)
setwd("../7.1_classes")
fileDist="Dist_JSD_Time0D-7D-4M.RDS"
data.dist=readRDS(fileDist) 
data.dist.list=list()

# Subset ------------
# --- First remove samples that did not pass the quality control from the md
matched=match(sample_md$sampleid,row.names(ASV.table)) #attr(data.dist,"Labels"))
sample_md=sample_md[!is.na(matched),]

# --- Extract ASV tables at different time points and replicates
#level=levels(sample_md$Experiment)[1]
experiments=levels(sample_md$Experiment)
labels=c()
for(level in experiments){
  idx=which(sample_md$Experiment == level)
  matched=match(sample_md$sampleid[idx],row.names(ASV.table))
  matched=matched[!is.na(matched)]
  data.dist.tmp=as.dist(as.matrix(data.dist)[matched, matched])
  ASV.tmp=ASV.table[matched,]
  #dim(ASV.tmp)
  ASV.tmp=ASV.tmp[,colSums(ASV.tmp) > 0]
  el.id=paste0("Time",level)
  labels=c(labels,el.id)
  data.dist.list[[el.id]]=data.dist.tmp
  ASV.list[[el.id]]=ASV.tmp
}
lapply(ASV.list,FUN = dim)
lapply(data.dist.list,FUN = length)

# Classification ----------
# --- Compute classification at different thresholds
#ASV.test=ASV.table[1:10,]
#
nclusters.list=list()
Kmax.list=list()
method="clusterSim" # One of "pamk" or  "clusterSim"
#label="Time0D"
Nclusters=90
for(label in labels){
  ASV.tmp=ASV.list[[label]]
  data.dist.tmp=data.dist.list[[label]]
  if(method == "pamk"){
    out_classes=pamk(data.dist.tmp,krange=1:Nclusters,
                     criterion = "ch")
    nclusters.list[[label]]=out_classes$crit
    Kmax.list[[label]]=out_classes$nc
    print_pamk_to_partition(out_classes,outputlabel = label)
  }else if(method == "clusterSim"){
    out_classes=otu_table_to_classes(t(ASV.tmp),Nclus = Nclusters,
                                     distance = data.dist.tmp,
                                     euclidean = TRUE) # centroids = FALSE, medoids = TRUE
    nclusters.list[[label]]=out_classes$nclusters
    Kmax.list[[label]]=which.max(nclusters.list[[label]])
  }
}
#Kmax.list[[1]]=2

# .... Plot results
setwd("../7.1_classes")
for(label in labels){
  nclusters.tmp=nclusters.list[[label]]
  if(method == "pamk"){
    filename="Plot_PAMk-CH_JSD"
  }else if(method == "clusterSim"){
    filename="Plot_PAM-IndexG1_orig-medoidsV2"
  }
  plot_calinski(nclusters.tmp,
                filename = filename,
                outputlabel = label) # check the maxima
}

# .... Extract partition
partition.list=list()
rownames(sample_md)=sample_md$sampleid
sample_md_exp=sample_md
if(method == "clusterSim"){ # if pamk it was already printed
  for(label in labels){
    Kmax.tmp=Kmax.list[[label]]
    data.dist.tmp=data.dist.list[[label]]
    partition=clustering_to_partition(data.dist.tmp,Kmax = Kmax.tmp,
                                      outputlabel = label,
                                      print.file = TRUE)
    label_part=paste("Part",label,Kmax.tmp,sep="_")
    sample_md_exp=add_metadata(sample_md_exp,partition,label_part)
    if((label == "Time0D")&(Kmax.tmp != 6)){ # we also want NatComm results if we do not find them
      Kmax.tmp=6
      partition=clustering_to_partition(data.dist.tmp,Kmax = Kmax.tmp,
                                        outputlabel = label,
                                        print.file = TRUE)
      label_part=paste("Part",label,Kmax.tmp,sep="_")
      sample_md_exp=add_metadata(sample_md_exp,partition,label_part)
    }
  }
}

# --- Complete "replicate" for samples != 4M and create a single 
#     "partition" and "replicate.partition" factors
sample_md$partition=NA
sample_md_exp$replicate.partition=NA
for(level in experiments){
  cols=grep(colnames(sample_md_exp),pattern = level)
  if(level == "0D"){
    cols = cols[2] # we have two partitions here
    replicate = 0
  }else if(level == "7D_rep1"){
    replicate = 1
  }else if(level == "7D_rep2"){
    replicate = 2
  }else if(level == "7D_rep3"){
    replicate = 3
  }else if(level == "7D_rep4"){
    replicate = 4
  }else if(level == "4M"){
    replicate = 0
  } 
  id=which(sample_md_exp$Experiment == level)
  sample_md_exp$replicate[id]=replicate # as.character(sample_md_exp$Experiment[id])
  sample_md_exp$partition[id]=sample_md_exp[id,cols]
  sample_md_exp$replicate.partition[id]=paste(sample_md_exp$replicate[id],
                                          sample_md_exp[id,cols],sep=".")
}

# --- Write new metadata table
fileOut=unlist(strsplit(fileMD,".",fixed=T))[1]
fileOut=paste0(fileOut,"_wJSDpart-split.csv")
write.table(sample_md_exp,file = fileOut,sep="\t",quote=F,row.names = F)

# Generate PCoA plots ---------

for(label in labels){
  id.part=grep(label,colnames(sample_md_exp))
  part.names=colnames(sample_md_exp)[id.part] # identify all partitions with that name
  data.dist.tmp=data.dist.list[[label]]
  for(part.name.tmp in part.names){
    dist_to_pcoa(data.dist.tmp,sample_md_exp,
                 factor1=part.name.tmp,factor2 = "Location",
                 vec.coor=c(1,2,3))
    dist_to_pcoa(data.dist.tmp,sample_md_exp,
                 factor1=part.name.tmp,factor2 = "Location",
                 vec.coor=c(2,3,1))
    dist_to_pcoa(data.dist.tmp,sample_md_exp,
                 factor1=part.name.tmp,factor2 = "Location",
                 vec.coor=c(1,3,2))
  }
}
