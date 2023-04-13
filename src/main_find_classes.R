##########################################
#  main_find_classes.R
##########################################
# Here I extend the pipeline suggested to analyse enterotypes 
# (http://enterotype.embl.de/enterotypes.html#sim) to analyse
# the existence of classes across experiments.
# 
# author = apascualgarcia.github.io
# date = May 11th, 2022. ETH-Zürich
#

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
source("dist_to_pcoa.R")

####### START EDITING
# Set parameters -------
nreads = 10000 # minimum number of reads to consider a sample
# ... First time should be set to true to compute the JSD matrix which is
#     computationally costly. Subsequent runs just fix to FALSE to skip it
#     and read from file. If set to TRUE, not includes or excludes are allowed.
compute.dist = FALSE
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 
              # "Time0D_7D_4M" = whole  dataset, no excludes
              # "Time0D_7D" = exclude 4M dataset
              # "Time0D_7D_matched" = exclude 4M and samples not resurrected

###### STOP EDITING

# Read files ------------
# --- Read ASVs table
setwd("../6_finalfiles")
fileOTU="seqtable_readyforanalysis.csv"
ASV.table=read.table(fileOTU,sep="\t")
#colnames(ASV.table)[1:5]
#rownames(ASV.table)[1:5]
#dim(ASV.table)

# ..... read metadata. Samples present in metadata were those passing the filtering
setwd("../4_dada2/")
fileMD="metadata_Time0D-7D-4M_May2022.csv"
sample_md<-read.table(fileMD,sep="\t",header=TRUE)

# .... read distance if already available
if(compute.dist == FALSE){
  setwd("../7.1_classes")
  fileDist="Dist_JSD_Time0D-7D-4M.RDS"
  data.dist=readRDS(fileDist) 
}

# Clean data   ------
ASV.table=ASV.table[,colSums(ASV.table) > 0]
ASV.table=ASV.table[rowSums(ASV.table) > nreads, ]
which(is.na(ASV.table))

# .... exclude samples that didn't pass quality control (not present in samples metadata)
matched=match(row.names(ASV.table),sample_md$sampleid) #attr(data.dist,"Labels"))
ASV.table=ASV.table[!is.na(matched),]
dim(ASV.table)
if(compute.dist == FALSE){
  # ... double check that is ordered as the ASV table
  matched = match(rownames(ASV.table), names(data.dist))
  which(is.na(matched))  # none
  dim(as.matrix(data.dist)) # 4 more elements
  data.dist.tmp=as.dist(as.matrix(data.dist)[matched, matched])
  data.dist = data.dist.tmp # same here
  dim(as.matrix(data.dist)) # 4 more elements
}
  
# ... rebuild the metadata, it may have more samples
matched=match(row.names(ASV.table),sample_md$sampleid)
matched=matched[!is.na(matched)]
sample_md = sample_md[matched, ]

# Exclude data ----
if((length(exclude_exp) != 0) | (match_exp == TRUE)){
  sample_md_red = sample_md
  if(length(exclude_exp) != 0){ # exclude experiments
    for(level in exclude_exp){
      idx=which(sample_md_red$Experiment != level)
      sample_md_red = sample_md_red[idx, ]
    }
  }
  
  if(match_exp == TRUE){ # exclude communities that were not resurrected
    nchild = table(sample_md_red$parent)
    idx = which(nchild == 5)
    matched = match(sample_md_red$parent, names(nchild)[idx]) 
    sample_md_red = sample_md_red[!is.na(matched), ]
  }
  sample_md = sample_md_red #
  
  # ... look for the positions of remaining samples
  matched=match(sample_md$sampleid,row.names(ASV.table))
  matched=matched[!is.na(matched)]
  
  # ... reshape objects 
  ASV.tmp=ASV.table[matched,] # this tmp file is unnecessary, but easier to debug
  dim(ASV.tmp)
  ASV.table = ASV.tmp # this tmp file is unnecessary, but easier to debug
  if(compute.dist == FALSE){
    data.dist.tmp=as.dist(as.matrix(data.dist)[matched, matched])
    data.dist = data.dist.tmp # same here
  }
  ASV.table=ASV.table[,colSums(ASV.table) > 0]
  # ... rebuild the metadata, it may have more samples
  matched=match(row.names(ASV.table),sample_md$sampleid)
  matched=matched[!is.na(matched)]
  sample_md = sample_md[matched, ]
}


# Classification ----------
# --- Compute classification at different thresholds
# ..... this is a computationally costly step
if(compute.dist == TRUE){
  out_classes=otu_table_to_classes(t(ASV.table),Nclus = 25) # no distance provided, it will return JSD
  data.dist = out_classes$data.dist # only exists if a distance was computed
  # .... Write data.dist (uncomment if it  is the first time you compute it)
  setwd("../7.1_classes")
  fileOut="Dist_JSD_Time0D-7D-4M.RDS"
  saveRDS(data.dist,file=fileOut)
}else{
  out_classes=otu_table_to_classes(t(ASV.table),Nclus = 25, 
                                   distance = data.dist, euclidean = T) 
}

# .... Extract objects
nclusters = out_classes$nclusters
Kmax = which.max(nclusters)

# .... Plot results
setwd("../7.1_classes")
filename="Plot_PAM-IndexG1_orig-medoidsV2"
plot_calinski(nclusters,
              filename = filename,
              outputlabel = output.label) # check the maxima

# .... Extract partition
partition=clustering_to_partition(data.dist,Kmax = Kmax,
                                  outputlabel = output.label,
                                  print.file = TRUE)
sample_md[, output.label] = partition[,"PartId"]

# .... Print output
fileOut=unlist(strsplit(fileMD,".",fixed=T))[1]
fileOut=paste0(fileOut,"_",output.label,"_wJSDpart-all.csv")
write.table(sample_md,file = fileOut,sep="\t",quote=F,row.names = F)

