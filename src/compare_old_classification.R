###################################
#    compare_classifications.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-05-27
# Script Name: compare_classifications.R  
# Script Description: This script takes two classifications
# and extract some statistics comparing them.
rm(list = ls())
# START EDITING ------------

file.Meta = "metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.Part = "Partition_PamClustering_SamplesTime0_ShannonJensen.vec"

# SET WORKING DIRECTORY -----------------------------
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirMeta=paste(this.dir,"/7.3_phyloseq/",sep="") # Dir of metadata
dirOutput=paste(this.dir,"/7.1_classes/",sep="") # Dir of output data CHANGE OUTPUT BY THE CORRECT DIRECTORY


# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse", "stringr", "ggplot2","psych","boot","gplots") # list of packages to load
n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}


# SOURCE FUNCTIONS ---------
# scripts <- c("script.R") # list of scripts to load
# 
# n_scripts <- length(scripts) # count how many scripts are required
# 
# for(n in 1:n_scripts){
#   cat("Loading script #", n, " of ", n_scripts, "... Currently Loading: ", scripts[n], "\n", sep = "")
#   lib_load <- paste("source(\"",scripts[n],"\")", sep = "") # create string of text for loading each library
#   eval(parse(text = lib_load)) # evaluate the string to load the script
# }

# STOP EDITING --------------------

# READ DATA --------------
# ..... read metadata. Samples present in metadata were those passing the filtering
setwd(dirMeta)
sample_md.all <-read.table(file = file.Meta, sep="\t", header=TRUE)
head(sample_md.all)[1:5,1:5]

setwd(this.dir)
setwd(dirOutput)
PartSJD.old = read.table(file.Part,row.names=1 )#,comment.char = "#",col.names = "PartId",row.names=1 )#,skip=1)#,comment.char = "#")
colnames(PartSJD.old)="PartId"

# MATCH DATA --------------

matched = match(sample_md.all$sampleid, rownames(PartSJD.old))
sample_md.time0 = sample_md.all[!is.na(matched),]
matched = matched[!is.na(matched)]
PartSJD.old = PartSJD.old[matched, ] 

# --- select partition
sel.part = "Part_Time0D_6"

PartSJD.new = sample_md.time0[, sel.part] 

matched.class.df = data.frame(PartSJD.old, PartSJD.new)
rownames(matched.class.df) = sample_md.time0$sampleid

# PAIRWISE DATA ---------

Nelem = dim(matched.class.df)[1]
row_names = c()
Nold = 0
Nnew = 0
Nboth = 0
for(i in 1:(Nelem - 1)){
  for(j in (i+1):(Nelem)){
    row_names = c(row_names, 
                  paste(rownames(matched.class.df)[i],
                        rownames(matched.class.df)[j]))
    row_names
    same_clus_old = 0
    same_clus_new = 0
    key = F
    if(matched.class.df$PartSJD.old[i] == matched.class.df$PartSJD.old[j]){
      same_clus_old = 1
      Nold = Nold + 1
      key = T
    }
    if(matched.class.df$PartSJD.new[i] == matched.class.df$PartSJD.new[j]){
      same_clus_new = 1
      Nnew = Nnew + 1
      if(key == T){Nboth = Nboth + 1}
    }
    class.tmp = cbind(same_clus_old, same_clus_new)
    if((i == 1)& (j == 2)){
      class.df = data.frame(class.tmp)
      next
    }
    class.df = rbind(class.df, class.tmp)
  }
}
rownames(class.df) = row_names

Ntot = Nelem * (Nelem - 1)/2
Ne = (Nold*Nnew + (Ntot - Nold)*(Ntot - Nnew))/Ntot
No = Nboth + (Ntot - Nold - Nnew + Nboth)

kappa = (No - Ne)/(Ntot - Ne)
(kappa)
(cohen.kappa(class.df))

heatmap.2(as.matrix(matched.class.df),
          Rowv = F,Colv = F,
          dendrogram = "none", trace = "none",
          col = c("1"="red","2" = "green4", "3"="gold",
                  "4"="pink", "5"="blue", "6" = "gray"))
          