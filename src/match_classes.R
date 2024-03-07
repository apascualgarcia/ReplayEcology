###################################
## match_classes.R
###################################
# This scripts takes a set of samples
# and the identifiers of two partitions defined on them
# and it computes the overlap between both.
# For comparison of classes defined in two
# sets (e.g. time 0 vs time 7) see the
# script TraceSamples_from
###################################
## Author: Dr. Alberto Pascual-García
## Date Created: 2023-04-16
## MIT Licensed.
## Contact: apascualgarcia.github.io
#################################### 
##
## INPUT:
## OUTPUT:  
## DEPENDENCIES:
## USAGE:
#####################################
## Source libraries if needed
#setwd("/home/apascual/APASCUAL/Research/Programs/libraries/R")  
## load up our functions into memory
# source("some_runction.R")

## load up the packages we will need:  (uncomment as required)
# require(ggplots2)
# require(tidyverse)
# require(data.table)
#####################################  

# START EDITING -----------
# --- Set the files needed
file.meta="metadata_Time0D-7D-4M_May2022_wJSDpart-merged.csv" # metadata with partitions

# --- Select category to classify the distances and prepare structures
part1 = "Part_Time0D_17" # "Time0D_7D" 
part2 =  "Time0D_7D" # "Time0D_7D"
# "partition" #"time" #"experiment.location" # Time0D_7D # Part_Time0D_17
#"replicate.partition"# "Experiment" # "parent"

# --- Experiments. Either both partitions are defined onto the same experiment or
#     the first partition (exp1) is defined for 0D and the second (exp2) on 7D
exp1 = "0D"
exp2 = "7D"
# STOP EDITING -----------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

#  --- Load metadata
setwd("../7.1_classes")
sample_metadata <-read.table(file.meta,sep="\t",header=TRUE)
sample_metadata$time = as.character(sample_metadata$Experiment)
id.7D = grep("7D",sample_metadata$time)
sample_metadata$time[id.7D] = "7D"

# --- Identify levels

if(exp1 == exp2){
  mode = "single"
  id.exp = which(sample_metadata$time == exp1)
  sample_metadata = sample_metadata[id.exp,]
  levels1 = levels(as.factor(sample_metadata[ ,part1])) # independently  of the mode
  levels1 = levels1[!is.na(levels1)]
  levels2 = levels(as.factor(sample_metadata[ ,part2]))
  levels2 = levels2[!is.na(levels2)]
}else{
  mode = "double"
  id.exp1 = which(sample_metadata$time == exp1)
  id.exp2 = which(sample_metadata$time == exp2)
  id.exp = c(id.exp1, id.exp2)
  levels1 = levels(as.factor(sample_metadata[id.exp1,part1])) # independently  of the mode
  levels1 = levels1[!is.na(levels1)]
  levels2 = levels(as.factor(sample_metadata[id.exp2,part2])) # independently  of the mode
  levels2 = levels2[!is.na(levels2)]
  replicates = levels(as.factor(sample_metadata[id.exp2, "replicate"]))
  sample_metadata = sample_metadata[id.exp,]
  new.levels2 = c() # for this mode, we need to include replicates in level2
  for(replicate in replicates){
    for(level in levels2){
      new.level = paste(replicate,level,sep = ".")
      new.levels2 = c(new.levels2, new.level)
    }
  }
  levels2 = new.levels2
}


# --- Build a matrix to compute the overlap
Nlevels1 = length(levels1); Nlevels2 = length(levels2)
overlap = matrix(0, nrow = Nlevels1, ncol = Nlevels2)
rownames(overlap) = levels1; colnames(overlap) = levels2
Nsamples = 0

# --- Compute the overlap, differentiate between modes
if(mode == "single"){
  for(i in 1:dim(sample_metadata)[1]){ # simply follow the list and compare partitions for each sample
    part1.tmp = sample_metadata[i, part1]
    part2.tmp = sample_metadata[i, part2]
    if(is.na(part1.tmp) | is.na(part2.tmp)){
      next
    }
    Nsamples = Nsamples + 1
    overlap[part1.tmp, part2.tmp] = overlap[part1.tmp, part2.tmp] + 1
  }
}else{
  id.exp1 = which(sample_metadata$time == exp1)
  id.exp2 = which(sample_metadata$time == exp2)
  Ndiv = 0 # number of samples with a child diverging
  div_1 = 0 # with only one diverging
  div_2 = 0 # with a split, half and half
  samples_div = data.frame()
  for(i in id.exp1){ # follow the list for exp1
    part1.tmp = sample_metadata[i, part1]
    if(is.na(part1.tmp)){next}
    parent = sample_metadata[i, "parent"] # identify the parent
    child.id = which(sample_metadata[id.exp2, "parent"] == parent)
    if(length(child.id) == 0){next} # no child
    all.part2.tmp = c()
    all.rep.part2.tmp = c()
    for(child in child.id){
      id.child = id.exp2[child]
      part2.tmp = sample_metadata[id.child, part2]
      if(is.na(part2.tmp)){next}
      replicate.tmp = sample_metadata[id.child, "replicate"]
      rep.part2.tmp = paste(replicate.tmp, part2.tmp, sep = ".")
      overlap[part1.tmp, rep.part2.tmp] = 
                            overlap[part1.tmp, rep.part2.tmp] + 1
      all.part2.tmp = c(all.part2.tmp, part2.tmp)
      all.rep.part2.tmp = c(all.rep.part2.tmp, rep.part2.tmp)
    }
    diverge = length(unique(all.part2.tmp))
    if(diverge > 1){
      Ndiv = Ndiv + 1
      sample.div = c(as.character(parent), all.rep.part2.tmp)
      count = table(all.part2.tmp)
      if(max(count) == 3){
        div_1 = div_1 + 1
      }else{
        div_2 = div_2 + 1
      }
      if(Ndiv == 1){
        samples_div = as.data.frame(t(sample.div))
      }else{
        if(length(sample.div) < dim(samples_div)[2]){
          sample.div = c(sample.div, NA)
        }
        samples_div = rbind(samples_div, t(sample.div))
      }
    }
  }
}

# Print outputs -------
setwd("../7.2_match")
labelOut = paste0("exp1-",exp1,"_exp2-",exp2,
                  "_part1-",part1,"_part2-",part2)
fileMat = paste0("Matrix_",labelOut,".tsv")
write.table(overlap, file = fileMat, sep = "\t", quote = F)

if(mode == "double"){
  if(Ndiv > 0){
    fileStat = paste0("Ndiv_samples_",labelOut,".txt")
    sink(file = fileStat)
    mes = paste("Number of divergent samples = ", Ndiv, "\n")
    cat(mes)
    mes = paste("Diverging only one = ", div_1, "\n")
    cat(mes)
    mes = paste("Diverging two = ", div_2, "\n")
    cat(mes)
    sink()
    fileList = paste0("List_div_samples_",labelOut,".tsv")
    write.table(samples_div, file = fileList, quote = F, 
                sep = "\t", row.names = F)
  }
}