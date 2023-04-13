#############################
#   reshape_metadata.R
#############################
# This script reads the metadata table, fills
# the "replicate" factor and creates a new 
# factor combining replicate and SJD partition. It
# is incorporated in main_find_classes_exp-split.R
# so it could be converted in a function.
####################
# author = apascualgarcia.github.io
# May 17th, 2022. ETH-Zürich
####################

rm(list=ls())
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
setwd(dirSrc)
setwd("../7.1_classes")
fileIn = "metadata_Time0D-7D-4M_May2022_wJSDpart-split.csv"
sample_md= read.table(fileIn,sep="\t")

# --- Complete "replicate" for samples != 4M and create a single 
#     "partition" and "replicate.partition" factors
sample_md$partition=NA
sample_md$replicate.partition=NA

experiments=levels(sample_md$Experiment) 
for(level in experiments){
  cols=grep(colnames(sample_md),pattern = level)
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
  id=which(sample_md$Experiment == level)
  sample_md$replicate[id]=as.numeric(replicate) #as.character(sample_md$Experiment[id])
  sample_md$partition[id]=sample_md[id,cols]
  sample_md$replicate.partition[id]=paste(sample_md$replicate[id],
                                         sample_md[id,cols],sep=".")
}

# end function
write.table(sample_md,file = fileIn,sep="\t",quote=F)

