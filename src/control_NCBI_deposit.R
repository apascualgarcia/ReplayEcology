###################################
## control_NCBI_deposit.R
###################################
# We control here that the all samples finally
# used are present in the NCBI deposition file.
#
###################################
## Author: Dr. Alberto Pascual-García
## Date Created: 2023-06-05
## MIT Licensed.
## Contact: apascualgarcia.github.io
#################################### 
##
## INPUT: ASV table, NCBI metadata
## OUTPUT: Summary of matching
## DEPENDENCIES: Should be run from rstudio
## USAGE:
#####################################
rm(list= ls())

# START EDITING -----------
dirData = "../6_finalfiles/" # Relative to the folder where the code is
dirMD = "../4_dada2/"

fileMD="metadata_Time0D-7D-4M_May2022.csv"
fileASV = "seqtable_readyforanalysis.csv"
fileNCBI = "NCBI_df.csv"

nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included

# STOP EDITING -----------

# --- Set the main directory, this must be run in rstudio to work
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

# ---  Read data
setwd(dirData)
ASV.table = read.table(fileASV)
head(ASV.table)[1:5,1:5]
NCBI.df = read.table(fileNCBI, sep = ",", header = T)
head(NCBI.df)[1:5,1:5]

# ..... read metadata
setwd(dirMD)
fileMD="metadata_Time0D-7D-4M_May2022.csv"
sample_md<-read.table(fileMD,sep="\t",header=TRUE)

# Clean data   ------
# check that indeed the table is already clean
dim(ASV.table)
ASV.table=ASV.table[,colSums(ASV.table) > 0]
ASV.table=ASV.table[rowSums(ASV.table) > nreads, ]
dim(ASV.table) # nothing changed
which(is.na(ASV.table))

# .... exclude samples that didn't pass previous controls (e.g. core communities, not present in samples metadata)
matched=match(row.names(ASV.table),sample_md$sampleid) #attr(data.dist,"Labels"))
ASV.table=ASV.table[!is.na(matched),]
dim(ASV.table) # nothing changed

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

  ASV.table=ASV.table[,colSums(ASV.table) > 0]
  # ... rebuild the metadata, it may have more samples
  matched=match(row.names(ASV.table),sample_md$sampleid)
  matched=matched[!is.na(matched)]
  sample_md = sample_md[matched, ]
}

# compare to NCBI table -----
dim(NCBI.df)

# --- create a sampleid matching the ones we use in the tables
NCBI.df$sampleid = sapply(NCBI.df$library_ID, 
                          FUN = function(x){sub("(.+?\\.)(.*)", "\\2", x)})
# ...... I observed some names where inconsistent
id.inconsistent = grep("day", NCBI.df$sampleid, value = F)
list.inconsistent = grep("day", NCBI.df$sampleid, value = T)
NCBI.df$library_ID[id.inconsistent] = paste(NCBI.df$library_ID[id.inconsistent],
                                            "X99", sep = ".")
# ...... Now repeat
NCBI.df$sampleid = sapply(NCBI.df$library_ID, 
                          FUN = function(x){sub("(.+?\\.)(.*)", "\\2", x)})

# ..... Finally match

matched = match(rownames(ASV.table), NCBI.df$sampleid)
list.notfound = rownames(ASV.table)[is.na(matched)]

# Report results ----
dim(NCBI.df) # 1426    7
dim(ASV.table) #  1375 1458
dim(sample_md) # 1375    9
length(which(sample_md$Experiment == "4M")) # 0
length(which(sample_md$Experiment == "0D")) # 275
length(which(sample_md$Experiment == "7D_rep1")) # 275
length(which(sample_md$Experiment == "7D_rep2")) # 275
length(which(sample_md$Experiment == "7D_rep3")) # 275
length(which(sample_md$Experiment == "7D_rep4")) # 275

setwd(dirData)
file.notfound = "Samples_notFound.list"
write.table(list.notfound, file = file.notfound, row.names = F, quote = F)

file.inconsistent = "Samples_inconsistent_libraryID.list"
write.table(list.inconsistent, file = file.inconsistent, row.names = F, quote = F)
