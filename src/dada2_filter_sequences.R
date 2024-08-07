# INSTALL DADA2 -----------------------------------------------------------

#NB TO USE DADA 3.14 AND REPRODUCE THE PIPELINE EXACTLY, YOU NEED TO FIRST INSTALL AND RUN R VERSION 4.1 

#unhash to install dada2 3.14
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.14")

# GETTING READY: LOAD DADA2 ---------------------------------------------------

#loads DADA2 and returns the package version you're running (should be 3.14 for complete reproducibility)
library(dada2); packageVersion("dada2")

# GETTING READY: READ IN FILES --------------------------------------------

path <- '2_demultiplexed/corrected' #directory containing the demultiplexed fastq files that were  deposited at NCBI under BioProject accession number PRJNA989519, alongside the other fastq files from other projects in the lab that were originally processed alongside them
path
fns <- list.files(path, pattern="fastq.gz", recursive = T) # CHANGE if different file extensions
length(fns)

# GETTING READY: EXRACT SAMPLE NAMES --------------------------------------

#write file names to sample.names
sample.names<-sub("^(.*)[.].*", "\\1", fns)

#remove the extension
sample.names<-sub("^(.*)[.].*", "\\1", sample.names)

#replace forward slash with _
sample.names<-gsub("\\/", "_", sample.names)

# # INSPECT READ QUALITY PROFILES -------------------------------------------
# 
#plot read quality profiles
qualplot_agg<-plotQualityProfile(file.path(path,fns), aggregate=T)
#save to 3_dada2 folder to save time later
saveRDS(qualplot_agg,'4_dada2/qualplot_agg.RDS')

# FILTER AND TRIM: PLACE FILTERED FILES IN SUB-DIRECTORY ---------------------------------------------------------

filts<-file.path('3_filtered', paste0(sample.names, "_filt.fastq.gz"))
filts
length(filts)

# FILTERING ---------------------------------------------------------------

out <- filterAndTrim(file.path(path,fns), filts,
                     truncLen=240,
                     maxN=0, maxEE=1, truncQ=11, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
#save output of filtering to RDS
saveRDS(out, '3_filtered/filteredandtrimmed.RDS')

#read output of filtering in (this allows you to avoid running filterAndTrim() again
out<-readRDS('3_filtered/filteredandtrimmed.RDS')
