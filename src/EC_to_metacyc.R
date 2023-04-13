##################################################
# EC_to_metacyc.R
##################################################
# In this script I read the table with all EC numbers predicted by PiCRUST, and
# I compute the mean proportion for the different classes. It takes a list
# of the most significant EC numbers and creates a heatmap. It also transforms the
# code to be identified by MetaCyc omics dashboard for further analysis. The
# script was generated from the script pathways_to_heatmap.R
#
##################################################
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list=ls())
library(gplots)

# START EDITING -----------
# --- Set the files needed
file.meta="metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv" # metadata with partitions
fileEC="pred_metagenome_unstrat_descrip_2R.tsv" # EC table with descriptors
fileSigPath="Pathways_significant_effSize0.1.list" # List identifying significant pathways
# STOP EDITING -----------

# Load data -------------
# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the Final.Repository
setwd(dirSrc)
source("output_matrices.R")
source("functionArgsList.R")

#  --- Load metadata
setwd("../7.1_classes")
sample_metadata = read.table(file.meta,sep="\t",header=TRUE)

#  --- Load table pathways
setwd("../8_PICRUSt2/EC_metagenome_out")
EC.df = read.table(fileEC,sep="\t",header = TRUE,row.names = 1)
rownames(EC.df)[dim(EC.df)[1]]

#  --- Load list significant pathways
setwd("Stats")
EC.sig = scan(fileSigPath,what="character",sep="\n")

# Start analysis ---------
setwd(dirSrc)
setwd("../8_PICRUSt2/EC_metagenome_out")

# --- Remove unwanted samples
EC.descrip=as.vector(EC.df[,1]) # Store descriptors
names(EC.descrip)=rownames(EC.df)
dim(EC.df)
exclude="4M"
id.keep=as.character(sample_metadata$sampleid[which(sample_metadata$Experiment != exclude)])
exclude="Rep0.Class6"
id.keep=c(id.keep,as.character(sample_metadata$sampleid[which(sample_metadata$replicate.partition != exclude)]))
id.keep=id.keep[duplicated(id.keep)] # we retain only those that fulfilled both conditions
EC.df=EC.df[,id.keep]
dim(EC.df)

# ..... double check
#matched=match(colnames(EC.df),sample_metadata$sampleid)
#View(sample_metadata[matched,])

# --- Compute proportions
col.tot=colSums(EC.df)
EC.df.prop=apply(EC.df,MARGIN = 2, FUN=function(x){x/sum(x)})
EC.df.prop=as.data.frame(EC.df.prop)
colSums(EC.df.prop)

# --- Aggregate samples by class and compute the mean
matched=match(colnames(EC.df.prop),sample_metadata$sampleid) # find samples position in the metadata
which(is.na(matched) == TRUE) # double check all were matched
grouping=as.character(sample_metadata$exp.replicate.partition[matched]) # grouping var in the correct order
EC.aggr=aggregate(t(EC.df.prop),list(Class = grouping),# note transposition
                    function(x){mean(x,na.rm=T)})
rownames(EC.aggr)=EC.aggr$Class
# **Warning** the following string was created manually, and it depends
#             on the order of aggregation. Check rownames(EC.aggr)
# --- Class 6 excluded
newNames=c("Final.Rep1.Class1","Final.Rep1.Class2","Final.Rep2.Class1","Final.Rep2.Class2",
           "Final.Rep3.Class1", "Final.Rep3.Class2", "Final.Rep4.Class1", "Final.Rep4.Class2",
            "Starting.Class1", "Starting.Class2","Starting.Class3" ,"Starting.Class4", "Starting.Class5" )
# --- Class 6 included
#newNames=c("Starting.Class1", "Starting.Class2","Starting.Class3" ,"Starting.Class4", "Starting.Class5" ,
#           "Starting.Class6" ,"Final.Rep1.Class1","Final.Rep1.Class2","Final.Rep2.Class1","Final.Rep2.Class2",
#           "Final.Rep3.Class1", "Final.Rep3.Class2", "Final.Rep4.Class1", "Final.Rep4.Class2")
rownames(EC.aggr)=newNames
EC.aggr=subset(EC.aggr,select = -c(Class))
EC.id.2metacyc=unlist(lapply(strsplit(colnames(EC.aggr),fixed=TRUE,split = ':'),
                             FUN=function(x){x[2]}))
EC.aggr.2metacyc=EC.aggr
colnames(EC.aggr.2metacyc)=EC.id.2metacyc
fileOut="EC_mean_proportion_byClass_2metacyc.tsv"
write.table(t(EC.aggr.2metacyc),file=fileOut,quote=FALSE,
            sep="\t",col.names=NA)

# --- Extract only pathways of interest
matched=match(EC.sig,colnames(EC.aggr))
which(is.na(matched) == TRUE) # double check
EC.aggr.sig=as.matrix(EC.aggr[,matched])
# rownames(EC.aggr.sig)
class(EC.aggr.sig)

# --- Transform into a z-score
EC.aggr.zscore=apply(EC.aggr.sig, MARGIN=2, 
                       FUN = function(x){(x-mean(x))/sd(x)})

# --- Rename
matched=match(colnames(EC.aggr.zscore),names(EC.descrip))
colnames(EC.aggr.zscore)=EC.descrip[matched]

# --- Finally plot
setwd("Figures")
labelOut="SigPathways_StartVsFinal_byClass" # label for output
lmat=rbind(3:4,2:1);lhei=c(1.5,4);lwid=c(4,1.5) # some graphical parameters
output_matrices(mat = t(EC.aggr.zscore),name.mat = labelOut,plot.heatmap = TRUE,
              par.heatmap = list(#xlab="Pathways",ylab="Classes",
                                 scale = "none",#"row",
                                 margins = c(16,30),
                                 cexCol = 2,cexRow = 1,
                                 #sepwidth = c(0.2,0.2),
                                 #lmat=lmat,lhei=lhei,lwid=lwid, # doesn't work ¿?
                                 key.xlab = "Z-score mean proportion",
                                 keysize = 0.75,
                                 key.par = list(cex.main=2,cex.axis=1.2,
                                                cex.lab=1.5),
                                col = "bluered"),
              par.pdf = list(width=14,height=20))

