##################################################
# pathways_to_heatmap.R
##################################################
# In this script I read the table with all KEGG pathways, a list identifying
# those that were significant in any comparison between community classes
# (performed in STAMP) and the metadata, and I compute the proportion of genes
# that each pathway Final.Represent for each sample, to then calculate the means across
# all samples belonging to the same community class. I finally Final.Represent 
# results in a heatmap.
#
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list=ls())
library(gplots)

# START EDITING -----------
# --- Set the files needed
file.meta="metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv" # metadata with partitions
fileKEGG="KEGGs_vs_samples.L3.tsv" # KEGG table
fileSigPath="Pathways_significant_FinalClass1VsClass2_effSize0.1.list" #"Pathways_significant_effSize0.2.list" # List identifying significant pathways
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
setwd("../8_PICRUSt2/KEGG_pathways_out")
kegg.df = read.table(fileKEGG,sep="\t",header = TRUE,row.names = 1)

#  --- Load list significant pathways
setwd("Stats")
kegg.sig = scan(fileSigPath,what="character",sep="\n")

# Start analysis ---------
setwd(dirSrc)
setwd("../8_PICRUSt2/KEGG_pathways_out")

# --- Remove unwanted samples
dim(kegg.df)
exclude="4M"
id.keep=as.character(sample_metadata$sampleid[which(sample_metadata$Experiment != exclude)])
exclude="0D"
id.keep=c(id.keep,as.character(sample_metadata$sampleid[which(sample_metadata$Experiment != exclude)]))
#exclude="Rep0.Class6"
#id.keep=c(id.keep,as.character(sample_metadata$sampleid[which(sample_metadata$replicate.partition != exclude)]))
id.keep=id.keep[duplicated(id.keep)] # we retain only those that fulfilled both conditions
kegg.df=kegg.df[,id.keep]
dim(kegg.df)

# ..... double check
#matched=match(colnames(kegg.df),sample_metadata$sampleid)
#View(sample_metadata[matched,])

# --- Compute proportions
col.tot=colSums(kegg.df)
kegg.df.prop=apply(kegg.df,MARGIN = 2, FUN=function(x){x/sum(x)})
kegg.df.prop=as.data.frame(kegg.df.prop)
colSums(kegg.df.prop)

# --- Aggregate samples by class and compute the mean
matched=match(colnames(kegg.df.prop),sample_metadata$sampleid) # find samples position in the metadata
which(is.na(matched) == TRUE) # double check all were matched
grouping=as.character(sample_metadata$exp.replicate.partition[matched]) # grouping var in the correct order
kegg.aggr=aggregate(t(kegg.df.prop),list(Class = grouping),# note transposition
                    function(x){mean(x,na.rm=T)})
rownames(kegg.aggr)=kegg.aggr$Class
# **Warning** the following string was created manually, and it depends
#             on the order of aggregation. Check rownames(kegg.aggr)
# --- 4M and Class 6 excluded
#newNames=c("Final.Rep1.Class1","Final.Rep1.Class2","Final.Rep2.Class1","Final.Rep2.Class2",
#           "Final.Rep3.Class1", "Final.Rep3.Class2", "Final.Rep4.Class1", "Final.Rep4.Class2",
#            "Starting.Class1", "Starting.Class2","Starting.Class3" ,"Starting.Class4", "Starting.Class5" )
# --- 4M excluded and Class 6 included
#newNames=c("Starting.Class1", "Starting.Class2","Starting.Class3" ,"Starting.Class4", "Starting.Class5" ,
#           "Starting.Class6" ,"Final.Rep1.Class1","Final.Rep1.Class2","Final.Rep2.Class1","Final.Rep2.Class2",
#           "Final.Rep3.Class1", "Final.Rep3.Class2", "Final.Rep4.Class1", "Final.Rep4.Class2")
# --- 4M and 0D excluded
#newNames=c("Final.Rep1.Class1","Final.Rep1.Class2","Final.Rep2.Class1","Final.Rep2.Class2",
#           "Final.Rep3.Class1", "Final.Rep3.Class2", "Final.Rep4.Class1", "Final.Rep4.Class2")
newNames=c("Rep. 1, Class 1","Rep. 1, Class 2","Rep. 2, Class 1","Rep. 2, Class 2",
           "Rep. 3, Class 1", "Rep. 3, Class 2", "Rep. 4, Class 1", "Rep. 4, Class 2")
rownames(kegg.aggr)=newNames
kegg.aggr=subset(kegg.aggr,select = -c(Class))

# --- Extract only pathways of interest
matched=match(kegg.sig,colnames(kegg.aggr))
which(is.na(matched) == TRUE) # double check
kegg.aggr.sig=as.matrix(kegg.aggr[,matched])
rownames(kegg.aggr.sig)
class(kegg.aggr.sig)

# --- Transform into a z-score
kegg.aggr.zscore=apply(kegg.aggr.sig, MARGIN=2, 
                       FUN = function(x){(x-mean(x))/sd(x)})
# --- Finally plot
setwd("Figures")
#labelOut="SigPathways_StartVsFinal_byClass" # label for output
labelOut="SigPathways_FinalOnly_byClass_effSize0.05_Ver" # label for output
lmat=rbind(c(4,3),c(2,1));lhei=c(1.5,4);lwid=c(4,1.5) # some graphical parameters
# .... vertical
output_matrices(mat = t(kegg.aggr.zscore),name.mat = labelOut,plot.heatmap = TRUE,
              par.heatmap = list(#xlab="Pathways",ylab="Classes",
                                 scale = "none",#"row",
                                 margins = c(14,40),
                                 #srtRow = 180,
                                 cexCol = 2,cexRow = 1.8,
                                 #sepwidth = c(0.2,0.2),
                                 #lmat=lmat, lhei=lhei,lwid=lwid, # doesn't work ¿?
                                 key.title = "Mean proportion",
                                 key.xlab = "Mean prop. Z-score",
                                 keysize = 0.8,
                                 key.par = list(cex.main=1,cex.axis=2.5,
                                                mar=c(14,6,6,0),
                                                mgp=c(4,2,0),
                                                cex.lab=2.5),
                                col = "bluered"),
              par.pdf = list(width=13,height=26))
# THis doesn't work:
# .... horizontal
labelOut="SigPathways_FinalOnly_byClass_effSize0.05_Hor" # label for output

output_matrices(mat = (kegg.aggr.zscore),name.mat = labelOut,plot.heatmap = TRUE,
                par.heatmap = list(#xlab="Pathways",ylab="Classes",
                  scale = "none",#"row",
                  margins = c(40,14),
                  cexCol = 1.6,cexRow = 2,
                  srtCol = 50,
                  #sepwidth = c(0.2,0.2),
                  #lmat=lmat, lhei=lhei,lwid=lwid, # doesn't work ¿?
                  key.title = "Mean proportion",
                  key.xlab = "Mean prop. Z-score",
                  keysize = 0.8,
                   key.par = list(cex.main=1,cex.axis=2,
                  #                mar=c(14,6,6,0), # should be within key.par
                  #                mgp=c(4,2,0),
                                  cex.lab=2.5),
                  col = "bluered"),
                par.pdf = list(width=26,height=13))

