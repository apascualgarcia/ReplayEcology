##################################################
# represent_attractors.R
##################################################
# In this script I take the starting communities
# and the replicates and I represent them in 
# an ordination plot as in phyloseq_analysis.R, but
# now focusing on matched communities and representing
# how starting communities end up in final ones.
# Zürich, December 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list=ls())
library(phyloseq)
library(reshape2)
#library(RDPutils) 
library(ggplot2)
library(usedist)

# START EDITING -----------
# --- Set the files needed
file.meta = "metadata_Time0D-7D-4M_May2022_wJSDpart-merged_ext.csv" # new metadata with combinations of columns
file.attr = "Distances_from0D_to_4.1_4.2.tsv" # additional info about the attractors
fileOTU="seqtable_readyforanalysis.csv" # ASVs table
rep.id = "Rep4" # choose the replicate to use to confront

# --- Prepare plot features
axesby=c(1,2) # axes to represent

# .... Color of points
#colorby="DominantClass"; colorlab = "Final class"
#labelsColor = c("Goes to 1", "Goes to 2","Goes to both",
#                "1", "2")
colorby="exp.partition"; colorlab = "Local class"
labelsColor = c("Final / 1", "Final / 2", "Starting / 1", "Starting / 2","Starting / 3",
                "Starting / 4", "Starting / 5")
# .... vectors of colors you can use
# colors=c("#F8766D","#7CAE00","#00B0F6","red","#00C1A3")
#colors = c(2, 3, 4, 5 ,6)
#colors = c("blue", "green4","red", "green", "cyan", "grey" , "gold")
colors = c("chocolate4","chartreuse",
  "red1", "green4", "magenta1", "gold2","deepskyblue1", "black")
valuesColor = colors

shapeby="Experiment"; shapelab = "Communities"
sizeby="MaxReplicates"
#library(scales)
#show_col(hue_pal()(20))
# --- Labels also are manually set
labelsShape = c("Starting", "Final")
valuesShape = c(16,17)

# ... Represent trajectories? 
trajectories = T # If T a random subset of trajectories will be represented
# STOP EDITING -----------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

#  --- Load metadata
setwd("../7.1_classes")
sample_metadata = read.table(file.meta, sep = "\t",header=TRUE)
rownames(sample_metadata) = sample_metadata$sampleid

# --- Load information about the realation between starting communities and final attractor
setwd("../7.5_attractors")
df_attractors = read.table(file.attr, sep = "\t",header=TRUE)
rownames(df_attractors) = df_attractors$sample

# --- Read ASVs table
setwd("../6_finalfiles")
otu.in=read.table(fileOTU,sep="\t")
otu.in.t=t(otu.in)
#fileOut="seqtable_readyforanalysis.t.csv"
#write.table(otu.in.t,file=fileOut,sep="\t",quote=FALSE)
dim(otu.in)
otu.in=as.matrix(otu.in)
head(otu.in)[1:5,1:5]


# --- Match samples to those present in the attractors file
# .... we need to incorporate df_attractor info to sample_metadata
# .... Initialize new factors
sample_metadata$DominantClass = NA
sample_metadata$MaxReplicates = NA
sample_metadata$ClosestCentroid = NA

# .... Match both dataframes and add info for starting communities (NA for final)
matched = match(df_attractors$sample, sample_metadata$sampleid)
which(is.na(matched) == TRUE) # double check
sample_metadata$DominantClass[matched] = df_attractors$Class
sample_metadata$MaxReplicates[matched] = df_attractors$Nsame.class
sample_metadata$ClosestCentroid[matched] = df_attractors$id.DminToC

# .... Copy the info of final communities class id in a new factor (NA for starting)
id_final = which(sample_metadata$replicate == rep.id)
sample_metadata$DominantClass[id_final] = sample_metadata$replicate.partition[id_final]
sample_metadata$MaxReplicates[id_final] = 1


# .... Extract from metadata desired samples
matched = c(matched, id_final)
sample_metadata = sample_metadata[matched, ]

# --- Build phyloseq objects with final samples
otu.pseq = otu_table(as.matrix(otu.in), taxa_are_rows = FALSE)
df.pseq = sample_data(sample_metadata)

# -- Build the phyloseq object
treeholes = merge_phyloseq(otu.pseq, df.pseq) # If all dataset

# Rarefaction ---------
# --- Create some plots for the minimum acceptable level of rarefaction
set.seed(15082022) # Today's date 15/08/2022. Stored for reproducibility
otu.pseq.rar = rarefy_even_depth(otu.pseq, sample.size = 10000) # 10k is the minimum sampling sites, and the alpha diversity pattern is already there
#otu.pseq.rar = otu.pseq.rar[,which(colSums(otu.pseq.rar)!=0)]
#plot(colSums(otu.pseq.rar))
treeholes.rar=merge_phyloseq(otu.pseq.rar, df.pseq)  

# Ordination analysis ---------- 
# .... Select method, distance was read from file
method="PCoA" 
dist="jsd"

# .... These are all # (or matched, depending on your choice)
setwd("../7.5_attractors")
fileOrd=paste0("Ordination_jsd_PCoA_StartVs",rep.id,".RDS")
#treeholes.ord = ordinate(treeholes.rar, method=method,distance=dist) #dist.rar)#, "unifrac")
#saveRDS(treeholes.ord,file=fileOrd)
treeholes.ord=readRDS(fileOrd)

#stop()

# ... Plot all replicates
treeholes.rar@sam_data$DominantClass = as.factor(treeholes.rar@sam_data$DominantClass)

# Plot -----------------
# ... prepare labels
labaxes=paste0("Axes", paste(axesby[1], axesby[2], sep = "-"))
explainX = round(treeholes.ord$values[axesby[1],2] * 100, digits = 2)
explainY = round(treeholes.ord$values[axesby[2],2]* 100, digits = 2)
xlab = paste0("PCoA, component ", axesby[1]," [",explainX,"%]")
ylab = paste0("PCoA, component ", axesby[2]," [",explainY,"%]")
labelOut=paste0("StartVs",rep.id)

if(trajectories == T){
  set.seed("20042023") # seed 19, 20, 21 stands for the first two digits
  # take some points from starting communities
  id.start = which(sample_metadata$ExpCompact == "Starting")
  # randomly take a subset
  id.start = sample(id.start, size = length(id.start)/10, replace = F)
  sample.start = sample_metadata$sampleid[id.start]
  # identify the child
  id.end = which(sample_metadata$ExpCompact == "Final")
  # match parent and child
  matched = match(sample.start, sample_metadata$parent[id.end])
  sample.end = sample_metadata$sampleid[id.end[matched]]
  # take coordinates
  ord.coor = treeholes.ord$vectors[,axesby]
  matched = match(sample.start, rownames(ord.coor))
  ord.coor.start = ord.coor[matched,]; 
  colnames(ord.coor.start) = c("x_start","y_start") 
  matched = match(sample.end, rownames(ord.coor))
  ord.coor.end = ord.coor[matched,]
  colnames(ord.coor.end) = c("x_end","y_end")
  ord.coor.samp = as.data.frame(cbind(ord.coor.start, ord.coor.end))
}

# ... prepare output file
plotTitle=paste("Plot",method,"_",dist,"_By",colorby,"_",labelOut,"_",labaxes,".pdf",sep="")
pdf(file=plotTitle,width=8,height=6)
title=paste(method," of ",dist," distance",sep="")

# ... plot
p = plot_ordination(treeholes.rar, treeholes.ord, color = colorby,
                    shape=shapeby,#size=sizeby,
                    axes=axesby) #,
p = p + geom_point(size = 3, alpha = 0.7) + #ggtitle(title)+ 
  labs(color = colorlab, shape = shapelab)+xlab(xlab)+ylab(ylab)+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        title=element_blank(), #element_text(size=16),
        axis.text=element_text(size=15),
        strip.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=16))+
  scale_shape_manual(labels = labelsShape,
                     values = valuesShape)+
  scale_color_manual(labels = labelsColor,
                     values = valuesColor) #colors)
if(trajectories == T){
  p = p + geom_segment(data = ord.coor.samp,
                       aes(x = x_start, y = y_start,
                           xend = x_end, yend = y_end, alpha=1),# colour = class),
                       color= "grey50", # "grey50",
                       lineend = "round", linejoin = "mitre",
                       #curvature = 2,
                       arrow = arrow(ends = "last", # first last or both
                                     angle = 15,
                                     length=unit(0.2,"cm"),
                                     type = "closed"),
                       inherit.aes = FALSE,
                       show.legend = FALSE)
}
print(p)
dev.off()

