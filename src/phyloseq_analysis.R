##################################################
# phyloseq_analysis.R
##################################################
# In this script I perform the typical pipeline of phyloseq analysis, to
# study phylogenetic distributions of taxa in the different samples, alpha diversity
# quantification, and dimensionality reduction, for treeholes samples.
#
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################
rm(list=ls())
library(phyloseq)
library(reshape2)
#library(RDPutils) 
library(ggplot2)
library(usedist)
library(plyr)
library(stringr)

# START EDITING -----------
# --- Set the files needed
file.meta = "metadata_Time0D-7D-4M_May2022_wJSDpart-merged.csv" # metadata with partitions
file.meta.out = "metadata_Time0D-7D-4M_May2022_wJSDpart-merged_ext.csv" # new metadata with combinations of columns
fileOTU = "seqtable_readyforanalysis.csv" # ASVs table
fileTaxonomy = "taxa_wsp_readyforanalysis.csv" # taxonomy
fileDist = "Dist_JSD_Time0D-7D-4M.RDS" # beta diversity distance

# select if matched or total dataset, will be processed below
set="all" #"matched" # or all

# select the name of the partition to create the barplots
sel_partition = "partition" # defaults to "partition", which is the clustering of each set made independently
                                    # the factor selected will override "partition"
# select how points in the pcoa will be differentiated
colorby="exp.partition" # "exp.partition" # "partition" "parent" "Location" replicate.partition
                      # if you selected in sel_partition a factor different than "partition" and
                      # you want to color by that factor, you can set here "partition"
shapeby="Time0D_7D" # keep empty if you don't want to change shapes
# STOP EDITING -----------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

#  --- Load metadata
setwd("../7.1_classes")
sample_metadata <- import_qiime_sample_data(file.meta) #read.table(file.meta,sep="\t",header=TRUE)
sample_metadata[, "partition"] = sample_metadata[, sel_partition]
sample_metadata$replicate = as.factor(paste("Rep",sample_metadata$replicate,sep=""))
sample_metadata$partition = as.factor(paste("Class",sample_metadata$partition,sep=""))
sample_metadata$replicate.partition = as.factor(paste(sample_metadata$replicate,
                                                      ".",sample_metadata$partition,sep=""))
experiments = unique(sample_metadata$Experiment)
sample_metadata$ExpCompact = NA
for(exp in experiments){
  id = which(sample_metadata$Experiment == exp)
  if(exp == "4M"){
    sample_metadata$ExpCompact[id] = "Evolved"
  }else if(exp == "0D"){
    sample_metadata$ExpCompact[id] = "Starting"
  }else{
    sample_metadata$ExpCompact[id] = "Final"
  }
}
sample_metadata$exp.replicate.partition=as.factor(paste(sample_metadata$ExpCompact,".",sample_metadata$replicate.partition,sep=""))
sample_metadata$exp.partition=as.factor(paste(sample_metadata$ExpCompact,".",sample_metadata$partition,sep=""))
write.table(sample_metadata,file = file.meta.out,sep="\t",quote=F,row.names = F) # overwrite old metadata

# --- Load beta diversity distance
#dist=readRDS(fileDist) # finally used the one computed internally in ordinate, see below

# --- Read ASVs table

setwd("../6_finalfiles")
otu.in=read.table(fileOTU,sep="\t")
otu.in.t=t(otu.in)
#fileOut="seqtable_readyforanalysis.t.csv"
#write.table(otu.in.t,file=fileOut,sep="\t",quote=FALSE)
dim(otu.in)
otu.in=as.matrix(otu.in)
head(otu.in)[1:5,1:5]

# --- Remove samples with less than 10K
otu.in = otu.in[which(rowSums(otu.in) >= 10000), ]
dim(otu.in)
otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = FALSE)

# --- Load taxonomy
taxonomy=read.table(fileTaxonomy,sep="\t",row.names = 1,header=TRUE) 
tax.pseq = tax_table(as.matrix(taxonomy))

# --- Finally build the phyloseq object
treeholes=merge_phyloseq(otu.pseq,sample_metadata,tax.pseq) # If all dataset
#str(treeholes)
#treeholes=merge_phyloseq(otu.pseq,sample_metadata.matched,tax.pseq) # If only matched

# Rarefaction ---------
# --- Create some plots for the minimum acceptable level of rarefaction
set.seed(15082022) # Today's date 15/08/2022. Stored for reproducibility
otu.pseq.rar = rarefy_even_depth(otu.pseq, sample.size = 10000) # 10k is the minimum sampling sites, and the alpha diversity pattern is already there
#otu.pseq.rar = otu.pseq.rar[,which(colSums(otu.pseq.rar)!=0)]
#plot(colSums(otu.pseq.rar))
treeholes.rar=merge_phyloseq(otu.pseq.rar,sample_metadata,tax.pseq) 
#str(treeholes.rar)

# Prune undesired samples -----
exclude="4M"
treeholes.rar=prune_samples(treeholes.rar@sam_data$Experiment != exclude, treeholes.rar)
#str(treeholes.rar)
# The following lines subset the distance computed with the in house jds function
# and read above. However, it works worse than the distance computed internally
# by the ordinate function, so there may be some bug when subsetting
#id.notexclude=which(sample_metadata$Experiment != exclude)
#samples_2keep = as.character(sample_metadata$sampleid[id.notexclude])
#dist.rar = dist_subset(dist,samples_2keep) 

# select now matched or all dataset
if(set == "all"){
  labelOut=set
}else{
  # ......These are only those matched and their replicas
  # APG The commented lines should be rewritten if matched samples are considered
  samples_list = treeholes.rar@sam_data$sampleid
  parents_list = treeholes.rar@sam_data$parent
  ntimes = count(parents_list)
  length(which(ntimes$freq == 5)) # 275
  parents_list_true = as.character(ntimes$x[ntimes$freq == 5])
  matched = match(parents_list, parents_list_true)
  keep.samples = as.character(samples_list[!is.na(matched)])
  
  #matched=match(treeholes.rar@sam_data$parent,Part.SJD.matched0$V1) # look those for those that have no match in a matched partition
  #keep.samples=!is.na(matched) # create a logical vector that will say true if the sample is to be kept
  #length(which(keep.samples==TRUE))
  treeholes.rar.matched=prune_samples(keep.samples, treeholes.rar)
  #print(set)
  labelOut=set
  treeholes.rar=treeholes.rar.matched
}
  
### From here double check what is matched and what is not
#dim(treeholes.rar@otu_table)

# Top Taxa ---------------
# --- Identify the top Ntop genus and families for each replica among all its samples 
#         (the Ntop+1 will be the remainder, labeled "others")
Ntop=20
lev.rep=as.list(levels(treeholes.rar@sam_data$replicate))
Nrep=length(lev.rep)
topESVs=c()
for (i in 1:Nrep){ # The top 15 in all replicates, and merge them together
  #rep=as.character(lev.rep[i])
  topTmp=names(sort(taxa_sums(subset_samples(treeholes.rar, replicate==lev.rep[[i]])), TRUE)[1:Ntop])
  topESVs=c(topESVs,topTmp)
}
topESVs=unique(topESVs)
taxTop = cbind(tax_table(treeholes.rar), family = NA, genus=NA, 
               OTUs=NA,ASV_genera=NA,ASV_genera_order=NA) #species = "other species") 
# I get an error in this line about the memory used, it may be a bug in the usage of tax_table, it happens
taxTop[topESVs, "family"] <- as(tax_table(treeholes.rar)[topESVs, "Family"], "character") # It happens when as.character is used
taxTop[topESVs, "genus"] <- as(tax_table(treeholes.rar)[topESVs, "Genus"],  "character") # solved updating packages in branca
taxTop[topESVs, "OTUs"] <- as(tax_table(treeholes.rar)[topESVs, "Species"],  "character") 
taxTop[topESVs, "ASV_genera"]  <- as(paste(topESVs,
                                           tax_table(treeholes.rar)[topESVs, "Genus"],sep="_"),
                                     "character")
taxTop[topESVs, "ASV_genera_order"]  <- as(paste(topESVs,
                                           tax_table(treeholes.rar)[topESVs, "Genus"],
                                           tax_table(treeholes.rar)[topESVs, "Order"],sep=" / "),
                                     "character") 
id.ANPR=grep("Allorhizobium",taxTop[,"ASV_genera_order"])
if(length(id.ANPR)>0){ # If the ANPR clade is in the list, we want to change the name
  strtmp=taxTop[id.ANPR,"ASV_genera_order"]
  strtmpsplit=unlist(strsplit(strtmp,split = "/",fixed = TRUE))
  newstr=paste(strtmpsplit[1],"ANRP clade",strtmpsplit[3],sep=" / ")
  taxTop[id.ANPR,"ASV_genera_order"]=newstr
}
tax_table(treeholes.rar) <- tax_table(taxTop)
rownames(tax_table(treeholes.rar))[1:5]
# Merge samples ------------
# Merge samples belonging to the same replica
treeholes.rar.byRep =merge_samples(treeholes.rar, "replicate")
sample_data(treeholes.rar.byRep)$replicate <- levels(sample_data(treeholes.rar)$replicate)
treeholes.rar.byRep = transform_sample_counts(treeholes.rar.byRep, function(x) 100 * x/sum(x))

# Same by partitions
treeholes.rar.byPar =merge_samples(treeholes.rar, "replicate.partition")
sample_data(treeholes.rar.byPar)$replicate.partition <- levels(sample_data(treeholes.rar)$replicate.partition)
treeholes.rar.byPar = transform_sample_counts(treeholes.rar.byPar, function(x) 100 * x/sum(x))

# Same without normalization to plot species
treeholes.rar.byParSp =merge_samples(treeholes.rar, "replicate.partition")
sample_data(treeholes.rar.byParSp)$replicate.partition <- levels(sample_data(treeholes.rar)$replicate.partition)

# .... abundances per replicate
treeholes.rar.byRep.19 = prune_taxa(topESVs, treeholes.rar.byRep)
treeholes.rar.byPar.19 = prune_taxa(topESVs, treeholes.rar.byPar)

# Plot palette ---------------------------
setwd("../7.3_phyloseq")

# :: Build a palette
library(RColorBrewer) # I'm gonna create a palette
selected=set
#n = length(levels(sample_metadata$Location)) # To select a number of colors (not sure where this is used)
qual_col_pals = brewer.pal.info[(brewer.pal.info$category == 'qual'),] # Get qualitative palettes

col_vector = c("deepskyblue1","green4","red","green","cyan","grey50","gold","black")
#col_vector2 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col_vector2 = c(brewer.pal(4,"Set2"),brewer.pal(4,"Dark2"))
col_vector2 = c(brewer.pal(6,"Spectral"),brewer.pal(4,"Accent"), brewer.pal(12,"Paired"),
                brewer.pal(8,"Set2"),brewer.pal(8,"Dark2")) # Select manually
#col_vector2=col_vector[c(1:5,7:length(col_vector))] # Select between the list jumping
#col_vector2[4]="cyan"
#col_vector2=col_vector[seq(from=1,to=length(col_vector),by=2)] # Select between the list jumping
col_vector = c(col_vector,col_vector2)
# ... And a vector for the partitions
colorCodes = c("Class1"="red1", "Class2"="green4", "Class3"="magenta1", "Class4"="grey50","Class5"="deepskyblue1", "Class6"="gold")
colorCodes2 = c("Rep0.Class1"="red", "Rep0.Class2"="green", "Rep0.Class3"="magenta1", "Rep0.Class4"="grey50",
               "Rep0.Class5"="deepskyblue1", "Rep0.Class6"="gold","Rep1.Class1"="red1","Rep1.Class2"="green1",
               "Rep2.Class1"="red2","Rep2.Class2"="green2","Rep3.Class1"="red3","Rep3.Class2"="green3",
               "Rep4.Class1"="red4","Rep4.Class2"="green4") # barplot bar species coloured by partition


# Plot bars biodiversity---------
plotTitle=paste0("BarPlot_diversity-replicates_",selected,"_Par",sel_partition,"_byGenus.pdf") # All samples
#plotTitle=paste0("BarPlot_diversity-replicates_",selected,"_bySpecies.pdf") # All samples
#plotTitle=paste("BarPlot_diversity-replicates_",selected,"_",labelOut,"_","Time0restrictPar_bySpecies.pdf",sep="") # Only matched

tmp.data.frame=psmelt(treeholes.rar.byPar) # this is the dataframe that the function plot_bar creates internally,
                                           # and then it selects the column to plot and hence the label, I
                                           # create a vector from that to add NA as a factor
species.vec=tmp.data.frame$ASV_genera_order # get the vector that will lead to the bars and the labels
#genus.vec=tmp.data.frame$genus
species.vec=addNA(species.vec) # add NA as a factor
#genus.vec=addNA(genus.vec)
levels(species.vec)=c(levels(species.vec),"other") # Add the new name, now other is at the end
#levels(genus.vec)=c(levels(genus.vec),"other")
newlabels=levels(species.vec)
#newlabels=levels(genus.vec)

pdf(file=plotTitle,width=18,height=10)
title = " "
#plot_bar(treeholes.rar.byRep,"replicate",fill="OTUs",title=title)+
plot_bar(treeholes.rar.byRep,"replicate",fill="ASV_genera_order",title=title)+
  ylab("Relative abundance (%)")+coord_flip()+ 
  scale_fill_manual(values = col_vector,labels=newlabels)+
  theme(axis.title = element_text(size=18),
        title=element_text(size=16),
        axis.text=element_text(size=15),
        strip.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=14))
  
dev.off()

# --- By partition
#plotTitle=paste("BarPlot_diversity-partitions_",selected,"_",labelOut,"Time0restrictPar_bySpecies.pdf",sep="") 
plotTitle=paste0("BarPlot_diversity-partitions_",selected,"_Par",sel_partition,"_byGenus.pdf") 

tmp.data.frame=psmelt(treeholes.rar.byRep) # this is the dataframe that the function plot_bar creates internally,
# and then it selects the column to plot and hence the label, I was tryinng to create a vector out of that to change NA
species.vec=tmp.data.frame$ASV_genera_order # get the vector that will lead to the bars and the labels
species.vec=addNA(species.vec) # add NA as a factor
levels(species.vec)=c(levels(species.vec),"Other") # Add the new name, now other is at the end
newlabels=levels(species.vec)

pdf(file=plotTitle,width=18,height=10)
title = " "
plot_bar(treeholes.rar.byPar,"replicate.partition",fill="ASV_genera_order",title=title)+ #,facet_grid=~partition)+
  ylab("Relative abundance (%)")+ 
  labs(fill = "ASV id. / Genera / Order")+
  scale_fill_manual(values = col_vector,labels=newlabels)+coord_flip()+
  theme(axis.title = element_text(size=18),
        title=element_text(size=16),
        axis.text=element_text(size=15),
        strip.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=14))

dev.off()

#plot_bar(treeholes.rar.byParSp,"species", fill="replicate.partition",title=title)+ #,facet_grid=~partition)+
#  ylab("Total abundance")+ scale_fill_manual(values = colorCodes2)+coord_flip()


# Ordination analysis ---------- 
# .... Select method, distance was read from file
method="PCoA" 
dist="jsd"

# .... These are all # (or matched, depending on your choice)
#treeholes.ord = ordinate(treeholes.rar, method=method,distance=dist) #dist.rar)#, "unifrac") # takes more than 1h for JSD
fileOrd=paste0("Ordination_jsd_PCoA_",selected,".RDS")
#saveRDS(treeholes.ord,file=fileOrd)
treeholes.ord=readRDS(fileOrd)

# ......These are only time 0
treeholes.rar.time0=subset_samples(treeholes.rar,replicate=="Rep0")
#treeholes.ord.time0=ordinate(treeholes.rar.time0, method=method,distance=dist)

# ... Plot all replicates

if(colorby=="partition"){
   Npart = length(levels(treeholes.rar@sam_data$partition))
   usecolor=colorCodes[1:Npart]
   if(length(shapeby) == 0){
     colorlab="Class"
   }else{
     colorlab = "Local class"
     shapelab = "Global class"
   }
}else{
  if(colorby == "replicate.partition"){
    if(length(shapeby) == 0){
      colorlab="Replicate / Class"
    }else{
      colorlab = "Set / Local class"
      shapelab = "Global class"
    }
    usecolor=col_vector
  }else if(colorby == "Location"){
    colorlab="Location"
    usecolor=col_vector
  }else if(colorby == "exp.partition"){
    if(length(shapeby) == 0){
      colorlab="Experiment / Class"
    }else{
      colorlab = "Set / Local class"
      shapelab = "Global class"
      treeholes.rar@sam_data$exp.partition = as.factor(str_replace(
        treeholes.rar@sam_data$exp.partition,".Class"," / "))
      Npart = length(levels(treeholes.rar@sam_data$exp.partition))
      usecolor=col_vector
      
    }
  }
}

# ... Set the plot you want
axesby=c(1,3)

# ... Create labels
labfacets=c("Starting communities","Final, Replicate 1","Final, Replicate 2",
            "Final, Replicate 3","Final, Replicate 4")
names(labfacets)=c("Rep0","Rep1","Rep2","Rep3","Rep4")
labaxes=paste0("Axes",paste(axesby[1],axesby[2],sep = "-"))
explainX = round(treeholes.ord$values[axesby[1],2] * 100, digits = 2)
explainY = round(treeholes.ord$values[axesby[2],2]* 100, digits = 2)
xlab = paste0("PCoA, component ", axesby[1]," [",explainX,"%]")
ylab = paste0("PCoA, component ", axesby[2]," [",explainY,"%]")

# ... prepare plot
title=paste(method," of ",dist," distance",sep="")
treeholes.ord$vectors[, "Axis.2"] = -1*treeholes.ord$vectors[, "Axis.2"]
if(length(shapeby) == 0){
  plotTitle=paste("Plot",method,"_",dist,"_Par",sel_partition,
                  "_By",colorby,"_",labelOut,"_",labaxes,".pdf",sep="")
  p = plot_ordination(treeholes.rar, treeholes.ord, color = colorby,
                      axes=axesby) #,label="parent")#,shape="partition")
}else{
  plotTitle=paste("Plot",method,"_",dist,"_Par",sel_partition,
                  "_By",colorby,"-shape",shapeby,"_",labelOut,"_",labaxes,".pdf",sep="")
  treeholes.rar@sam_data$Time0D_7D = as.factor(treeholes.rar@sam_data$Time0D_7D)
  p = plot_ordination(treeholes.rar, treeholes.ord, color = colorby,
                      shape=shapeby,
                      axes=axesby) #,label="parent")#,shape="partition")
}

p = p + geom_point(size = 4.0, alpha = 0.7) + ggtitle(title)
if(length(shapeby) == 0){
    p = p + labs(color = colorlab)
  }else{
    p = p + labs(color = colorlab, shape = shapelab)
}
p = p + xlab(xlab)+ylab(ylab)+
  theme_bw()+
  theme(axis.title = element_text(size=24),
        title=element_blank(), #element_text(size=16),
        axis.text=element_text(size=20),
        strip.text=element_text(size=22),
        legend.title = element_text(size=22),
        legend.text=element_text(size=20))
p = p + facet_wrap(~replicate,nrow=1,ncol=5,
               labeller=labeller(replicate=labfacets))+
  scale_color_manual(values = usecolor) #scale_color_hue()# scale_color_brewer(palette="Accent")


# .... print
pdf(file=plotTitle,width=22,height=6)
print(p)
dev.off()


# ... Plot only time0
#colorby="Location" # "partition" "parent" "Location"
if(colorby=="partition"){
  usecolor=colorCodes
}else{
  usecolor=col_vector
}
axesby=c(1,2)
labaxes="Axes1-2"
plotTitle=paste("Plot",method,"_",dist,"_Par",sel_partition,"_By",colorby,"_",labaxes,"_",labelOut,"_Time0.pdf",sep="")
pdf(file=plotTitle,width=13,height=8)
title=paste(method," of ",dist," distance",sep="")
p = plot_ordination(treeholes.rar.time0, treeholes.ord.time0, color = colorby,
                    axes=axesby) #,label="parent")#,shape="partition")
p = p + geom_point(size = 2.5, alpha = 0.7) + ggtitle(title) +
  theme(axis.title = element_text(size=18),
        title=element_text(size=16),
        axis.text=element_text(size=15),
        strip.text=element_text(size=18),
        legend.title = element_text(size=18),
        legend.text=element_text(size=16))

p + facet_wrap(~replicate,nrow=1,ncol=5)+ scale_color_manual(values = usecolor) #scale_color_hue()# scale_color_brewer(palette="Accent")
dev.off()
