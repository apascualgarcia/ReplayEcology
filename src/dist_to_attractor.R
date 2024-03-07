##########################################
#  dist_to_attractor.R
##########################################
# Here I compute the distance between the starting
# communities and different points of the final 
# communities, trying to characterize the relationship
# between the initial conditions and the final attractors.
# I consider the centroids of the final classes and
# the closest community of the final classes to the
# each starting community (excluding its own) as the references to 
# explain the fate of starting communities.
# Therefore,the representation desired is the distance between the
# starting community and the centroid or closest point of
# the attractor against the maximum number of replicates that
# end up in the same final class.
##########################################
# author = apascualgarcia.github.io
# date = Dec 6th, 2022. ETH-Zürich
##########################################
rm(list=ls())

# Set libraries and functions needed--------
library(ggplot2)
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)
source("extract_subset.R")
source("extract_centroid.R")
source("extract_dist_vec_feat.R")
source("KLD.R")
source("JSD.R")

#source("")
# START EDITING ###########
# --- Set input directories and files
# ..... ASV table
dirOTU = "../6_finalfiles"
fileOTU="seqtable_readyforanalysis.csv"

# ..... Metadata
dirMD = "../7.1_classes"
fileMD="metadata_Time0D-7D-4M_May2022_wJSDpart-merged.csv"

# .... describe here the factors and levels needed to extract each subset
#      you should create lists if needed although current implementation considers
#      vectors. The first factor should be the reference subset, while 
#      subsequent subsets are those against which the reference factor will
#      be compared, finding for them their centroid and closest member.
sample.id = "sampleid" # column for the samples
factor.vec.subsets = c("Experiment", "replicate.partition", "replicate.partition") # one element (or list) for each subset
level.vec.subsets = c("0D", "4.1", "4.2") # the correspondent level

# ... Output directory
dirOut = "../7.5_attractors"
# STOP EDITING #########

# Read files ------------
# --- Read ASVs table
setwd(dirOTU)
ASV.table=read.table(fileOTU,sep="\t")
dim(ASV.table)
ASV.table=as.matrix(ASV.table)
colnames(ASV.table)[1:5]
rownames(ASV.table)[1:5]

# ..... read metadata
setwd(dirMD)
sample_md<-read.table(fileMD,sep="\t",header=TRUE)

# Subset ------------
# --- First remove samples that did not pass the quality control (not in metadata)
matched=match(sample_md$sampleid,row.names(ASV.table)) #attr(data.dist,"Labels"))
sample_md=sample_md[!is.na(matched),]

# --- Remove 4M samples
id.keep = which(sample_md$Experiment != "4M")
sample_md = sample_md[id.keep, ]


# --- Extract subsets of samples
ASV.sub.rel.list = list()

i = 0
for(factor.vec in factor.vec.subsets){
  i = i + 1
  level.vec = level.vec.subsets[i]
  ASV.sub.tmp = extract_subset(ASV.table, sample.md = sample_md,
                               factor.vec = factor.vec, level.vec = level.vec,
                               sample.id = sample.id)  
  
  # ..... Compute relative abundances and introduce pseudocount for JSD 
  ASV.sub.rel.tmp = apply(ASV.sub.tmp, MARGIN = 1,FUN = function(x)(x/sum(x))) # normalize
  ASV.sub.rel.tmp = t(ASV.sub.rel.tmp)
  pseudocount = 1e-07
  ASV.sub.rel.tmp = apply(ASV.sub.rel.tmp, 1:2, 
                       function(x) ifelse (x==0, pseudocount, x))
  ASV.sub.rel.list[[i]] = ASV.sub.rel.tmp
}

# Compute the centroids --------
centroids = list()
for(i in 2:length(factor.vec.subsets)){ # note that we exclude the reference subset
  k=i-1
  centroids[[k]] = extract_centroid(ASV.table, sample.md = sample_md, 
                                    factor.vec = factor.vec.subsets[i], 
                                    level.vec = level.vec.subsets[i],
                                    sample.id = sample.id)
  if(k == 1){
    df.centroids = as.data.frame(t(centroids[[k]]))
  }else{
    df.centroids = rbind(df.centroids, centroids[[k]])
  }
}
names.centroids = level.vec.subsets[2:length(level.vec.subsets)]
names.centroids = paste0("Centroid", names.centroids)
names(centroids) = names.centroids
rownames(df.centroids) = names.centroids
# Compute the distances --------

ASV.ref.relAb = ASV.sub.rel.list[[1]] # reference subset
key_zero = T
nchild = 0
nchild_min = 0
rank_child = vector("numeric", length = dim(ASV.ref.relAb)[1])
for(i in 1:dim(ASV.ref.relAb)[1]){
  sample = rownames(ASV.ref.relAb)[i]
  # ... First compute the distance against the centroids
  dist.vec.tmp = vector("numeric", length = length(centroids))
  for(k in 1:length(centroids)){
    s1 = ASV.ref.relAb[i, ]
    s2 = centroids[[k]]
    s2[s2 == 0] = pseudocount
    dist.vec.tmp[k] = JSD(s1, s2)
  }
  names(dist.vec.tmp) = names(centroids)
  features.list = extract_dist_vec_feat(dist.vec.tmp, "ToC")
  dist.vec = features.list$dist.vec
  char.vec = features.list$char.vec
  dist.vec = c(dist.vec.tmp, dist.vec)
  dist.df.tmp = as.data.frame(t(dist.vec))
  char.df.tmp = as.data.frame(t(char.vec))
  all.df.tmp = cbind(sample, dist.df.tmp, char.df.tmp)
  
  # ... Second, compute the distance against the subsets

  for(k in 2:length(ASV.sub.rel.list)){ # for each subset
    ASV.sub.tmp = ASV.sub.rel.list[[k]]
    dist.vec.tmp = vector("numeric", length = dim(ASV.sub.tmp)[1]) # reinitialize
    key_par = F
    for(u in 1:dim(ASV.sub.tmp)[1]){ # compare "sample" against all the communities in the subset 
      sample_sub = row.names(ASV.sub.tmp)[u]
      id.s.sub = which(sample_md$sampleid == sample_sub)
      sample_par = as.character(sample_md$parent[id.s.sub])
      s1 = ASV.ref.relAb[i, ]
      s2 = ASV.sub.tmp[u, ]
      JSD.tmp = JSD(s1, s2)
      if(JSD.tmp == 0){ # control that there are not weird possibilities
        idtcl.df.tmp = rbind(s1, s2)
        rownames(idtcl.df.tmp) = c(sample, sample_sub)
        if(key_zero == T){
          idtcl.df = idtcl.df.tmp
          key_zero = F
        }else{
          idtcl.df = rbind(idtcl.df, idtcl.df.tmp)
        }
        mes = "identical samples found --> inspect idtcl.df for details"
        warning(mes)
      }
      if(sample == sample_par){ # control that comparison with own daughters is not included
          u_par = u
          name_child = sample_sub
          key_par = T
          #browser()
      } # skip comparison against its own daughters
      dist.vec.tmp[u] = JSD.tmp
    } # for each sample in the subset
    names(dist.vec.tmp) = rownames(ASV.sub.tmp)
    
    # ... check if a child community was the closest
    id.max = which.max(dist.vec.tmp)
    if(key_par == T){ # if a daughter community is found
      nchild = nchild + 1 # count that it exist
      sort_dist = sort(dist.vec.tmp)
      rank_child[i] = which(names(sort_dist) == name_child) 
      if(id.max  == u_par){ # check if it was the one with minimum distance
        nchild_min = nchild_min + 1 # count how many child have the min dist
        dist.vec.tmp[u_par] = NA # remove from further computation
      }
    }
    label = paste0("ToBdC",level.vec.subsets[k])
    features.list = extract_dist_vec_feat(dist.vec.tmp, label = label)
    dist.vec = features.list$dist.vec
    char.vec = features.list$char.vec
    dist.df.tmp = as.data.frame(t(dist.vec))
    char.df.tmp = as.data.frame(t(char.vec))
    all.df.tmp = cbind(all.df.tmp, dist.df.tmp, char.df.tmp)
  } # For each subset

  # ... Look for the closest border
  minDistToBd.vec = grep("minDistToBd",colnames(all.df.tmp))
  minDistToBd = min(all.df.tmp[, minDistToBd.vec])
  id.class.minDistToBd = which.min(all.df.tmp[, minDistToBd.vec])
  id.DminToBd.vec = grep("id.DminToBd",colnames(all.df.tmp))
  id.minDistToBd = all.df.tmp[,id.DminToBd.vec[id.class.minDistToBd]]
  id.class.minDistToBd = paste0("C",id.class.minDistToBd)
  # ... Look for the furthest border among the closest ones
  maxDistToBd = max(all.df.tmp[, minDistToBd.vec])
  id.class.maxDistToBd = which.max(all.df.tmp[, minDistToBd.vec])
  #id.DmaxToBd.vec = grep("id.DmaxToBd",colnames(all.df.tmp))
  id.maxDistToBd = all.df.tmp[,id.DminToBd.vec[id.class.maxDistToBd]]
  id.class.maxDistToBd = paste0("C",id.class.maxDistToBd)
  # ... compute other metrics
  meanDistToBd  = (maxDistToBd + minDistToBd) /2
  diffDistToBd = maxDistToBd - minDistToBd
  diffRelDistToBd = diffDistToBd / minDistToBd
  
  # ... final wrap up
  df.bd = data.frame(maxDistToBd, minDistToBd, meanDistToBd,
                     diffDistToBd, diffRelDistToBd)
  df.bd.id = data.frame(id.minDistToBd, id.maxDistToBd, 
                        id.class.minDistToBd, id.class.maxDistToBd)
  rownames(df.bd.id) = NULL
  all.df.tmp = cbind(all.df.tmp, df.bd, df.bd.id)
  if(i == 1){
    dist.df = all.df.tmp
  }else{
    dist.df = rbind(dist.df, all.df.tmp)
  }
}
# Report the number of times a child sample was the minimum
rank_child = rank_child[which(rank_child != 0)]
nchild_min = length(which(rank_child == 1))
q = quantile(rank_child, probs = c(0,0.05,0.25,0.5,0.75,0.95,1))


#browser()
# Compute the maximum number of replicates --------
# falling within the same class at time 7
Nsame.part = c()
part.id = c()
sample.keep = c()
for(i in 1:dim(ASV.ref.relAb)[1]){
  sample = rownames(ASV.ref.relAb)[i]
  id.sample = which(sample_md$parent == sample)
  if(length(id.sample) != 5){next}else{
    sample.keep = c(sample.keep, i)
    #browser()
  }
  id.keep = which(sample_md$replicate[id.sample] != 0)
  partitions = sample_md$partition[id.sample]
  partitions = partitions[id.keep]
  # ... new version using table to account for new level ids.
  #length.part = table(partitions)
  #N.tmp = max(length.part)
  #part.tmp = names(which.max(length.part))
  # ... commented older version that was using other level ids
   length.part = rle(partitions)
   N.tmp = max(length.part$lengths)
   id.tmp = which.max(length.part$lengths)
   part.tmp = length.part$values[id.tmp]
   Nsame.part = c(Nsame.part, N.tmp)
   part.id = c(part.id, part.tmp)
}
dist.df = dist.df[sample.keep, ]
dist.df$Nsame.class = Nsame.part
dist.df$Class = part.id
id.tie = which(Nsame.part == 2)
dist.df$Class[id.tie] = "none"
dist.df$Converge = as.character(dist.df$Nsame.class)
dist.df$Converge[which(dist.df$Converge == 4)] = 1
dist.df$Converge[which(dist.df$Converge != 1)] = 0
dist.df$Converge = as.factor(dist.df$Converge)

# Inspect ---------------

# ... Extract some df for inspection
# ......dist.df for those having samples split half and half
dist.df.tie = dist.df[id.tie, ]
samples.tie = dist.df$sample[id.tie]
# ......sample_md for them
matched = match(samples.tie, sample_md$sampleid)
sample_md.tie = as.data.frame(sample_md[matched, ])
# ......look for columns of interest
id.id = grep("id.*min.*", colnames(dist.df),perl = T,)
id.id = c(id.id, length(colnames(dist.df)))
dist.df.ids = dist.df[, id.id]
dist.df.tie.ids = dist.df.tie[, id.id]

# Write output --------
setwd(dirOut)
label_ref = level.vec.subsets[1]
label_subsets = c()
for(i in 2:length(level.vec.subsets)){
  label_subsets = paste0(label_subsets,"_",level.vec.subsets[i])
}
  
fileSummary = paste0("Distances_from",
                     label_ref,"_to",label_subsets,".tsv")
fileSummIDs = paste0("Distances-IDs_from",
                     label_ref,"_to",label_subsets,".tsv")
fileCentroid = paste0("ASV_table_centroids",label_subsets,".tsv")
write.table(df.centroids, fileCentroid, sep = "\t", quote = F)
write.table(dist.df, fileSummary, sep = "\t", quote = F, row.names = F)
write.table(dist.df.ids, fileSummIDs, sep = "\t", quote = F, row.names = F)

fileChild = paste0("RankChildCommDist_from",
                   label_ref,"_to",label_subsets,".txt")
# ... divert info childs to file
sink(file = fileCentroid)
mes1 = paste0("The number of times where the minimum distance is found at the child was ", nchild_min,
              "out of ", nchild, " (",nchild_min/nchild*100,"%)")
cat(mes1)
mes2 = paste0("The ranking quantiles are: ")
cat(mes2)
(q)
cat("Representing the top X%")
(q/nchild*100)
sink()
# .... stop sink
# Plot results ---------
#dist.df[] <- lapply(dist.df, function(x) if(is.factor(x)) as.numeric(x)
#               else x)
dist.df$Class = as.factor(dist.df$Class)
class(dist.df$Nsame.class) # integer, convert to factor to control better
dist.df$Nsame.class = as.factor(dist.df$Nsame.class)
# for(i in 1:length(colnames(dist.df))){
#   check.class = class(dist.df[, i])
#   if(check.class == "numeric"){
#     varx = colnames(dist.df)[i]
#     file.plot = paste0("Plot_",varx,"_VS_NsameClass",label_subsets,".pdf")
#     pdf(file.plot, width = 6, height = 6)
#     gg = ggplot(dist.df, aes_string(x = varx,
#                              y = "Nsame.class", 
#                               color = "Class", 
#                              shape = "id.DminToC"))+
#       geom_point() +
#       #geom_jitter()+
#       geom_smooth()
#     print(gg)
#     dev.off()
#   }
# }
# ... Here we plot combinations of variables with common properties
varx.vec = c("meanDistToBd","minDistToBd","minDistToBd","meanDistToBd")
vary.vec = c("meanDistToC","minDistToC","meanDistToC","minDistToC")
xlab.vec = c("Mean","Minimum","Minimum","Mean")
ylab.vec = c("Mean","Minimum","Mean","Minimum")
xlab.label = "distance to borders"
ylab.label = "distance to centroids"
alpha = 0.75
size.Nconverge = FALSE # If true, it will represent three shapes, one for each category 
# of replicates converging. If false, it will represent two categories, = 4 and otherwise
for(i in 1:length(varx.vec)){
  varx = varx.vec[i]
  vary = vary.vec[i]
  xlab = paste(xlab.vec[i],xlab.label)
  ylab = paste(ylab.vec[i],ylab.label)
  file.plot = paste0("Plot_",varx,"_VS_",vary,label_subsets,".pdf")
  pdf(file.plot, width = 7, height = 5)
  gg = ggplot(dist.df, aes_string(x = varx,
                                  y = vary))
  if(size.Nconverge == TRUE){
    gg = gg +
      geom_point(aes_string(color = "Class", 
                            shape = "id.DminToC",
                            size = "Nsame.class"),alpha = alpha) +
      scale_size_manual(values=c("2" = 1.5,"3" = 2.3,"4" = 3.2))
  }else{
    gg = gg +
      geom_point(aes_string(color = "Class", 
                            shape = "id.DminToC",
                            size = "Converge"),alpha = alpha) +
      scale_size_manual(values=c("1" = 3.2,"0" = 1.5),
                        labels=c("1" = "= 4","0" = "< 4"))
  }

  gg = gg + scale_shape_manual(labels = c("1", "2"),
                     values = c(16, 17))+
    scale_color_manual(labels = c("1", "2","both"),
                       values = c("chocolate4","chartreuse","blue"))+ #c("red", "blue", "gold"))+
  xlab(xlab)+ylab(ylab)+
  labs(color = "Converge final class",
      shape = "Closest centroid",
      size = "Num. rep. converge")+
    theme_bw()+
   theme(axis.title = element_text(size = 16),
         axis.text = element_text(size = 12),
         legend.title = element_text(size = 12),
         legend.text = element_text(size = 12))
  #geom_jitter()+
  #geom_smooth()
  print(gg)
  dev.off()
}

