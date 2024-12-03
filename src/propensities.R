###################################
#     propensities.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-17
# Script Name: propensities.R   
# Script Description: This script first classifies communities in different
# trajectories, taking into account if each starting communities has 
# its four replicates ending up in the  same final class (convergent trajectory) or not (divergent).
# Then, for each ASV, it estimates the statistical propensity of being observed
# in each type of trajectory. This estimation is made independently for
# starting and final communities.

rm(list = ls())
# START EDITING ----------------

# --- Input and output files
file.ASV = "seqtable_readyforanalysis.csv"
file.taxa = "taxa_wsp_readyforanalysis.csv"
file.Meta = "metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"
file.RF = "varTop_Class-partition_ntree20000_mtry38.tsv"

# --- Directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirASV=paste(this.dir,"/6_finalfiles/",sep="") # Dir of ASV table
dirMeta=paste(this.dir,"/7.3_phyloseq/",sep="") # Dir of metadata
dirRF=paste(this.dir,"/7.6_keystone/random_forest",sep="") # Dir of output data 
dirOut=paste(this.dir,"/7.6_keystone/propensities",sep="") # Dir of output data 

# --- Packages & scripts

scripts <- c("clean_ASV_table.R",
             "propensities_fun.R","compare_propensities.R",
             "sigboot.R","sigclass.R",
             "heatmap.2.mod.R") # list of functions to load

packages <- c("tidyverse", "stringr", "ggplot2","gplots","UpSetR","venn")#,
              #"ggVennDiagram")# "VennDiagram") # list of packages to load

git.packages <- c("yanlinlin82/ggvenn")

###### STOP EDITING -------------

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")

n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}
# install github packages
new.pkg <- git.packages[!(git.packages %in% installed.packages())] # determine which packages aren't installed

if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
setwd(dirSrc)
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}
#library(ggvenn)

# SOURCE FUNCTIONS ---------

n_scripts <- length(scripts) # count how many packages are required

for(n in 1:n_scripts){
  cat("Loading script #", n, " of ", n_scripts, "... Currently Loading: ", scripts[n], "\n", sep = "")
  lib_load <- paste("source(\"",scripts[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}
# READ INPUT FILES ----------

dir.create(dirOut)

# --- Read ASVs table
setwd(dirASV)
ASV.table.all = read.table(file = file.ASV, sep="\t")
colnames(ASV.table.all)[1:5]
rownames(ASV.table.all)[1:5]
dim(ASV.table.all)
head(ASV.table.all)[1:5,1:5]

# --- Dir taxa
taxa.df = read.table(file = file.taxa, sep = "\t", header = T)

# ..... read metadata. Samples present in metadata were those passing the filtering
setwd(dirMeta)
sample_md.all <-read.table(file = file.Meta, sep="\t", header=TRUE)
head(sample_md.all)[1:5,1:5]

# ..... read top ASVs inferred from RF
setwd(dirRF)
top.vars = read.table(file = file.RF)

# Clean data   ----------
# --- Extract a set with only time 7
exclude_exp = c("4M","0D")
clean.data.list = clean_ASV_table(ASV.table.all, sample_md.all,
                                  match_exp = T, exclude_exp = exclude_exp)

ASV.table.7D = clean.data.list$ASV.table
sample_md.7D = clean.data.list$sample_md

# --- Another set with only time 0
exclude_exp = c("4M","7D_rep1","7D_rep2","7D_rep3","7D_rep4")

clean.data.list = clean_ASV_table(ASV.table.all, sample_md.all,
                                  match_exp = T, exclude_exp = exclude_exp)

ASV.table.0D = clean.data.list$ASV.table
sample_md.0D = clean.data.list$sample_md

# ... create frequencies tables
ASV.table.0D.freq = ASV.table.0D/rowSums(ASV.table.0D)
ASV.table.7D.freq = ASV.table.7D/rowSums(ASV.table.7D)
stages = c("0D","7D")
ASV.table.freq.list = list()
ASV.table.freq.list[["0D"]] = ASV.table.0D.freq
ASV.table.freq.list[["7D"]] = ASV.table.7D.freq

# COMPUTE PROPENSITIES --------

# --- Extract parents
parents.list = unique(sample_md.7D$parent)
#parent = parents.list[1] # debug

# --- Compute the probability that each parent has of having
#     their childs in a pure (all in class 1 or 2) or mixed classification.
child_to_parents = vector(mode = "character", length = dim(ASV.table.7D)[1])
names(child_to_parents) = rownames(ASV.table.7D)
n_class1 = 0; class1_id = c(); parent1_id = c()
n_class2 = 0; class2_id = c(); parent2_id = c()
n_classmix = 0; classmix_id = c(); parentmix_id = c()
for(parent in parents.list){
  #browser()
  child_id = which(sample_md.7D$parent == parent) # extract childs
  child_classes = as.character(sample_md.7D$partition[child_id]) # and their classifications
  nclass = unique(child_classes) # check if unique
  child_to_parents[child_id]=parent
  if(length(nclass) == 1){ # pure class
    if(nclass[1] == "Class1"){ # either class 1
          n_class1 = n_class1 + 1
          class1_id = c(class1_id, child_id)
          parent1_id = c(parent1_id, parent)
    }else{                    # or class 2
      n_class2 = n_class2 + 1
      class2_id = c(class2_id, child_id)
      parent2_id = c(parent2_id, parent)
    }
  }else{  # mixed classes
    n_classmix = n_classmix + 1
    classmix_id = c(classmix_id, child_id)
    parentmix_id = c(parentmix_id, parent)
  }
}
n_class_tot = n_class1 + n_class2 + n_classmix

P_class1 = n_class1 / n_class_tot
P_class2 = n_class2 / n_class_tot
P_classmix = n_classmix / n_class_tot
P_class.list = list("class1" = P_class1,
                    "class2" = P_class2,
                    "classmix" = P_classmix)

# START ANALYSIS ---------------------
experiment = c("starting","final")
P_summary_all.list = list()
ASV_sig_prop.list = list()
for(exp in experiment){
  
  if(exp == "starting"){
    ASV.table = ASV.table.0D
    sample_md = sample_md.0D
    classes_id.list = list("class1" = parent1_id, 
                           "class2" = parent2_id, 
                           "classmix" = parentmix_id)
  }else{
    ASV.table = ASV.table.7D
    sample_md = sample_md.7D
    classes_id.list = list("class1" = class1_id, 
                           "class2" = class2_id, 
                           "classmix" = classmix_id)
  }
  # --- Extract ASV tables for each class and binarize
  
  class_names = names(classes_id.list)
  ASV.table.class.bin.list = list()
  for(class in class_names){
    class_id = classes_id.list[[class]]
    # ... identify samples for each type of class/mix
    if(exp == "starting"){
      samples_class = class_id
    }else{
      samples_class = sample_md$sampleid[class_id]
    }
    # ... extract the ASV table
    matched = match(samples_class, rownames(ASV.table))
    ASV.table.class = ASV.table[matched, ]
    # ... binarize
    ASV.table.class.bin = ASV.table.class
    ASV.table.class.bin[ASV.table.class.bin > 0] = 1
    ASV.table.class.bin.list[[class]]=ASV.table.class.bin
  }
  ASV.table.bin = ASV.table
  ASV.table.bin[ASV.table.bin > 0] = 1
  ASV.table.class.bin.list[["all"]] = ASV.table.bin
  
  class_names_ext = c(class_names, "all")
  
  # MAIN CALL ------------
  # --- Compute observed propensities
  Prop_freq = propensities_fun(ASV.table.class.bin.list, P_class.list,
                               class_names)
  Propensities = Prop_freq$Prop.list
  freq_ASV.list = Prop_freq$freq_ASV.list
  test = unlist(lapply(freq_ASV.list,function(x){x["ASV_49"]}))
  log(test[2:4])-log(test[1])-log(unlist(P_class.list))
  proptest = unlist(lapply(Propensities,function(x){x["ASV_49"]}))
  # --- Compute bootstrapped propensities
  nboot = 1000
  cat("--- Bootstrapping propensities")
  for(k in 0:nboot){
    ASV.boot.list = list()
    for(class in class_names_ext){
      ASV.boot = ASV.table.class.bin.list[[class]]
      dim_boot = dim(ASV.boot)[1]
      id_boot = sample.int(dim_boot, size = dim_boot, replace = T)
      ASV.boot = ASV.boot[id_boot, ]
      ASV.boot.list[[class]] = ASV.boot
      # extract ASV table
      # randomize
      # call propensities_fun
      # store, etc.
    }
    Prop_freq.boot = propensities_fun(ASV.boot.list, P_class.list,
                                      class_names)
    Propensities.boot = Prop_freq.boot$Prop.list
    if(k == 0){
      P_obs_class1 = Propensities[[class_names[1]]]
      P_obs_class2 = Propensities[[class_names[2]]]
      P_obs_classmix = Propensities[[class_names[3]]]
      P_boot_class1 = as.data.frame(P_obs_class1) # we also include observed value in the first line
      P_boot_class2 = as.data.frame(P_obs_class2) # of bootstrapped dfs
      P_boot_classmix = as.data.frame(P_obs_classmix)
    }else{
      P_boot_class1 = cbind(P_boot_class1, as.data.frame(Propensities.boot[[class_names[1]]]))
      P_boot_class2 = cbind(P_boot_class2, as.data.frame(Propensities.boot[[class_names[2]]]))
      P_boot_classmix = cbind(P_boot_classmix, as.data.frame(Propensities.boot[[class_names[3]]]))
    }
  }
  
  #P_boot_class1[P_boot_class1 == NaN] = 0 # aware of this
  # .... rename columns
  bootnames = paste0("boot_",seq(1:nboot))
  colnames(P_boot_class1)[2:(nboot+1)] = bootnames
  colnames(P_boot_class2)[2:(nboot+1)] = bootnames
  colnames(P_boot_classmix)[2:(nboot+1)] = bootnames
  
  # EXTRACT STATISTICS ---------
  
  # --- Bootstrapped propensities for each ASV (alpha = 0.05)
  q=c(0.05,0.95)
  qntl = apply(P_boot_class1[,2:(nboot+1)], 
               MARGIN = 1, FUN = function(x){quantile(x, na.rm = T, p = q)})
  P_boot_summary_class1 = cbind(P_boot_class1[,1], t(qntl))
  sig = apply(P_boot_summary_class1, MARGIN = 1,
              FUN = function(x){sigboot(x[1],x[2],x[3])})
  P_boot_summary_class1 = cbind(P_boot_summary_class1, sig)
  colnames(P_boot_summary_class1)[1] = c("P_obs")
  colnames(P_boot_summary_class1) = paste0(colnames(P_boot_summary_class1),
                                           "_class1")
  dim(P_boot_summary_class1)
  ##
  qntl = apply(P_boot_class2[,2:(nboot+1)], 
               MARGIN = 1, FUN = function(x){quantile(x, na.rm = T, p = q)})
  P_boot_summary_class2 = cbind(P_boot_class2[,1], t(qntl))
  sig = apply(P_boot_summary_class2, MARGIN = 1,
              FUN = function(x){sigboot(x[1],x[2],x[3])})
  P_boot_summary_class2 = cbind(P_boot_summary_class2, sig)
  colnames(P_boot_summary_class2)[1] = c("P_obs")
  colnames(P_boot_summary_class2) = paste0(colnames(P_boot_summary_class2),
                                           "_class2")
  dim(P_boot_summary_class2)
  ##
  qntl = apply(P_boot_classmix[,2:(nboot+1)], 
               MARGIN = 1, FUN = function(x){quantile(x, na.rm = T, p = q)})
  P_boot_summary_classmix = cbind(P_boot_classmix[,1], t(qntl))
  sig = apply(P_boot_summary_classmix, MARGIN = 1,
              FUN = function(x){sigboot(x[1],x[2],x[3])})
  P_boot_summary_classmix = cbind(P_boot_summary_classmix, sig)
  colnames(P_boot_summary_classmix)[1] = c("P_obs")
  colnames(P_boot_summary_classmix) = paste0(colnames(P_boot_summary_classmix),
                                        "_classmix")
  dim(P_boot_summary_classmix)
  
  P_summary_all = as.data.frame(cbind(P_boot_summary_class1,P_boot_summary_class2,P_boot_summary_classmix))
  
  # COMPARE PROPENSITIES AMONG CLASSES ---------
  # --- Identify significant propensities for each class
  
  columns_propensities = list()
  columns_propensities[["class1"]] = c(1,2,3)
  columns_propensities[["class2"]] = c(5,6,7)
  columns_propensities[["classmix"]] = c(9,10,11)
  Prop.compare = compare_propensities(P_summary_all, columns_propensities)
  
  P_summary_all = Prop.compare$P_summary
  ASV.sig = Prop.compare$ASV_sig
  #stop()
  var = paste(exp,"signif","_")
  ASV_sig_prop.list[[exp]][var]
  P_summary_all.list[[exp]] = P_summary_all
  
  # --- Print some results and save
  
  setwd(dirOut)
  fileOut = paste0("Propensities_",exp,"-ASVs_forFinalClass.tsv")
  write.table(P_summary_all, file = fileOut, quote = F, sep = "\t")
  
  # SUMMARY SIGNIFICANT PROPENSITIES ---------
  for(class in class_names){
    var = paste("sig",class,sep="_")
    id.pos.sig = which(P_summary_all[, var] == 1)
    id.neg.sig = which(P_summary_all[, var] == -1)
    ASV.pos.sig = rownames(P_summary_all)[id.pos.sig]
    ASV.neg.sig = rownames(P_summary_all)[id.neg.sig]
    var = paste(exp,class,sep = "_")
    ASV_sig_prop.list[[var]] = ASV.pos.sig
  }
}
#stop()
# VENN DIAGRAMS ------------
# --- Plot Venn diagram
#
setwd(dirOut)

str(ASV_sig_prop.list)
reorder.vec = c(1,3,2,5,6,4)
reorder.names = names(ASV_sig_prop.list)[reorder.vec]
ASV_venn.list = ASV_sig_prop.list[reorder.names]
names(ASV_venn.list) = c("converge class 1 (start)","diverge (start)",
           "converge class 2 (start)","converge class 2 (end)",
           "diverge (end)","converge class 1 (end)")
plotOut = "Plot_VennDiagram_all.pdf"
pdf(file = plotOut,width = 9, height = 9)
venn::venn(ASV_venn.list,
           ilabels = "counts", 
           ilcs = 1, sncs = 1.2,
           zcolor = "style",
           box = F,
           opacity = 0.4)

dev.off()


# plotOut = "Plot_VennDiagram_startingVsfinal.pdf"
# pdf(file = plotOut,width = 9, height = 9)
# ggvenn(
#   ASV_sig_prop.list,
#   columns = c("starting_class1","starting_class2",
#               "final_class1","final_class2"),
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4,
#   show_outside = "auto"
# )
# dev.off()
# 
# plotOut = "Plot_VennDiagram_startingVslost.pdf"
# pdf(file = plotOut,width = 9, height = 9)
# ggvenn(
#   ASV_sig_prop.list,
#   columns = c("starting_class1","starting_class2","starting_classmix",
#               "lost"),
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4,
#   show_outside = "auto"
# )
# dev.off()
# 
# plotOut = "Plot_VennDiagram_mixVsfinal.pdf"
# pdf(file = plotOut,width = 9, height = 9)
# ggvenn(
#   ASV_sig_prop.list,
#   columns = c("final_class1","final_class2",
#               "starting_classmix",
#               "final_classmix"),
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4,
#   show_outside = "auto"
# )
# dev.off()

# TEST CORE SPECIES HYPOTHESIS ---------

# --- In Venn diagrams we observed a set of ASVs with propensity for class1 (st and end)
# ..... extract all intersections
inters_list = venn(ASV_sig_prop.list,show.plot = F)
fileOut = "Intersections_PropensitySets.RDS"
saveRDS(inters_list,file = fileOut)

# --- Create plots for all intersections:
setwd(dirSrc)
source("propensities_intersections.R")

# --- Explore a bit more specific subsets
# .... select  subsets
select = c("starting_class1:final_class1",
           "starting_class1", "final_class1",
           #"starting_class2:final_class2",
           "starting_class2", "final_class2",
           "starting_classmix","final_classmix",
           #"starting_class1:final_class1",
           "starting_class1:starting_classmix:final_class1:final_class2:final_classmix")
# .... choose new labels for plots below
new_names = c( "class 1 (starting and final)",
               "class 1 (starting only)", "class 1 (final only)",
               #"class 2 (starting and final)",
               "class 2 (starting only)", "class 2 (final only)",
               "divergent (starting only)", "divergent (final only)",
               "cosmopolitan")

all_inters = names(attributes(inters_list)$intersections) # all subsets

prob_q = c(0,0.05,0.25, 0.5,0.75,0.95,1) # quantiles inspected
fileOut = "Stats_ASVpropensitySetsFreqs_byTimeAndTrajectory.txt"
sink(file = fileOut)
cat("# Quantiles of the ASVs frequencies sums of each ASV group, defined by their propensities \n ")
cat("# ... across all communities belonging to specific time point  (starting/final) and trajectory. \n ")

for(set in all_inters){
 cat(paste0("*** \n"))
 cat(paste0("*** Extracting mean frequencies quantiles for set: >>>",set,"<<< \n"))
 cat(paste0("*** \n"))
 for(stage in stages){
   if(stage == "0D"){
     classes_id.list = list("class1" = parent1_id, 
                            "class2" = parent2_id, 
                            "classmix" = parentmix_id)
   }else{
     classes_id.list = list("class1" = class1_id, 
                            "class2" = class2_id, 
                            "classmix" = classmix_id)
   }
   ASV.table.freq.tmp = ASV.table.freq.list[[stage]]
   cat(paste0("==== in ASV table: ",stage,"\n"))
   # ..... identify the one between starting and final class 1
   inters_class1 = unlist(attributes(inters_list)$intersections[set]) #$`starting_class1:final_class1`
   
   # ..... match the labels in the ASV table
   matched = match(colnames(ASV.table.freq.tmp), inters_class1)
   N_ASVs = length(which(!is.na(matched) == T))
   if(N_ASVs < 2){
     cat(paste0("*** WARNING: Less than two elements found, skipping... \n"))
     next
   }else{
     cat(paste0("====> Number of ASVs is: ",N_ASVs, "\n"))
   }
   freq_core_in_class1 = rowSums(ASV.table.freq.tmp[classes_id.list[["class1"]],
                                                    !is.na(matched)])
   freq_core_in_class2 = rowSums(ASV.table.freq.tmp[classes_id.list[["class2"]],
                                                    !is.na(matched)])
   freq_core_in_classmix = rowSums(ASV.table.freq.tmp[classes_id.list[["classmix"]],
                                                      !is.na(matched)])
   #mean_core_in_class1 = mean(freq_core_in_class1)
   qq = quantile(freq_core_in_class1)
   cat(paste0("....... communities converging to class 1: \n"))
   print(qq)   
   #0%         25%         50%         75%        100% 
   #0.001721253 0.196251300 0.410655010 0.602931633 0.966434855 
   #mean_core_in_class2 = mean(freq_core_in_class2)
   qq = quantile(freq_core_in_class2, probs = prob_q)
   cat(paste0("....... communities converging to class 2: \n"))
   print(qq)
   #0%         25%         50%         75%        100% 
   #0.001136283 0.004529859 0.011546921 0.037273554 0.933253795  
   #quantile(freq_core_in_class2, probs = prob_q) # highly skewed, explore further
   #80%        85%        90%        95%       100% 
   #0.04959852 0.19113577 0.28976649 0.59320835 0.93117695 
   #mean_core_in_classmix = mean(freq_core_in_classmix)
   qq = quantile(freq_core_in_classmix, probs = prob_q)
   cat(paste0("....... communities diverging: \n"))
   print(qq)
   #25%        50%        75%        95%       100% 
   #0.01084323 0.05912757 0.39224138 0.71098713 0.85598910 
   
   #*** This part is no longer used. If new columns are added to include the
   # experiment and stage, you can get a df useful for several plots and computations
   # freq_core_in_class1 = data.frame(freq_core_in_class1,
   #                                  rep("convergent, class 1",
   #                                      times = length(freq_core_in_class1)))
   # colnames(freq_core_in_class1) = c("freq","Trajectory")
   # freq_core_in_class2 = data.frame(freq_core_in_class2,
   #                                  rep("convergent, class 2",
   #                                      times = length(freq_core_in_class2)))
   # colnames(freq_core_in_class2) = c("freq","Trajectory")
   # freq_core_in_classmix = data.frame(freq_core_in_classmix,
   #                                    rep("divergent",
   #                                        times = length(freq_core_in_classmix)))
   # colnames(freq_core_in_classmix) = c("freq","Trajectory")
   # 
   # # ... how many communities have freqs lower than 0.1? 
   # length(which(freq_core_in_class1 < 0.1)) # 31
   # length(which(freq_core_in_class2 < 0.1)) # 36
   # length(which(freq_core_in_classmix < 0.1)) # 36  
 }
}

sink()

#**** The following should be fixed after recent changes in the code
#* (a more general version coded in the function propensities_intersections.R)
# freq_core.df = rbind(freq_core_in_class1,freq_core_in_class2,freq_core_in_classmix)
# freq_core.df = freq_core.df[order(-freq_core.df$freq), ]
# freq_core.df$rank = seq(1,dim(freq_core.df)[1])
# 
# # ...... Plot distribution of frequencies
# plotOut = "Plot_Frequency_CoreClass1Starting_vs_trajectory.pdf"
# pdf(plotOut, width = 10)
# gg = ggplot(freq_core.df,
#             aes(x = rank, y = freq, color = Trajectory))+
#   geom_point(alpha = 0.3, size = 4) + #, 
#   #           position=position_jitter(height=.3, width=.3))+
#   geom_hline(yintercept = 0.1)+
#   #geom_line()+
#   #scale_y_log10()+
#   ylab("Sum of ASVs relative abundances")+
#   ggtitle("Starting communities, core ASVs in final class 1")+
#   theme_bw() +
#   theme(title = element_text(size = 18),
#         axis.title = element_text(size = 22),
#         axis.text =  element_text(size = 18),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 22),
#         legend.position = "inside",
#         legend.position.inside = c(0.8,0.8))
# 
# gg
# dev.off()

# IDENTIFY ASVs ----------
#stop()

Ninters = length(attributes(inters_list)$intersections)
Nsel = 25
taxa_prop.list = list()
all_taxa = c()
for(i in 1:Ninters){
  taxa.tmp = unlist(attributes(inters_list)$intersections[i])
  set.tmp = names(attributes(inters_list)$intersections[i])
  matched = match(taxa.tmp, taxa.df$ASV_names)
  taxa_info = as.character(taxa.df[matched,"Family"])
  taxa_prop.list[[set.tmp]] = sort(table(taxa_info),decreasing = T)
  all_taxa = c(all_taxa,taxa_prop.list[[set.tmp]][1:Nsel])
}

all_taxa = unique(names(all_taxa))
prop.matrix = matrix(0, nrow = Ninters, ncol = length(all_taxa))
colnames(prop.matrix) = all_taxa
rownames(prop.matrix) = names(taxa_prop.list)

for(set in names(taxa_prop.list)){
  taxa_tmp = taxa_prop.list[[set]]
  Ntax = min(Nsel,  length(taxa_tmp))
  taxa_tmp = taxa_prop.list[[set]][1:Ntax]
  prop.matrix[set,names(taxa_tmp)] = taxa_tmp
}
#stop()

prop.matrix = prop.matrix[select, ]
prop.matrix = prop.matrix[, colSums(prop.matrix) > 0]
rownames(prop.matrix) = new_names

library("RColorBrewer")
my_palette <- colorRampPalette(c("darkgray", "yellow", "red"))(n = (max(prop.matrix)+1))
plotOut = "heatmap_numASVsXfamilyVsTrajectoryPropensity.pdf"
pdf(file = plotOut, width = 12, height = 15)
heatmap.2.mod(t(prop.matrix),
          dendrogram = "none",scale = "none",
          trace = "none",density.info = "none", 
          Colv = F,
          cex.lab = 3, cexCol = 2, cexRow =1,
          offsetRow = 0, srtCol = 50,
          xlab = "Propensity for trajectory (experiment)", ylab = "Family",
          key = T,key.xlab = "Number of ASVs",
          key.par = list(cex.lab = 2),
          keysize = 0.7, key.title = "",
          margins = c(25, 20),
          col = my_palette)
dev.off()

# RULE OUT OTHER HYPOTHESIS  --------
# --- Identify ASVs that were not resurrected
ASV_match.id = match(colnames(ASV.table.0D),colnames(ASV.table))
ASV_lost = which(is.na(ASV_match.id))
ASV_present = which(!is.na(ASV_match.id))
ASV_lost.id = colnames(ASV.table.0D)[ASV_lost]
ASV_present.id = colnames(ASV.table.0D)[ASV_present]


ASV_sig_prop.list[["present"]] = ASV_present.id
ASV_sig_prop.list[["lost"]] = ASV_lost.id


rowSums(ASV.table.0D.freq)[1:10]
ASV_lost.freq.mean = colMeans(ASV.table.0D.freq)[ASV_lost]
ASV.freq.0D.mean = colMeans(ASV.table.0D.freq)
ASV.freq.7D.mean = colMeans(ASV.table.7D.freq)
ASV.freq.0D.mean.df = data_frame(ASV.freq.0D.mean)
ASV.freq.0D.mean.df$resurrected = "Yes"
ASV.freq.0D.mean.df$resurrected[ASV_lost] = "No"
rownames(ASV.freq.0D.mean.df) = names(ASV.freq.0D.mean)
ASV.freq.0D.mean.df = 
  ASV.freq.0D.mean.df[order(-ASV.freq.0D.mean.df$ASV.freq.0D.mean), ]
ASV.freq.0D.mean.df$rank = seq(1,dim(ASV.freq.0D.mean.df)[1])

# ..... extract some statistics
# ....... quantiles of the frequencies
quantile(ASV.freq.0D.mean)
quantile(ASV_lost.freq.mean)
#0%          25%          50%          75%         100% 
#1.103712e-07 8.601029e-06 4.581363e-05 2.003135e-04 1.011489e-01 # all
#1.153473e-07 1.628884e-05 8.895504e-05 2.000281e-04 1.974997e-03 # lost

# ...... proportion of reads lost
total_reads = sum(ASV.table.0D)
total_reads_present = sum(ASV.table.0D[,ASV_present])
total_reads_lost = sum(ASV.table.0D[,ASV_lost])
(prop_lost = total_reads_lost / total_reads) # 3.50%
(prop_present = total_reads_present / total_reads) # 96.49%

#### RESULTS
# 26% ASVs were not resurrected
# highest mean abundance of lost is 1.97e.-03
# median mean abundance of lost is 9.89e-05
# ASVs lost contributed to 3.45% of the reads observed in
# starting communities
######## 



# ...... Plot distribution of frequencies
plotOut = "Plot_Frequency_ASV0D_ObservedVsLost.pdf"
pdf(plotOut)
gg = ggplot(ASV.freq.0D.mean.df,
            aes(x = rank, y = ASV.freq.0D.mean, color = resurrected))+
  geom_point(alpha = 0.3, 
             position=position_jitter(height=.3, width=.3))+
  geom_line()+
  scale_y_log10()+
  ylab("ASV mean frequency")+
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))
  
gg
dev.off()

# ..... Match lost with Random Forest results
match_lost_RF = match(top.vars$x,ASV_lost.id) # none found
match_present_RF = match(top.vars$x,ASV_present.id) # all found

matched = match(top.vars$x, names(ASV.freq.0D.mean))
freq_top_RF_0D = ASV.freq.0D.mean[matched]
matched = match(top.vars$x, names(ASV.freq.7D.mean))
freq_top_RF_7D = ASV.freq.7D.mean[matched]

plot(freq_top_RF_0D, freq_top_RF_7D)

#### RESULT
# None of the ASVs top-ranked by RF were lost
# Some of the ASVs top-ranked have very low frequencies
# --- 0D
#ASV_61       ASV_16       ASV_25       ASV_97      ASV_555 
#0.0064836595 0.0093587634 0.0201178332 0.0012974220 0.0002438752 
#ASV_10      ASV_326      ASV_128        ASV_7       ASV_34 
#0.0115788597 0.0000811271 0.0043386977 0.0324894101 0.0055128769 
#ASV_14       ASV_45        ASV_8        ASV_2        ASV_6 
#0.0033413257 0.0077317591 0.0622384628 0.0443696927 0.0218279558 
# --- 7D
#ASV_61       ASV_16       ASV_25       ASV_97      ASV_555 
#2.872850e-03 2.030276e-02 6.087088e-03 6.547903e-04 5.824622e-06 
#ASV_10      ASV_326      ASV_128        ASV_7       ASV_34 
#2.949281e-02 2.010420e-06 2.217810e-05 3.154162e-02 4.698859e-03 
#ASV_14       ASV_45        ASV_8        ASV_2        ASV_6 
#2.020553e-02 3.100777e-03 1.630755e-02 5.081778e-02 3.694864e-02 
########

# ...... Match top-ranked RF with high propensities for split


for(class in names(ASV_sig_prop.list)){
  if(class == "lost"){next}
  ASV_class = ASV_sig_prop.list[[class]]
  freq_highprop_0D = ASV.freq.0D.mean[ASV_class]
  cat("** The sum of mean frequencies of ",class," members is:",sum(freq_highprop_0D))
  cat("** The intersection with non-resurrected ASVs \n")
  matched = match(ASV_class, ASV_lost.id)
  ASV_inters =  ASV_class[!is.na(matched)]
  freq_inters = ASV.freq.0D.mean[ASV_inters]
  Nint = length(ASV_inters)
  ratio_inters = sum(freq_inters) / sum(freq_highprop_0D)
  cat("** has ",Nint," members with a sum of mean frequencies",sum(freq_inters))
  cat("  representing a fraction of the class: ",ratio_inters,"\n")
}
#** The sum of mean frequencies of  starting_class1  members is: 0.9215132** The intersection with non-resurrected ASVs 
#** has  214  members with a sum of mean frequencies 0.03910438  representing a fraction of the class:  0.04243496 
#** The sum of mean frequencies of  starting_class2  members is: 0.3030728** The intersection with non-resurrected ASVs 
#** has  28  members with a sum of mean frequencies 0.004164248  representing a fraction of the class:  0.01374009 
#** The sum of mean frequencies of  starting_classmix  members is: 0.2929441** The intersection with non-resurrected ASVs 
#** has  35  members with a sum of mean frequencies 0.007297948  representing a fraction of the class:  0.02491242

ASV_sig_prop_classmix = ASV.sig.list[[class_names[3]]]
matched = match(top.vars$x, ASV_sig_prop_classmix)
ASV_sig_prop_classmix[!is.na(matched)]

matched = match(ASV_sig_prop_classmix, names(ASV.freq.0D.mean))
freq_highprop_0D = ASV.freq.0D.mean[matched]
freq_highprop_0D[is.na(freq_highprop_0D)] = 0
matched = match(ASV_sig_prop_classmix,  names(ASV.freq.7D.mean))
freq_highprop_7D = ASV.freq.7D.mean[matched]
(freq_highprop_7D - freq_highprop_0D)
(freq_highprop_0D)
plot(freq_highprop_0D, freq_highprop_7D)

#### RESULT
# None of the ASVs top-ranked by RF had a high propensity
# of being found in the mixed class
# ASVs having a high propensity had the following starting freqs
#  ASV_75       ASV_83       ASV_96      ASV_130      ASV_212 
# 2.271083e-03 2.502503e-03 1.181271e-03 2.477562e-04 2.268792e-04 
# ASV_223         ASV_810     ASV_1445 
# 4.594386e-04        NA      1.180231e-05 
#
# Interestingly, ASV_810 was not detected
#
# Have a look at the top RF ranked

matched = match(top.vars$x, rownames(P_summary_all))
View(P_summary_all[matched, ])
# top ranked ASV tend to have a significantly positive propensity
# for being observed in final class 1





################# LEFTOVERS:
# First dirty version of propensities comparison, some relevant summaries
# to keep before I further modify the code
# # --- Individual positive significant propensities (CI does not include 0)
# id.pos.sig.class1 = which(P_summary_all$sig_class1 == 1) # 848
# id.neg.sig.class1 = which(P_summary_all$sig_class1 == -1) # none
# #848 significant propensity for being observed in class 1
# id.pos.sig.class2 = which(P_summary_all$sig_class2 == 1) # 396
# id.neg.sig.class2 = which(P_summary_all$sig_class2 == -1) # 34
# #396 positive and 34 negative significant propensity for being observed in class 2
# id.pos.sig.classmix = which(P_summary_all$sig_classmix == 1) # 648
# id.neg.sig.classmix = which(P_summary_all$sig_classmix == -1) # 1
# #648 positive and 1 negative significant propensity for being observed in class mix
# # --- Individual positive significant propensities in more than one class
# id.sig.class1_and_2 = c(id.pos.sig.class1,id.pos.sig.class2)
# freq.class1_and_2 = as.data.frame(table(id.sig.class1_and_2))
# id.tmp = which(freq.class1_and_2$Freq == 2)
# id.sig.class1_and_2 = 
#   as.numeric(freq.class1_and_2$id.sig.class1_and_2[id.tmp])
# # 321 significant (positive) propensity in both 1 and 2
# id.sig.class1_and_classmix = c(id.pos.sig.class1,id.pos.sig.classmix)
# freq.class1_and_classmix = as.data.frame(table(id.sig.class1_and_classmix))
# id.tmp = which(freq.class1_and_classmix$Freq == 2)
# id.sig.class1_and_classmix = 
#   as.numeric(freq.class1_and_classmix$id.sig.class1_and_classmix[id.tmp])
# # 540 siginificant propensity in both 1 and mix
# id.sig.class2_and_classmix = c(id.pos.sig.class2,id.pos.sig.classmix)
# freq.class2_and_classmix = as.data.frame(table(id.sig.class2_and_classmix))
# id.tmp = which(freq.class2_and_classmix$Freq == 2)
# id.sig.class2_and_classmix = 
#   as.numeric(freq.class2_and_classmix$id.sig.class2_and_classmix[id.tmp])
# # 333 significant propensity in 2 and mix
# id.sig.class_all =  c(id.pos.sig.class1,id.pos.sig.class2,id.pos.sig.classmix)
# freq.class_all = as.data.frame(table(id.sig.class_all))
# id.tmp = which(freq.class_all$Freq == 3)
# id.sig.class_all = 
#   as.numeric(freq.class_all$id.sig.class_all[id.tmp])
# # 295 significant propensity in 3