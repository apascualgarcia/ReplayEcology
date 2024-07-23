 ###################################
## communities_superposition.R
###################################
# Given two sets of starting and final
# communities encoded in matrices P and Q0 of
# dimensions n x d, with n the number of communities
# and d the number of taxa, we look for a matrix U
# finding the optimal rotation between both
# matrices. Once we recover the matrix
# U we test the prediction by taking
# random samples of P: Pr1, Pr2, ... and transforming
# them with U to get Qr1, Qr2, ..., and we plot results 
# against one or more (new) sets of
# Q matrices, Q1, Q2, ...
###################################
## Author: Dr. Alberto Pascual-García
## Date Created: 2022-12-13
## MIT Licensed.
## Contact: apascualgarcia.github.io
#################################### 
##
## INPUT:
#  * An ASV vs samples matrix containing all communities.
#  * A metadata file to split the ASV matrix in subsets P, Q1, Q2, etc.
## OUTPUT:  
#  * The transformation matrix U and the RMSD between P and Q0
#  * The mean and standard deviation of Pr %*% U
## DEPENDENCIES:
# "extract_subset.R", "kabsch.R","aitchinson.R",
# "KLD.R","JSD.R","reshape_df_to_ggplot.R" and functions therein
#
#####################################
## Source libraries if needed
rm(list=ls())

## START EDITING ###########

# Set libraries and functions needed--------
library(ggplot2)
library(boot) # bootstrap rmsd
library(gtools) # draw dirichlet distributions
library(reshape2)
library(phyloseq)
library(ggpubr)

function.vec = c("extract_subset.R", "kabsch.R","aitchinson.R",
                 "KLD.R","JSD.R","reshape_df_to_ggplot.R" ) 

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
factor.vec.subsets = as.character(rep("Experiment", 5)) # one element (or list) for each subset
level.vec.subsets = c("0D", "7D_rep1", "7D_rep2", "7D_rep3", "7D_rep4") # the correspondent level
Rep = 4 # Replicate used for the transformation

# .... Compute confidence intervals for rmsd and cross-covariance
compute_CI = T # it takes several minutes, computed only once

# ... Output directory
dirOut = "../9_predictions"
## STOP EDITING ###########

# Source functions --------
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)
for(i in 1:length(function.vec)){
  function.name = function.vec[i]
  mes = paste("** Sourcing script", function.name,"\n")
  cat(mes)
  source(function.name)
}

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

# --- Extract subsets of samples
ASV.sub.rel.list = list()
ASV.sub.list = list()
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
  rownames(ASV.sub.rel.tmp) = rownames(ASV.sub.tmp)
  colnames(ASV.sub.rel.tmp) = colnames(ASV.sub.tmp)
  
  ASV.sub.rel.list[[i]] = ASV.sub.rel.tmp
  ASV.sub.list[[i]] = ASV.sub.tmp
  
}
ASV.ref.relAb = ASV.sub.rel.list[[1]]
ASV.ref.sub = ASV.sub.list[[1]]

lapply(ASV.sub.list, FUN = dim)

# --- Work only with samples that were resurrected
sample.keep = c()
for(i in 1:dim(ASV.ref.relAb)[1]){
  sample = rownames(ASV.ref.relAb)[i]
  id.sample = which(sample_md$parent == sample)
  if(length(id.sample) != 5){next}else{
    sample.keep = c(sample.keep, i)
    #browser()
  }
}

ASV.ref.relAb = ASV.ref.relAb[sample.keep, ]
ASV.ref.sub = ASV.ref.sub[sample.keep, ]
dim(ASV.ref.relAb)

# --- Ensure that matrices samples are in the same order:
for(i in 2:length(ASV.sub.rel.list)){
  vec2match = sapply(rownames(ASV.sub.rel.list[[i]]), FUN = function(x){
    unlist(strsplit(x, split = ".", fixed = TRUE))[1]})
  matched = match(rownames(ASV.ref.relAb), vec2match)
  test = which(is.na(matched))
  if(length(test) > 0){stop()}
  ASV.sub.rel.list[[i]] = ASV.sub.rel.list[[i]][matched, ]
  ASV.sub.list[[i]] = ASV.sub.list[[i]][matched, ]
}
lapply(ASV.sub.list, dim)
dim(ASV.ref.relAb)

# --- Work only with ASVs observed in both datasets:
ASV.ref.relAb.tmp = ASV.ref.relAb
ASV.ref.sub.tmp = ASV.ref.sub
for(i in 2:length(ASV.sub.rel.list)){
 matched = match(colnames(ASV.ref.relAb.tmp), colnames(ASV.sub.rel.list[[i]]))
 ASV.ref.relAb.tmp = ASV.ref.relAb.tmp[, !is.na(matched)]
 ASV.ref.sub.tmp = ASV.ref.sub.tmp[, !is.na(matched)]
}
dim(ASV.ref.relAb.tmp)
for(i in 2:length(ASV.sub.rel.list)){
  matched = match(colnames(ASV.ref.relAb.tmp), colnames(ASV.sub.rel.list[[i]]))
  ASV.sub.rel.list[[i]] = ASV.sub.rel.list[[i]][, matched]
  ASV.sub.list[[i]] = ASV.sub.list[[i]][, matched]
}
ASV.ref.relAb = ASV.ref.relAb.tmp
ASV.ref.sub = ASV.ref.sub.tmp
lapply(ASV.sub.rel.list, dim)
lapply(ASV.sub.list, dim)
dim(ASV.ref.relAb)

# Start computation -------
# --- Compute the transformation matrix
Pbas =  ASV.ref.sub
Qbas = ASV.sub.list[[Rep+1]]
P = ASV.ref.relAb # Starting community relative abundance
Q = ASV.sub.rel.list[[Rep+1]] # Against replicate 4, final one
dim(P); dim(Q)
superimp = kabsch_R(Q, P) # Find optimal rotation

# ..... extract quantities
(rmsd = superimp$rmsd) # 0.4860083
U = superimp$U
Prot = superimp$Prot 
rownames(Prot) = rownames(P)
colnames(Prot) = colnames(P)
rownames(U) = colnames(P)
colnames(U) = colnames(P)

# --- Compute the covariance between starting and final matrix
# ... here we can't consider relative abundances so we resample
# ... the data considering the relative abundances
Q.mult = matrix(NA, nrow = nrow(Q), ncol = ncol(Q))
P.mult = matrix(NA, nrow = nrow(P), ncol = ncol(P))
colnames(Q.mult) = colnames(P.mult) =  colnames(Q)
rownames(Q.mult) = rownames(Q)
rownames(P.mult) = rownames(P)
Q.rel.rnd = matrix(NA, nrow = nrow(Q), ncol = ncol(Q))
P.rel.rnd = matrix(NA, nrow = nrow(P), ncol = ncol(P))
colnames(Q.rel.rnd) = colnames(P.rel.rnd) =  colnames(Q.rel.rnd)
rownames(Q.rel.rnd) = rownames(Q)
rownames(P.rel.rnd) = rownames(P)
size = 10000
for(i in 1:dim(Q)[1]){
  alpha = Qbas[i, ]
  Q.dir = rdirichlet(1, alpha = alpha)
  Q.rel.rnd[i, ] =  Q.dir
  Q.mult[i, ] = rmultinom(n = 1, size = size, prob = Q.dir)
  alpha = Pbas[i, ]
  P.dir = rdirichlet(1, alpha = alpha)
  P.rel.rnd[i, ] = P.dir
  P.mult[i, ] = rmultinom(n = 1, size = size, prob = P.dir)
}
probs = c(0.9,0.95,0.97,1)
quantile(P.mult, probs = probs)
quantile(Q.mult, probs = probs)
QP.mult = rbind(Q.mult, P.mult)
QP.mult = QP.mult[, which(colSums(QP.mult) > 0)]

# ... Reshape the original Q and P
matched = match(colnames(QP.mult), colnames(Q))
Q.mult = Q.mult[, matched]
P.mult = P.mult[, matched]
colnames(Q.mult) = paste0(colnames(Q.mult),"_s")
colnames(P.mult) = paste0(colnames(P.mult),"_e")

# ... Rescale
#Cross = aitchinson(Q, P) # rescales internally, comment rescale if you use it
#Q.mult = log10(Q.mult + pseudocount)
#P.mult = log10(P.mult + pseudocount)
#Cross.cov = cov(Q.mult, P.mult)
#Cross  = cor(Q.mult, P.mult)
#stop()
# quantile(Cross, na.rm = T)
# quantile(Cross, na.rm = T, probs = probs)

# .... Write output
setwd(dirOut)
fileProt = paste0("Rotated_communities_Rep",Rep,".tsv")
write.table(Prot, file=fileProt, sep = "\t", quote = FALSE)

fileTrans = paste0("Transformation_matrix_Starting2Rep",Rep,".tsv")
write.table(U, file=fileTrans, sep = "\t", quote = FALSE)

# fileCross = paste0("CrossCovariance_matrix_Starting2Rep",Rep,".tsv")
# write.table(Cross, file=fileCross, sep = "\t", quote = FALSE)

# Bootstrap ----------
if(compute_CI == T){ # Estimate confidence intervals, this routine may take several minutes
  # --- Estimate a confidence interval for rmsd bootstrapping
  Nrnd = 50
  set.seed(131222)
  # ... Create a function to use boot function, we want to directly recover the rmsd
  rmsd_fun_boot = function(data, indices, Q){
    d = data[indices, ] # bootstrap samples in both matrices (paired)
    Q = Q[indices, ]
    superimp = kabsch_R(Q, d)
    rmsd = superimp$rmsd
    return(rmsd)
  }
  rmsd_fun_halfboot = function(data, indices, Q){
    d = data[indices, ] # bootstrap samples only in starting
    superimp = kabsch_R(Q, d)
    rmsd = superimp$rmsd
    return(rmsd)
  }
  rmsd_fun_rnd = function(data, Q){
    # randomize all cells so the result should be fully random
    indices = sample(seq(1, dim(data)[1]*dim(data)[2]))
    d = matrix(data[indices], nrow = dim(data)[1], 
               ncol = dim(data)[2]) # 
    superimp = kabsch_R(Q, d)
    rmsd = superimp$rmsd
    return(rmsd)
  }
  rmsd_fun_fullrnd = function(data, Q){
    # generate a fully random matrix
    Nrand = dim(data)[1]*dim(data)[2]
    d = matrix(rnorm(Nrand, mean = 0.5, sd = 0.1), nrow = dim(data)[1], 
               ncol = dim(data)[2]) # 
    superimp = kabsch_R(Q, d)
    rmsd = superimp$rmsd
    return(rmsd)
  }
  # -- Bootstrap and randomize (this takes several minutes)
  reps_boot = boot(data=P, statistic=rmsd_fun_boot, R=Nrnd, Q = Q)
  reps_halfboot = boot(data=P, statistic=rmsd_fun_halfboot, R=Nrnd, Q = Q)
  reps_rand = vector(mode = "numeric", length = Nrnd)
  reps_full_rand = vector(mode = "numeric", length = Nrnd)
  for(i in 1:Nrnd){
    reps_rand[i] = rmsd_fun_rnd(P, Q)
    reps_full_rand[i] = rmsd_fun_fullrnd(P, Q)
  }

  # ... write summary
  fileBootRMSD = paste0("Summary_Bootstrap_rmsd_summary_StartingVsRep",Rep,".txt")
  sink(file = fileBootRMSD)
  (quantile(reps_boot$t,probs =  c(0.95,0.5,0.05))) # [0.47, 0.44]
  (rmsd) # 0.48
  sink()
  fileHalfBootRMSD = paste0("Summary_BootstrapHalf_rmsd_summary_StartingVsRep",Rep,".txt")
  sink(file = fileHalfBootRMSD)
  (quantile(reps_halfboot$t,probs =  c(0.95,0.5,0.05))) # [0.54, 0.51]
  (rmsd) # 0.48
  sink()
  fileRandRMSD = paste0("Summary_Rand_rmsd_summary_StartingVsRep",Rep,".txt")
  sink(file = fileRandRMSD)
  (quantile(reps_rand,probs =  c(0.95,0.5,0.05))) # [0.559, 0.542]
  (rmsd) # 0.48
  sink()
  fileFullRandRMSD = paste0("Summary_FullRand_rmsd_summary_StartingVsRep",Rep,".txt")
  sink(file = fileFullRandRMSD)
  (quantile(reps_full_rand, probs =  c(0.95,0.5,0.05))) # [2.39, 2.38]
  (rmsd) # 0.48
  sink()
  # ... plot results
  fileBootPlot = paste0("Plot_bootstrap_rmsd_summary_StartingVsRep",Rep,".txt")
  pdf(file = fileBootPlot)
  plot(reps_boot)
  dev.off()
  fileHalfBootPlot = paste0("Plot_halfbootstrap_rmsd_summary_StartingVsRep",Rep,".txt")
  pdf(file = fileHalfBootPlot)
  plot(reps_halfboot)
  dev.off()
  #stop()
  # --- Estimate a confidence interval for the cross-covariance
  # Cross.boot.mean = matrix(0, nrow = nrow(Cross), ncol = ncol(Cross))
  # Cross.boot.var = matrix(0, nrow = nrow(Cross), ncol = ncol(Cross))
  # for(i in 1:Nrnd){
  #   id.boot = sample(nrow(Q), replace = TRUE)
  #   Q.boot = Q.mult[id.boot,]
  #   P.boot = P.mult[id.boot,]
  #   Cross.boot = cor(Q.boot, P.boot)
  #   Cross.boot.mean = Cross.boot.mean + Cross.boot
  #   Cross.boot.var = Cross.boot.var + Cross.boot^2
  # }
  # Cross.boot.mean = Cross.boot.mean/Nrnd
  # quantile(Cross.boot.mean, na.rm = T) # double check within expected limits
  # Cross.boot.sd = sqrt(Cross.boot.var/Nrnd - Cross.boot.mean^2)
  # quantile(Cross.boot.sd, na.rm = T)
  # Cross.sign = sign(Cross)
  # Cross.boot.CI = Cross.boot.mean + Cross.sign * 2 * Cross.boot.sd
  # quantile(Cross.boot.CI, na.rm = T)
  # Cross[which(abs(Cross) < abs(Cross.boot.CI))] = 0 # remove non-significant
  # length(which(Cross != 0)) # still too many, 1M
}
#quantile(Cross, probs = c(0,0.25, 0.5, probs), na.rm = TRUE)

# .... modify cross covariance to work with a network
# Cross[lower.tri(Cross, diag = T)] = 0
# dim(Cross)
# Cross_long = melt(Cross)
# id.nonzero = which(abs(Cross_long$value) > 0)
# Cross_long = Cross_long[id.nonzero, ]
# quantile(Cross_long$value)
# qqnorm(Cross_long$value)
# hist = hist(Cross_long$value, breaks = 50, )
# mean = mean(Cross_long$value)
# stdv = sd(Cross_long$value)
# z.score = (Cross_long$value - mean)/stdv
# quantile(z.score)
# z.score.sig = which(abs(z.score) > 5)
# #stop()
# Cross_long = Cross_long[z.score.sig, ]
# Cross_long$Type = NA
# Cross_long$Type[which(Cross_long$value > 0)] = "Pos"
# Cross_long$Type[which(Cross_long$value < 0)] = "Neg"
# 
# #fileTrans = paste0("Network_transformation_matrix_Starting2Rep",Rep,".tsv") # This was U_long, no longer computed
# fileCross = paste0("Network_CrossCovariance_matrix_Starting2Rep",Rep,".tsv")
# write.table(Cross_long, file=fileCross, sep = "\t", quote = FALSE, row.names = F)

# Independent SVD to extract eigenvectors  --------
#ord.list = list()
k = 0
# ... ordination for the replicates that were not used in the transformation
fileSumm = "Summary_SVD_replicates_to_startRotated.txt"
sink(fileSumm) # start sink to file
cat("# ** SVD of independent replicates \n")
Qc.list = list()
for(i in 2:length(ASV.sub.rel.list)){ #replicates 1, 2, 3, 4
  k = k+1
  cat(paste0("# --- Final replicate: ",k,"\n"))
  Qtmp = ASV.sub.rel.list[[i]]
  Qc = scale(Qtmp, center = TRUE, scale = FALSE)
  Qc.list[[i]] = Qc 
  cat("# RMSD between final replicate and transformed starting communities \n")
  Dsr.sq = (Prot - Qc)**2
  rmsd = sqrt(sum(apply(Dsr.sq, 1, sum)) / dim(Dsr.sq)[1]);
  print(rmsd)
  ord = svd(Qc)
  variance.explained = prop.table(ord$d^2)
  cat("# variance explained with the SVD (5 first components) \n")
  print(variance.explained[1:5])
  sv = as.data.frame(ord$u[, c(1,2,3,4)]) # first and second singular vector
  col_names = c("sv1", "sv2","sv3","sv4")
  col_names = paste0("Rep",k,"_",col_names)
  colnames(sv) = col_names 
  if(k == 1){
    sv.df = sv
  }else{
    sv.df = cbind(sv.df,sv)
  }
}

# .... the following was done for the representations below
sv.df$Rep3_sv1 = (-1)*sv.df$Rep3_sv1
sv.df$Rep3_sv2 = (-1)*sv.df$Rep3_sv2

# ... transformed data
cat("# SVD of transformed starting community \n")
ord.Prot = svd(Prot)
variance.explained = prop.table(ord.Prot$d^2)
cat("# variance explained with the SVD (5 first components) \n")
(variance.explained[1:5])
sv = as.data.frame(ord.Prot$u[, c(1,2,3,4)]) # first and second singular vector
col_names = c("sv1", "sv2","sv3", "sv4")
col_names = paste0("Rot_",col_names)
colnames(sv) = col_names
sv.df = cbind(sv.df, sv)
colnames(sv.df) 

# ... Compare if eigenvectors are in the same order are the same
Neig = dim(sv.df)[2]
prod = matrix(data = 0, nrow = Neig, ncol = Neig) 
for(i in 1:Neig){
  for(j in i:Neig){
    prod[i, j] = sv.df[, i] %*% sv.df[, j]
  }
}
colnames(prod) = colnames(sv.df)
rownames(prod) = colnames(sv.df)
cat("# All-against-all scalar product between left eigenvectors \n")
options(max.print=500)
(prod)
sink() # end sink
options(max.print=50)

# --- Plots
sv.df.gg = reshape_df_to_ggplot(sv.df, x.vec = c("Rot_sv1"),
                     y.vec = c("Rep1_sv1","Rep2_sv1","Rep3_sv1"),
                      char.list = list(c("Rep. 1","Rep. 2", "Rep. 3")))

colnames(sv.df.gg)
xlab = "SVD, component 1 (prediction)"
ylab = "SVD, component 1 (experiments)"
filePlot1 = "Plot_SVD1_predictionVsExperiments_tmp.pdf"
pdf(file = filePlot1, width = 10)
g = ggplot(sv.df.gg,aes(x = x.out, y = y.out, color = new_factor1))+
  geom_point(alpha = 0.7, size = 5)+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.y.npc = c(1.00,0.99, 0.97),
           size = 6,show.legend = F)+ 
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #         label.y.npc = c(0.99, 0.94, 0.90))+ 
  #stat_regline_equation(label.y.npc = c(0.95, 0.91, 0.87))+ 
  scale_color_manual(labels = c("1", "2", "3", "4"),
                     values = c(2, 3, 4, 5))+
  xlab(xlab)+ylab(ylab)+
  labs(color = "Replicates")+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
print(g)
dev.off()

# .... Same for component 2
sv.df.gg = reshape_df_to_ggplot(sv.df, x.vec = c("Rot_sv1"),
                                y.vec = c("Rep1_sv2","Rep2_sv2","Rep3_sv2"),
                                char.list = list(c("Rep. 1","Rep. 2", "Rep. 3")))

colnames(sv.df.gg)
xlab = "SVD, component 1 (prediction)"
ylab = "SVD, component 2 (experiments)"
filePlot2 = "Plot_SVD1_prediction_Vs_SVD2_Experiments.pdf"
pdf(file = filePlot2, width = 10)
g = ggplot(sv.df.gg,aes(x = x.out, y = y.out, color = new_factor1))+
  geom_point(alpha = 0.7, size = 5)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.y.npc = c(1.00,0.99, 0.97),
           size = 6,show.legend = F)+ 
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #         label.y.npc = c(0.99, 0.94, 0.90))+ 
  #stat_regline_equation(label.y.npc = c(0.95, 0.91, 0.87))+ 
  scale_color_manual(labels = c("1", "2", "3", "4"),
                     values = c(2, 3, 4, 5))+
  xlab(xlab)+ylab(ylab)+
  labs(color = "Replicates")+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
print(g)
dev.off()


# SVD of all samples together and facet -----
# Start creating a single centered object, this can be done because
# the function scale centers each column with respect to each own mean,
# so concatenating centered matrices is harmless. This is a better strategy
# because performing independent SVDs does not allow for pairing
# eigenvectors across them, since the order could be different

Mc = Prot
for(i in 1:length(Qc.list)){
  Mc = rbind(Mc, Qc.list[[i]])
}
dim(Mc)

# --- Perform SVD
ord = svd(Mc)
variance.explained = prop.table(ord$d^2)
cat("# variance explained with the SVD (5 first components) \n")
print(variance.explained[1:5])
sv = as.data.frame(ord$u[, c(1,2)]) 

# ...... reshape for ggplot
Nsubset = 5 # 4 replicates plus transformed data
Lsubset = dim(sv)[1]/Nsubset
# ..... extract first eigenvector 
sv.long.x = rep(sv$V1[1:Lsubset], 3)# 4) # transformed data (x axis)
sv.long.y = c(sv$V1[(Lsubset+1):(2*Lsubset)], # 4 exp replicates (y axis)
              sv$V1[(2*Lsubset+1):(3*Lsubset)],
              sv$V1[(3*Lsubset+1):(4*Lsubset)])#,
              #sv$V1[(4*Lsubset+1):(5*Lsubset)])
sv.long.fac = c(rep("1", Lsubset),
                rep("2", Lsubset),
                rep("3", Lsubset))#,
                #rep("4", Lsubset))
sv.long = data.frame(sv.long.x, sv.long.y, sv.long.fac)

xlab = "SVD, component 1 (prediction)"
ylab = "SVD, component 1 (experiments)"
filePlot3 = "Plot_SVD1_prediction_Vs_SVD1_Experiments_sameSVD.pdf"
pdf(file = filePlot3, width = 10)
g = ggplot(sv.long, aes(x = sv.long.x, y =sv.long.y, color = sv.long.fac))+
  geom_point(alpha = 0.7, size = 5)+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.y.npc = c(1.00,0.99, 0.97),
           size = 6,show.legend = F)+ 
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #         label.y.npc = c(0.99, 0.94, 0.90))+ 
  #stat_regline_equation(label.y.npc = c(0.95, 0.91, 0.87))+ 
  scale_color_manual(labels = c("1", "2", "3", "4"),
                     values = c(2, 3, 4, 5))+
  xlab(xlab)+ylab(ylab)+
  labs(color = "Replicates")+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
print(g)
dev.off()

# ..... extract second eigenvector 
sv.long.x = rep(sv$V2[1:Lsubset], 3)# 4) # transformed data (x axis)
sv.long.y = c(sv$V2[(Lsubset+1):(2*Lsubset)], # 4 exp replicates (y axis)
              sv$V2[(2*Lsubset+1):(3*Lsubset)],
              sv$V2[(3*Lsubset+1):(4*Lsubset)])#,
#sv$V1[(4*Lsubset+1):(5*Lsubset)])
sv.long.fac = c(rep("1", Lsubset),
                rep("2", Lsubset),
                rep("3", Lsubset))#,
#rep("4", Lsubset))
sv.long = data.frame(sv.long.x, sv.long.y, sv.long.fac)

xlab = "SVD, component 2 (prediction)"
ylab = "SVD, component 2 (experiments)"
filePlot4 = "Plot_SVD2_prediction_Vs_SVD2_Experiments_sameSVD.pdf"
pdf(file = filePlot4, width = 10)
g = ggplot(sv.long, aes(x = sv.long.x, y =sv.long.y, color = sv.long.fac))+
  geom_point(alpha = 0.7, size = 5)+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.y.npc = c(1.00,0.99, 0.97),
           size = 6,show.legend = F)+ 
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #         label.y.npc = c(0.99, 0.94, 0.90))+ 
  #stat_regline_equation(label.y.npc = c(0.95, 0.91, 0.87))+ 
  scale_color_manual(labels = c("1", "2", "3", "4"),
                     values = c(2, 3, 4, 5))+
  xlab(xlab)+ylab(ylab)+
  labs(color = "Replicates")+
  theme_bw()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
print(g)
dev.off()

