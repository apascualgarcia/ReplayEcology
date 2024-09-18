##################################################
# match_sets_and_filter.R
##################################################

# This script performs an initial filtering of the sequence table down 
# to those samples that contain the names of communities used in this study. 
# It then filters out the least abundant ASVs 
# (to reduce number of ASVs/spurious ASVs), and removes samples with less 
# than 10K sequences, in line with our previous work.

# Zürich/Cardiff, April 2023/August 2024
# Theoretical Biology, ETH/E uropean Centre for the Environment & Human Health, University of Exeter
# Alberto Pascual Garcia/ Matt Lloyd Jones
# apascualgarcia.github.io / https://github.com/befriendabacterium/

rm(list=ls())
#setwd("/home/apascual/Nextcloud/Research/Projects/FunctionalGroups/Repositories/convergence/Partial_Matt_pipeline/src")
library(stringi)
library(stringr)
library(seqinr)

# --- Directories
this.dir = strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
dirSrc = paste(this.dir,"/src/",sep="") # Directory where the code is
dirDada2 = paste(this.dir,"/4_dada2/",sep="") # Dir of ASV table
dirOut = paste(this.dir,"/6_finalfiles/",sep="") # Dir of ASV table
setwd(dirSrc)

# --- Define objects

file.in.seq = 'seqtab.nochim.RDS'
file.in.taxa = 'taxa_wsp.csv'
file.in.meta = 'metadata_Time0D-7D-4M_May2022.csv'
file.out.taxa = 'taxa_wsp_matchedandfiltered.csv'
file.out.fasta = 'seqtable_matchedandfiltered.fasta'
file.out.seqtable = 'seqtable_matchedandfiltered.csv'

# --- Load objects from Dada

setwd(dirDada2)
seqtab_treeholes <-readRDS(file.in.seq)
#sample df table froM DADA2
taxa_wsp<-read.csv(file.in.taxa)

#  --- Load metadata
sample_md<-read.table(file.in.meta,sep="\t",header=TRUE)
head(sample_md)[1:5,]
nrow(seqtab_treeholes)
nrow(sample_md) # 2188

# MATCH, Remove undesired samples ----
# --- First, we isolate the relevant bit in the name of the file matching the metadata
names.tmp=sapply(rownames(seqtab_treeholes)
                 ,FUN=function(x){unlist(stri_split(x,fixed="_"))[2]})
rownames(seqtab_treeholes)=names.tmp

# --- Samples at 4M should give an exact match after this transformation
sample_md_4M=sample_md[which(sample_md$Experiment == "4M"),]
matched = match(sample_md_4M$sampleid,names.tmp)
length(matched[!is.na(matched)]) # 96 is correct
id_4M=matched[!is.na(matched)]

# .... we create 4M df
seqtab_treeholes_4M=seqtab_treeholes[id_4M,]
rownames(seqtab_treeholes_4M)[c(1:10,90:96)]

# --- Now we keep the remainder for further transformation
pat=regex("^([A-Z]+)([0-9]+)(\\.)")
names.tmp.tmp=gsub(pattern=pat,x=names.tmp,replacement = "")
rownames(seqtab_treeholes)=names.tmp.tmp

sample_md_rest=sample_md[which(sample_md$Experiment != "4M"),]
matched = match(sample_md_rest$sampleid,names.tmp.tmp)
length(matched[!is.na(matched)]) # 2165 is correct
id_rest=matched[!is.na(matched)]

# .... we create the remainder df
seqtab_treeholes_rest=seqtab_treeholes[id_rest,]
rownames(seqtab_treeholes_rest)[c(1:10,1000:1010)]

# --- Merge both 4M and remainder dfs
seqtab_treeholes=rbind(seqtab_treeholes_4M,seqtab_treeholes_rest)
dim(seqtab_treeholes)

# FILTER SEQUENCES --------------------------------------------------------

#remove OTUs with less than 100 reads across all samples
seqtab_treeholes<-seqtab_treeholes[,-which(colSums(seqtab_treeholes)<100)]
ncol(seqtab_treeholes)

#remove OTUs occurring in less than 10 samples
seqtab_treeholes<-seqtab_treeholes[,-which(colSums(seqtab_treeholes>0)<10)]
ncol(seqtab_treeholes)

#remove communities with less than 10,000 sequences
seqtab_treeholes<-seqtab_treeholes[which(rowSums(seqtab_treeholes[,-1])>10000),]
ncol(seqtab_treeholes)
nrow(seqtab_treeholes)

# FILTER THE SAMPLE METADATA AND TAXA TABLE ACCORDING TO REMAINING COMMUNITIES AND ASVS, respectively -------------------------------------------------------------------------

#list sequences/ASVs in taxa_wsp which are still in seqtab_treeholes after filtering
sequences_to_keep<-which(taxa_wsp$X%in%colnames(seqtab_treeholes))
#filter them out
taxa_wsp<-taxa_wsp[sequences_to_keep,]

# WRITE NEW SEQTABLE AS CSV AND FASTA -----------------------------------------------------------
setwd(dirOut)
ASV_sequences<-as.list(colnames(seqtab_treeholes))
ASV_names<-paste('ASV',1:ncol(seqtab_treeholes), sep='_')

taxa_wsp<-cbind(ASV_names,taxa_wsp)
colnames(taxa_wsp)[which(colnames(taxa_wsp)=='X')]<-'sequence'
write.table(taxa_wsp,file.out.taxa, row.names=F,quote = FALSE, sep = "\t")

#write a fasta file
seqinr::write.fasta(sequences=ASV_sequences, names=ASV_names,file.out.fasta, open = "w", nbchar = 60, as.string = FALSE)

#rename the sequence colnames in the seqtab to shorthand ASV names
colnames(seqtab_treeholes)<-ASV_names
# rename the samples to match the metadata
rownames(seqtab_treeholes)[1:10]
write.table(seqtab_treeholes,file.out.seqtable,sep="\t",quote=FALSE)

