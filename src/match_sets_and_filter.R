rm(list=ls())
#setwd("/home/apascual/Nextcloud/Research/Projects/FunctionalGroups/Repositories/convergence/Partial_Matt_pipeline/src")
library(stringi)
library(stringr)
# --- Directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
setwd(dirSrc)

# --- Load objects from Dada
seqtab_treeholes <-readRDS('../4_dada2/seqtab.nochim.RDS')
#sample df table froM DADA2
samdf_treeholes <-read.csv('../4_dada2/samdf.csv')
taxa_wsp<-read.csv('../4_dada2/taxa_wsp.csv')

#  --- Load metadata
sample_md<-read.table('../4_dada2/metadata_Time0D-7D-4M_May2022.csv',
                      sep="\t",header=TRUE)
head(sample_md)[1:5,]
nrow(seqtab_treeholes)
nrow(samdf_treeholes) # 2843
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
samdf_treeholes_4M=samdf_treeholes[id_4M,]
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
samdf_treeholes_rest=samdf_treeholes[id_rest,]
rownames(seqtab_treeholes_rest)[c(1:10,1000:1010)]

# --- Merge both 4M and remainder dfs
seqtab_treeholes=rbind(seqtab_treeholes_4M,seqtab_treeholes_rest)
samdf_treholes=rbind(seqtab_treeholes_4M,seqtab_treeholes_rest)

dim(seqtab_treeholes)
dim(samdf_treholes)

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

#list communities in samdf which are still in seqtab_treeholes after filtering
samples_to_keep<-which(samdf_treeholes$filename%in%rownames(seqtab_treeholes))
#filter them out
samdf_treeholes<-samdf_treeholes[samples_to_keep,]

#list sequences/ASVs in samdf which are still in seqtab_treeholes after filtering
sequences_to_keep<-which(taxa_wsp$X%in%colnames(seqtab_treeholes))
#filter them out
taxa_wsp<-taxa_wsp[sequences_to_keep,]

# WRITE NEW SEQTABLE AS CSV AND FASTA -----------------------------------------------------------

ASV_sequences<-as.list(colnames(seqtab_treeholes))
ASV_names<-paste('ASV',1:ncol(seqtab_treeholes), sep='_')

taxa_wsp<-cbind(ASV_names,taxa_wsp)
colnames(taxa_wsp)[which(colnames(taxa_wsp)=='X')]<-'sequence'
write.table(taxa_wsp,'../6_finalfiles/taxa_wsp_readyforanalysis.csv', row.names=F,quote = FALSE)

#write a fasta file
seqinr::write.fasta(sequences=ASV_sequences, names=ASV_names,'../6_finalfiles/seqtable_readyforanalysis.fasta', open = "w", nbchar = 60, as.string = FALSE)

#rename the sequence colnames in the seqtab to shorthand ASV names
colnames(seqtab_treeholes)<-ASV_names
# rename the samples to match the metadata
rownames(seqtab_treeholes)[1:10]
write.table(seqtab_treeholes,'../6_finalfiles/seqtable_readyforanalysis.csv',sep="\t",quote=FALSE)

