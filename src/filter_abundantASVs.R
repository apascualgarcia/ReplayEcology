here::here()
library(stringi)

# READ IN DATA ------------------------------------------------------------

#sequence table from DADA2
seqtab_treeholes<-readRDS('5_treeholes/seqtab_treeholes.RDS')
samdf_treeholes<-read.csv('5_treeholes/samdf_treeholes.csv')
taxa_wsp<-read.csv('4_dada2/taxa_wsp.csv')

nrow(seqtab_treeholes)
nrow(samdf_treeholes)

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
write.table(taxa_wsp,'6_finalfiles/taxa_wsp_readyforanalysis.csv',sep="\t" ,row.names=F,quote = FALSE)

#write a fasta file
seqinr::write.fasta(sequences=ASV_sequences, names=ASV_names,'6_finalfiles/seqtable_readyforanalysis.fasta', open = "w", nbchar = 60, as.string = FALSE)

#rename the sequence colnames in the seqtab to shorthand ASV names
colnames(seqtab_treeholes)<-ASV_names
# rename the samples to match the metadata
rownames(seqtab_treeholes)[1:10]
names.tmp=sapply(rownames(seqtab_treeholes)
       ,FUN=function(x){unlist(stri_split(x,fixed="_"))[2]})
pat=regex("^([A-Z]+)([0-9]+)(\\.)")
names.tmp.tmp=gsub(pattern=pat,x=names.tmp,replacement = "")
length(names.tmp.tmp) #
length(unique(names.tmp.tmp))
rownames(seqtab_treeholes)=names.tmp.tmp
write.table(seqtab_treeholes,'6_finalfiles/seqtable_readyforanalysis.csv',sep="\t",quote=FALSE)
