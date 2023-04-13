here::here()

# FILTER BASED ON QUALITY CONTROL OF SEQUENCES --------------------------------------------------------

#remove ASVs with less than 100 reads across all samples
seqtab_nochim_filtered<-seqtab_nochim_filtered[,-which(colSums(seqtab_nochim_filtered)<100)]
ncol(seqtab_nochim_filtered)

#remove ASVs occurring in less than 10 samples
seqtab_nochim_filtered<-seqtab_nochim_filtered[,-which(colSums(seqtab_nochim_filtered>0)<10)]
ncol(seqtab_nochim_filtered)

#remove ASVs with less than 10,000 sequences
seqtab_nochim_filtered<-seqtab_nochim_filtered[which(rowSums(seqtab_nochim_filtered[,-1])>10000),]
ncol(seqtab_nochim_filtered)

rownames(seqtab_nochim_filtered)