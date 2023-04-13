here::here()

# READ IN DATA ------------------------------------------------------------

#sequence table from DADA2
seqtab_nochim<-readRDS('4_dada2/seqtab.nochim.RDS')
#sample df table froM DADA2
samdf<-read.csv('4_dada2/samdf.csv')

#check same row numbers
nrow(seqtab_nochim)
nrow(samdf)

#sample metadata
sample_md<-read.csv('4_dada2/metadata_Time0-Time7-SJD_merged_2021.4stamp.csv')

# FILTER OUT RELEVANT SAMPLES FROM SEQUENCE TABLE ----------------------------------------------------------

stringstofind<-paste(sample_md$sampleid, collapse = '|')
samples_needed<-grep(stringstofind, rownames(seqtab_nochim))

#filter sequence table to relevant samples
seqtab_nochim_filtered<-seqtab_nochim[samples_needed,]
nrow(seqtab_nochim_filtered)

#filter sequence table to relevant samples
samdf_filtered<-samdf[samples_needed,]
nrow(samdf_filtered)

# WRITE DATAFRAMES --------------------------------------------------------

write.csv(samdf_filtered,'5_treeholes/samdf_treeholes.csv', row.names=F)
saveRDS(samdf_filtered,'5_treeholes/samdf_treeholes.RDS')

write.csv(seqtab_nochim_filtered,'5_treeholes/seqtab_treeholes.csv', row.names=F)
saveRDS(seqtab_nochim_filtered,'5_treeholes/seqtab_treeholes.RDS')

