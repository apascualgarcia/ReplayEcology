##################################################
# remove_problematicreads.R
##################################################

# In this script I remove a very small number of problematic 
# reads in the fastq files that were preventing downstream processing
# For context, a very small minority of the lines (i.e. 4 lines) in each 
# of the demultiplexed FASTQ files had quality scores that were a different 
# length to the sequences. This is odd but seems to be how the data arrived 
# from the sequencing company. The problem seems restricted to day 0 #
#data/first sequencing run. Online discussions suggest it is due to some 
# file corruption (https://www.biostars.org/p/180310/; 
# https://www.biostars.org/p/231090/; 
# https://forum.qiime2.org/t/fastq-gz-and-quality-score-length-do-not-match-using-type-emppairedendsequences-from-miseq/14142/7), 
# perhaps occurring at the sequencing center, and that these records can be removed. 
#Therefore, we removed this very small minority of problematic reads from the demultiplexed 
# FASTQs in order to enable downstream analysis with these FASTQs.

# Cardiff, August 2024
# European Centre for the Environment & Human Health, University of Exeter
# Matt Lloyd Jones
# https://github.com/befriendabacterium/

# READ IN FASTQ FILE -----------------------------------------------------------

#make a vector of the folders containing the seqs in the demultiplexed directory
folders<-c('day0','day7')

for (d in 1:2){

  seq.files<-list.files(paste('2_demultiplexed',folders[d],sep='/'))
  
for (i in 1:length(seq.files)){

fastq_temp<-readLines(paste('2_demultiplexed',folders[d],seq.files[i], sep='/'))

# COMPARE SEQ AND QUAL LENGTHS -----------------------------------------------------------
#code adapted from https://github.com/benjjneb/dada2/issues/654

nc<-sapply(fastq_temp,nchar)
seqs.len <- nc[seq(2, length(nc), 4)]
quals.len <- nc[seq(4, length(nc), 4)]
diffs<-as.numeric(seqs.len-quals.len)
table(diffs)
nomatch<-which(diffs!=0)*4

fastq_temp[nomatch]

# REMOVE PROBLEMATIC READS ------------------------------------------------

#if errors are found, remove them (if not, the uncorrected fastq will be written)
if(length(nomatch)!=0){
toremove_entries<-sort(as.numeric(sapply(nomatch,function(x)seq(x,x-3))))
fastq_temp<-fastq_temp[-toremove_entries]
}

#write the corrected fastq
write.table(fastq_temp,
            paste('2_demultiplexed/corrected',folders[d],seq.files[i], sep='/'),
            row.names = F, col.names = F, quote=F)

  }
}
  
