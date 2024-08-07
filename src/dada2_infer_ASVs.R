#list the files in the 8 sub-folders
filts<-list.files(c('3_filtered/thisstudy_day0', 
                  '3_filtered/thisstudy_day7r1',
                  '3_filtered/thisstudy_day7r2',
                  '3_filtered/thisstudy_day7r3',
                  '3_filtered/thisstudy_day7r4',
                  '3_filtered/scheuerlstudy',
                  '3_filtered/mombrikotbstudy',
                  '3_filtered/misc'),
                  pattern="fastq.gz", recursive = F) 

# LEARN THE ERROR RATE (takes a long time) ----------------------------------------------------

#learn error rates
err <- learnErrors(filts, multithread=FALSE)
#save output of learning error rates i
saveRDS(err, '4_dada2/err.RDS')
#read output of learning error rates in (this allows you to avoid running learnErrors() again
err<-readRDS('4_dada2/err.RDS')

#hashed out cos takes ages and can crash pipeline
plotErrors(err, nominalQ=TRUE)

# SAMPLE INFERENCE --------------------------------------------------------

#run DADA2's main sample/ASV inference function
dada <- dada(filts, err=err, multithread=FALSE)
#save output of sample inference
saveRDS(dada,'4_dada2/dada_output.RDS')
#read output of sample inference in (this allows you to avoid running dada() again
dada<-readRDS('4_dada2/dada_output.RDS')

# CONSTRUCT SEQUENCE TABLE ------------------------------------------------

#construct sequence table
seqtab <- makeSequenceTable(dada)
#check dimensions to get number of samples and ASVs
dim(seqtab)
#save sequence table as RDS
saveRDS(seqtab, '4_dada2/seqtab.RDS')
#save sequence table as CSV for more human-readable format
write.csv(seqtab,'4_dada2/seqtab.csv')

#read in sequence table (this allows you to avoid running makeSequenceTable again)
seqtab<-readRDS('4_dada2/seqtab.RDS')

# REMOVE CHIMERAS https://github.com/benjjneb/dada2/issues/1430 ---------------------------------------------------------

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
#check dimensions to get number of samples and ASVs
dim(seqtab.nochim)
#save sequence table as RDS
saveRDS(seqtab.nochim, '4_dada2/seqtab.nochim.RDS')
#save sequence table as CSV for more human-readable format
write.csv(seqtab.nochim,'4_dada2/seqtab.nochim.csv')

#read in sequence table (this allows you to avoid running removeBimeraDenovo again)
seqtab.nochim<-readRDS('4_dada2/seqtab.nochim.RDS')

# TRACK READS THROUGH THE PIPELINE ----------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)

# ASSIGN TAXONOMY ---------------------------------------------------------

#assign taxonmy to ASVs using SILVA database
taxa <- assignTaxonomy(seqtab.nochim, "4_dada2/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
#save taxonomy table as RDS
saveRDS(taxa, '4_dada2/taxa.RDS')
#save taxonomy table as CSV for more human-readable format
write.csv(taxa, '4_dada2/taxa.csv')
#read in taxonomy table (this allows you to avoid running assignTaxonomy() again
taxa<-readRDS('4_dada2/taxa.RDS')

# ASSIGN SPECIES-LEVEL TAXONOMY ------------------------------------------------------

#read in species level SILVA database
species_fasta<-"4_dada2/silva_species_assignment_v138.1.fa.gz"

# There is some memory issue, causing the addSpecies to fail (https://github.com/benjjneb/dada2/issues/733)
# So we cut the table in chunks:
chunk.size <- 2000
chunks <- split(c(1:nrow(taxa)),
                sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size)))
chunks.species <- lapply(chunks,
                         function(x){
                           return(addSpecies(taxa[x,],
                                             refFasta = species_fasta, verbose = TRUE))
                         })
taxa_wsp <- do.call(rbind, chunks.species)

#save species-level taxonomy table as RDS
saveRDS(taxa_wsp, '4_dada2/taxa_wsp.RDS')
#save species-level taxonomy table as CSV for more human-readable format
write.csv(taxa_wsp, '4_dada2/taxa_wsp.csv')

# Removing sequence rownames for display only
taxa.print <- taxa 
rownames(taxa.print) <- NULL
#inspect the output
head(taxa.print)

# WRITE METADATA ----------------------------------------------------------
# 
# #get the final list of sample names for the samples after chimera removal (sample as lines 22-32)
# samples.out <- rownames(seqtab.nochim)
# samples.out <- gsub("\\.", "_", samples.out)
# samples.out <- sub("\\-.*", "", samples.out)
# 
# #get the sampling days for each sample
# day <- sapply(strsplit(samples.out, "_"), `[`, 1)
# day <- as.numeric(gsub("day","",day))
# 
# #get the 96 well plate well for each sample
# well <- sapply(strsplit(samples.out, "_"), `[`, 2)
# 
# #get the community name for each sample
# community <- sapply(strsplit(samples.out, "_"), `[`, 3)
# 
# #get the experimental replicate (1-4) for each sample
# rep <- as.numeric(sapply(strsplit(samples.out, "_"), `[`, 4))
# 
# #get the sequence run ID for each sample
# seq.run<-c()
# seq.run[grep('052214DR16s', samples.out)]<-'052214DR16s'
# seq.run[grep('021216DR515F', samples.out)]<-'021216DR515F'
# levels(as.factor(seq.run)) #check number of runs (should be 2)
# 
# #build the metadata dataframe based on above
# samdf <- data.frame(community=community, day=day, rep=rep, seq.run=seq.run, well=well)
# #save it as a CSV
# write.csv(samdf,'4_dada2/samdf.csv')
# #save it as a CSV
# samdf<-read.csv('4_dada2/samdf.csv', row.names = 1)
