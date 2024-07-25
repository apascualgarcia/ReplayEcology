# INSTALL DADA2 -----------------------------------------------------------

#NB TO USE DADA 3.14 AND REPRODUCE THE PIPELINE EXACTLY, YOU NEED TO FIRST INSTALL AND RUN R VERSION 4.1 

#unhash to install dada2 3.14
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.14")

# GETTING READY: LOAD DADA2 ---------------------------------------------------

#loads DADA2 and returns the package version you're running (should be 3.14 for complete reproducibility)
library(dada2); packageVersion("dada2")

# GETTING READY: READ IN FILES --------------------------------------------

path <- '1_demultiplexed' #directory containing the demultiplexed fastq files that were  deposited at NCBI under BioProject accession number PRJNA989519
path
fns <- list.files(path, pattern="fastq.gz", recursive = T) # CHANGE if different file extensions
length(fns)

# GETTING READY: EXRACT SAMPLE NAMES --------------------------------------

#write file names to sample.names
sample.names<-sub("^(.*)[.].*", "\\1", fns)

#remove the extension
sample.names<-sub("^(.*)[.].*", "\\1", sample.names)

#replace forward slash with _
sample.names<-gsub("\\/", "_", sample.names)

# INSPECT READ QUALITY PROFILES -------------------------------------------

#plot read quality profiles
qualplot_agg<-plotQualityProfile(file.path(path,fns), aggregate=T)
#save to 3_dada2 folder to save time later
saveRDS(qualplot_agg,'4_dada2/qualplot_agg.RDS')

# FILTER AND TRIM: PLACE FILTERED FILES IN SUB-DIRECTORY ---------------------------------------------------------

filts<-file.path('3_filtered', paste0(sample.names, "_filt.fastq.gz"))
filts

# FILTERING ---------------------------------------------------------------

out <- filterAndTrim(file.path(path,fns), filts,
                     truncLen=240, 
                     maxN=0, maxEE=1, truncQ=11, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
#save output of filtering to RDS
saveRDS(out, '3_filtered/filteredandtrimmed.RDS')

#read output of filtering in (this allows you to avoid running filterAndTrim() again
out<-readRDS('3_filtered/filteredandtrimmed.RDS')

#subset three objects for only samples hat were actually filtered - see https://github.com/benjjneb/dada2/issues/1279
filts<-filts[file.exists(filts)]
out<-out[which(file.exists(filts)),]
fns<-fns[file.exists(filts)]

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

#get the final list of sample names for the samples after chimera removal (sample as lines 22-32)
samples.out <- rownames(seqtab.nochim)
samples.out <- gsub("\\.", "_", samples.out)
samples.out <- sub("\\-.*", "", samples.out)

#get the sampling days for each sample
day <- sapply(strsplit(samples.out, "_"), `[`, 1)
day <- as.numeric(gsub("day","",day))

#get the 96 well plate well for each sample
well <- sapply(strsplit(samples.out, "_"), `[`, 2)

#get the community name for each sample
community <- sapply(strsplit(samples.out, "_"), `[`, 3)

#get the experimental replicate (1-4) for each sample
rep <- as.numeric(sapply(strsplit(samples.out, "_"), `[`, 4))

#get the sequence run ID for each sample
seq.run<-c()
seq.run[grep('052214DR16s', samples.out)]<-'052214DR16s'
seq.run[grep('021216DR515F', samples.out)]<-'021216DR515F'
levels(as.factor(seq.run)) #check number of runs (should be 2)

#build the metadata dataframe based on above
samdf <- data.frame(community=community, day=day, rep=rep, seq.run=seq.run, well=well)
#save it as a CSV
write.csv(samdf,'4_dada2/samdf.csv')
