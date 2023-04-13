# INSTALL DADA2 -----------------------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.14")

# GETTING READY: LOAD DADA2 ---------------------------------------------------

library(dada2); packageVersion("dada2")

# GETTING READY: READ IN FILES --------------------------------------------

path <- '2_demultiplexed/corrected' #CHANGE ME to the directory containing your demultiplexed fastq files
path
fns <- list.files(path, pattern="fastq.gz", recursive = T) # CHANGE if different file extensions
length(fns)

# GETTING READY: EXRACT SAMPLE NAMES --------------------------------------

sample.names<-sub("^(.*)[.].*", "\\1", fns)
sample.names<-sub("^(.*)[.].*", "\\1", sample.names)
sample.names

#sample.names<-stringr::str_extract(fns, "[^_]+")
#replace forward slash with _
sample.names<-gsub("\\/", "_", sample.names)
#replace period with _
#sample.names<-gsub("\\.", "_", sample.names)

# INSPECT READ QUALITY PROFILES -------------------------------------------

#qualplot_agg<-plotQualityProfile(file.path(path,fns), aggregate=T)
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
saveRDS(out, '3_filtered/filteredandtrimmed.RDS')

out<-readRDS('3_filtered/filteredandtrimmed.RDS')

#subset three objects for only samples hat were actually filtered - see https://github.com/benjjneb/dada2/issues/1279
filts<-filts[file.exists(filts)]
out<-out[which(file.exists(filts)),]
fns<-fns[file.exists(filts)]

# LEARN THE ERROR RATE (takes a long time) ----------------------------------------------------

err <- learnErrors(filts, multithread=FALSE)
saveRDS(err, '4_dada2/err.RDS')

#hashed out cos takes ages and can crash pipeline
plotErrors(err, nominalQ=TRUE)

# SAMPLE INFERENCE --------------------------------------------------------

dada <- dada(filts, err=err, multithread=FALSE)
saveRDS(dada,'4_dada2/dada_output.RDS')

# CONSTRUCT SEQUENCE TABLE ------------------------------------------------

seqtab <- makeSequenceTable(dada)
dim(seqtab)
saveRDS(seqtab, '4_dada2/seqtab.RDS')
seqtab<-readRDS('4_dada2/seqtab.RDS')
write.csv(seqtab,'4_dada2/seqtab.csv')

table(nchar(getSequences(seqtab)))

# REMOVE CHIMERAS https://github.com/benjjneb/dada2/issues/1430 ---------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, '4_dada2/seqtab.nochim.RDS')
seqtab.nochim<-readRDS('4_dada2/seqtab.nochim.RDS')
write.csv(seqtab.nochim,'4_dada2/seqtab.nochim.csv')

# TRACK READS THROUGH THE PIPELINE ----------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)

# ASSIGN TAXONOMY ---------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, "4_dada2/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
saveRDS(taxa, '4_dada2/taxa.RDS')
write.csv(taxa, '4_dada2/taxa.csv')
#taxa<-readRDS('4_dada2/taxa.RDS')

# ASSIGN SPECIES-LEVEL TAXONOMY ------------------------------------------------------

taxa<-readRDS('4_dada2/taxa.RDS')
species_fasta="4_dada2/silva_species_assignment_v138.1.fa.gz"

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

saveRDS(taxa_wsp, '4_dada2/taxa_wsp.RDS')
write.csv(taxa_wsp, '4_dada2/taxa_wsp.csv')

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# WRITE METADATA ----------------------------------------------------------

seqtab.nochim<-readRDS('4_dada2/seqtab.nochim.RDS')
samples.out <- rownames(seqtab.nochim)
samples.out <- gsub("\\.", "_", samples.out)
samples.out <- sub("\\-.*", "", samples.out)

day <- sapply(strsplit(samples.out, "_"), `[`, 1)
day <- as.numeric(gsub("day","",day))
well <- sapply(strsplit(samples.out, "_"), `[`, 2)
community <- sapply(strsplit(samples.out, "_"), `[`, 3)
rep <- as.numeric(sapply(strsplit(samples.out, "_"), `[`, 4))
seq.run<-c()
seq.run[grep('052214DR16s', samples.out)]<-'052214DR16s'
seq.run[grep('021216DR515F', samples.out)]<-'021216DR515F'
levels(as.factor(seq.run))

seq.run

samdf <- data.frame(community=community, day=day, rep=rep, seq.run=seq.run, well=well)
write.csv(samdf,'4_dada2/samdf.csv')


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
