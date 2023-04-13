







# ALTERNATIVE TAXONOMY ----------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("4_dada2/DECIPHER_2.20.0/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
saveRDS(taxid,'4_dada2/taxid.RDS')
#taxid<-readRDS('4_dada2/taxid.RDS')
write.csv(taxid,'4_dada2/taxid.csv')

# BONUS: HANDOFF TO PHYLOSEQ ----------------------------------------------

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")

library(Biostrings); packageVersion("Biostrings")

library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())