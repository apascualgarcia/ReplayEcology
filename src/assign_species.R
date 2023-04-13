here::here()

# READ IN TAXA TABLE ------------------------------------------------------

taxa<-readRDS('4_dada2/taxa.RDS')

# ASSIGN SPECIES-LEVEL TAXONOMY ------------------------------------

taxa_wsp <- addSpecies(taxa, "4_dada2/silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa_wsp, '4_dada2/taxa_wsp.RDS')
write.csv(taxa_wsp, '4_dada2/taxa_wsp.csv')

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
