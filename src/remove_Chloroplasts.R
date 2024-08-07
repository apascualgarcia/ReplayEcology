##################################################
# remove_Choroplasts.R
##################################################
# In this script I read the OTU table and the list of
# ASVs classified as Chloroplasts and I remove them.
# Zürich, August 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################

this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(dirSrc)

# --- Read ASVs table
setwd("../6_finalfiles")
fileOTU="seqtable_matchedandfiltered.csv"
ASV.table=read.table(fileOTU,sep="\t")
dim(ASV.table)
ASV.table=as.matrix(ASV.table)
colnames(ASV.table)[1:5]
rownames(ASV.table)[1:5]

# --- Read list  of Chloroplasts
fileChloro="ASVs_Chloroplasts.list"
ASV.chloro=read.table(fileChloro)
dim(ASV.chloro)

# --- Remove them from the table
matched=match(colnames(ASV.table),ASV.chloro$V1)
ASV.table.clean=ASV.table[,is.na(matched)]
dim(ASV.table.clean)

# --- Read list  of Mitochondria
fileMito="ASVs_Mitochondria.list"
ASV.mito=read.table(fileMito)
dim(ASV.mito)

# --- Remove them from the table
matched=match(colnames(ASV.table.clean),ASV.mito$V1)
ASV.table.clean=ASV.table.clean[,is.na(matched)]
dim(ASV.table.clean)

# --- Write final ASV table
fileOTU="seqtable_readyforanalysis.csv"
write.table(ASV.table.clean,file=fileOTU,quote=FALSE,sep="\t")

fileOTUt="seqtable_readyforanalysis.t.csv"
write.table(t(ASV.table.clean),file=fileOTUt,quote=FALSE,sep="\t")
