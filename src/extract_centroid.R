##########################################
#  extract_centroid.R
##########################################
# This function extracts for a given ASV
# table a centroid of the table or of a subset
# of samples determined by a vector of
# factor levels present in a metadata file.
# The centroid is the median value of the
# relative abundance of each ASV in the subset.
# INPUT
#    ASV.table = A ASV table, with ASV in XX and samples
#         in YY
#    sample.md = A df with the metadata for each sample
#    factor.vec = A char vector with the factors that should be
#        considered to extract the subset of samples
#    level.vec = A char vector with the level that should be
#        considered for each factor in factor.vec to
#        subset samples. Levels should be given in the
#        same order than factors (i.e. both vectors are paired).
#    sample.id = A string with the name of the column where the
#        names of the samples are found in the metadata table
# OUTPUT
#    A vector with the relative abundance of the centroid.
# REQUIRES extract_subset.R
########################################
# author = apascualgarcia.github.io
# date = Dec 6th, 2022. ETH-Zürich
##########################################

extract_centroid = function(ASV.table, sample.md = data.frame(), 
                            factor.vec = c(), level.vec = c(),
                            sample.id = c()){
  # --- Extract the desired subset of samples
  ASV.sub = extract_subset(ASV.table, sample.md = sample.md, 
                 factor.vec = factor.vec, level.vec = level.vec,
                 sample.id = sample.id)
 # --- Compute the centroid
 # ..... normalize the table
  ASV.relAb = apply(ASV.sub, MARGIN = 1,FUN = function(x)(x/sum(x))) # Normalize the table
  ASV.relAb = t(ASV.relAb)
 # ..... resample the table with those probabilities
  set.seed(1323)
  
  ASV.sub.sam = matrix(0,nrow=nrow(ASV.relAb),ncol=ncol(ASV.relAb))

  for(i in 1:nrow(ASV.relAb)){
    sub.sample= table(sample(size=10000, x=seq(from=1, to=ncol(ASV.relAb)),
                             prob= ASV.relAb[i,], replace=T))
    ASV.sub.sam[i,as.numeric(names(sub.sample))] = sub.sample
  }
  centroid = apply(ASV.sub.sam, MARGIN = 2, FUN = median)
  centroid = centroid / sum(centroid)
  names(centroid) = colnames(ASV.relAb)
  return(centroid)
}