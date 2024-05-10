###################################
#    clean_ASV_table.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-17
# Script Name: name_of_script.R  
# Script Description: This script takes the ASV table
# and metadata, cleans the ASV table removing low
# frequency ASVs and low-depth samples, and then
# gives two options two further select samples:
# --- match.exp = FALSE # if TRUE it will keep only
#                       # samples that were resurrected
# --- exclude.exp = c("4M") # a character vector indicating
#                     # the experiments to exclude, following
#                     # the levels defined in the samples' metadata
#                     # field "Experiment". Defaults to "4M" and it will
#                     # always be removed

clean_ASV_table = function(ASV.table, sample_md,
                           match_exp = FALSE,
                           exclude_exp = c("4M"),
                           nreads = 10000){
  # CLEAN DATA ---------- 
  
  exclude_exp_default = c("4M")
  
  # Clean data   ------
  ASV.table=ASV.table[,colSums(ASV.table) > 0]
  ASV.table=ASV.table[rowSums(ASV.table) > nreads, ]
  which(is.na(ASV.table))
  
  # .... exclude samples that didn't pass quality control (not present in samples metadata)
  matched=match(row.names(ASV.table),sample_md$sampleid) #attr(data.dist,"Labels"))
  ASV.table=ASV.table[!is.na(matched),]
  dim(ASV.table)
  
  # ... rebuild the metadata, it may have more samples
  matched=match(row.names(ASV.table),sample_md$sampleid)
  matched=matched[!is.na(matched)]
  sample_md = sample_md[matched, ]
  
  # Exclude data ----
  if((length(exclude_exp) != 0) | (match_exp == TRUE)){
    sample_md_red = sample_md
    
    # --- Always exclude 4M    
    idx = which(sample_md_red$Experiment != exclude_exp_default)
    sample_md_red = sample_md_red[idx, ]
    
    # --- Continue  matching resurrected communities
    if(match_exp == TRUE){ # exclude communities that were not resurrected
      nmatched = 5 # parent community id should appear 5 times
      start_in_exclude = grep("0D", exclude_exp)
      final_in_exclude = grep("7D", exclude_exp)
      nexcluded = length(start_in_exclude)-length(final_in_exclude)
      if(nexcluded == nmatched){
        stop("It doesn't make sense excluding 0D and 7D communities and match_exp = T") 
      }
      nchild = table(sample_md_red$parent)
      idx = which(nchild == nmatched)
      matched = match(sample_md_red$parent, names(nchild)[idx]) 
      sample_md_red = sample_md_red[!is.na(matched), ]
    }
    
    # --- Finally, remove specific experiments or replicates
    if(length(exclude_exp) > 1){ # exclude experiments beyond 4M
      #browser()
      for(level in exclude_exp){
        idx = which(sample_md_red$Experiment != level)
        sample_md_red = sample_md_red[idx, ]
      }
    }
    sample_md = sample_md_red #
    
    # ... look for the positions of remaining samples
    matched=match(sample_md$sampleid,row.names(ASV.table))
    matched=matched[!is.na(matched)]
    
    # ... reshape objects 
    ASV.tmp=ASV.table[matched,] # this tmp file is unnecessary, but easier to debug
    dim(ASV.tmp)
    ASV.table = ASV.tmp # this tmp file is unnecessary, but easier to debug
    ASV.table=ASV.table[,colSums(ASV.table) > 0]
    # ... rebuild the metadata, it may have more samples
    matched=match(row.names(ASV.table),sample_md$sampleid)
    matched=matched[!is.na(matched)]
    sample_md = sample_md[matched, ]
  }
  return(list("ASV.table" = ASV.table,"sample_md" = sample_md))
}


