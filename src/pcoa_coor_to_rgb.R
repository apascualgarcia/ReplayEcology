###################################
## script.R
###################################
## Purpose of script:
#
#
###################################
## Author: Dr. Alberto Pascual-García
## Date Created: 2023-06-16
## MIT Licensed.
## Contact: apascualgarcia.github.io
#################################### 
# This function takes a phyloseq object and an ordination and
# creates three color vectors with a gradient mapping the
# three coordinates of the ordination (one for each combination
# 1+2, 1+3 and 2+3). The match of colors to datapoints is specific
# of the treeholes data, but  it can be easily modified for
# other datets. In this case, we want the same colors for each
# dot according to the identity of the parent community (i.e., the
# four daughter replicates will have the same color than the parent)
## INPUT: pyloseq object and ordination
## OUTPUT:  A list with the three colors in the positions,
#         k = coordinate1 + coordinate2 (3, 4 and 5)
#####################################  

pcoa_coor_to_rgb = function(phyl.df, phyl.ord){

  # ... identify parents
  id.parents = which(phyl.df@sam_data$ExpCompact == "Starting")
  parents = as.character(phyl.df@sam_data$sampleid[id.parents])
  
  # ... look for their coordinates
  matched = match(parents, rownames(phyl.ord$vectors))
  parents.coor = phyl.ord$vectors[matched, ]
  
  # ... rescale coordinates to match [0,1] rgb interval
  for(i in 1:3){
    min.coor = min(parents.coor[, i])
    if(min.coor < 0){
      parents.coor[, i] = parents.coor[, i] + abs(min.coor)
    }else{
      parents.coor[, i] = parents.coor[, i] - min.coor
    }
    max.coor = max(parents.coor[, i])
    s = 1/max.coor
    parents.coor[, i] = s* parents.coor[, i]
  }
  color_list = list()
  parents.names = rownames(parents.coor)
  # We generate a different color vector for each of the three
  # combinations of PCoA axes (1,2) (1,3) and (2,3)
  for(u in 1:2){
    for(v in (u+1):3){
      k = u+v
      if(k == 3){ # three possibilities, we fix to 0 the free coordinate
        color_vec = rgb(parents.coor[,u],
                        parents.coor[,v],0)
      }else if((u+v) == 4){
        color_vec = rgb(parents.coor[,u],0,
                        parents.coor[,v])
      }else{
        color_vec = rgb(0,parents.coor[,u],
                        parents.coor[,v])
      }
      names(color_vec) = parents.names
      color_list[[k]] = color_vec
    }
  }
  
  # ... finally, extend the vector beyond parents
  
  matched = match(as.character(phyl.df@sam_data$parent),
                  parents.names)
  for(k in c(3,4,5)){
    color_list[[k]] = color_list[[k]][matched]
  }
  
  return(color_list)
}