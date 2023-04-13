################################
# extract_dist_vec_feat.R
################################
#
# This function takes a named numeric vector 
# and it returns a number of metrics, such as
# the maximum, minimum, mean, and identity (cell name)
# of these quantities. It returns a list with
# a numeric vector (for values) and a character
# vector (for identities). It additionally has
# a variable "label" to include a character term
# for relabeling the names of the vector for
# multiple calls.
# NOTE: NAs are removed
#
extract_dist_vec_feat = function(dist.vec,label = ""){
  #browser()
  maxDist = max(dist.vec,na.rm = T)
  id.max = names(which.max(dist.vec))
  minDist = min(dist.vec,na.rm = T)
  id.min = names(which.min(dist.vec))
  meanDist = mean(dist.vec,na.rm = T)
  diffDist = maxDist - minDist
  diffDistRel = diffDist / meanDist
  names_vec = c("maxDist","minDist","meanDist","diffDist","diffDistRel")
  names_vec = paste0(names_vec,label)
  dist.vec = c(maxDist, minDist, 
               meanDist, diffDist, diffDistRel)
  names(dist.vec) = names_vec
  names_vec = c("id.Dmax", "id.Dmin")
  names_vec = paste0(names_vec,label)
  char.vec = c(id.max, id.min)
  names(char.vec) = names_vec
  return(list("dist.vec" = dist.vec, "char.vec" = char.vec))
}