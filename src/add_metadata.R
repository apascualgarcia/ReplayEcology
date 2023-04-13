add_metadata=function(sample_md,partition,label){
  sample_md[[label]]=NA
  matched=match(rownames(partition),rownames(sample_md))
  sample_md[[label]][matched]=partition
  return(sample_md)
}