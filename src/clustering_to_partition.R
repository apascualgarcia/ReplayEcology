clustering_to_partition = function(data.dist,Kmax=2,
                                        outputfile='Partition_PamClustering_indexG1',
                                        outputlabel=c(),
                                        print.file = FALSE){
  # this function performs a pam clustering with data.dist,
  # joining elments until Kmax cluster are left data.dist and writes
  # a file with the identity of the cluster each element belongs to. If print.file
  # is TRUE, the partition is printed in a file. If provided,
  # outputlabel is used to complement the name of the file outputfile 
  # (the extension will always be .vec)
  # Fix here the optimal number of clusters
  data.cluster=pam.clustering(data.dist, Kmax)
  data.cluster.matrix=as.matrix(data.cluster)
  rownames(data.cluster.matrix)=attr(data.dist,"Labels")
  colnames(data.cluster.matrix)=c("PartId")
  if(print.file == TRUE){
    if(length(outputlabel)>0){outputlabel=paste0(outputlabel,"_")}
    fileOut=paste0(outputfile,"_",outputlabel,"Kmax-",Kmax,".vec")  
    write.table(data.cluster.matrix,
                file=fileOut,
                quote = FALSE,
                row.names = TRUE)
  }
  return(data.cluster.matrix)
}
