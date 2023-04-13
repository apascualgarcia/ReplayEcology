print_pamk_to_partition = function(out_classes,
                                   outputfile='Partition_PamK_indexCH',
                                   outputlabel=c(),
                                   print.file = FALSE){
  # this function performs extracts the partition obtained with pamk 
  # of package fpc and prints it to a file
  data.cluster=as.data.frame(out_classes$pamobject$clustering)
  data.cluster.matrix=as.matrix(data.cluster)
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
