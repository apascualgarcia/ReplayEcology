plot_calinski=function(nclusters,filename="Plot_PAM-IndexG1medoids_JSD",
                       outputlabel=c(),
                      xlab="k clusters",ylab="Calinski-Harabasz index",
                      main="Jensen Shannon Divergence",...){
  # This function plots the output of otu_table_to_classes.R 
  # in a pdf file
  if(length(outputlabel)>0){filename=paste0(filename,"_")}
  fileOut=paste0(filename,outputlabel,".pdf")  
  pdf(file=fileOut, width = 10,heigh=10)
  par(cex.axis=3,cex.lab=3,cex.main=3,mar=c(11,11,8,8),mgp=c(6,2,0))
  plot(nclusters, type="h", col="red", lwd=20,lend=1,
       xlab=xlab, ylab=ylab,main=main,...)
  dev.off()
}


