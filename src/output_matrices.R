output_matrices = function(mat,name.mat="",pathOut=".",
                           plot.heatmap=FALSE,
                           par.heatmap = list(),par.pdf = list()){
  # This script takes a matrix and prints it to file, further representing
  # it with a heatmap if the option is selected.
  # mat = matrix to plot
  # name.mat = character, additional name for output files
  # pathOut = optional directory for output files
  # par.heatmap = parameters for heatmap2, it will override defaults.
  # par.pdf = parameters for pdf, it will override defaults.
  
  library(reshape2)
  
  current.dir=getwd()
  setwd(pathOut)
  fileOut=paste0("matrix_",name.mat,".tsv")
  write.table(mat,file=fileOut,sep="\t",quote = FALSE)
  #browser()
  if(plot.heatmap == TRUE){
    if(any(dim(mat) < 2)){
      mes=paste("** I can hardly generate a heatmap for a solution
              with a single value, please review dimensions of the matrix. 
                I leave the function...")
      warning(mes)
      return()
    }
    min.mat=min(mat)
    max.mat=max(mat)
    if((min.mat < 0) & (max.mat > 0)){
      # The commented lines cannot be used with do.call
      #my_palette <- colorRampPalette(c("blue","white","red"))(n = 99) # 100 percentcol_breaks = c(seq(-1,0,length=1)),  
      #col_breaks = c(seq(min.mat,-0.5,length=45),
      #               seq(-0.49,0.49,length=10),
      #               seq(0.5,max.mat,length=45))              
                     #seq(31,60,length=45))
      #my_color_fun=function(){my_palette}
      palette="hcl.colors"
      #additional.args=list(col = my_color_fun, breaks = col_breaks)
      additional.args=list(col = palette)
    }else{
      palette="heat.colors"
      additional.args=list(col = palette)
    }
    
    default.args=list(density.info="none",
                      dendrogram="both",
                      trace="none",
                      keysize = 1,
                      key.xlab = "Mean JSD",
                      key.title = "",
                      cexRow = 1.25,cexCol = 1.25, cex.lab=3,
                      key.par = list(cex.main=2,cex.axis=1.2,
                                     cex.lab=1.5), #usr=c(0,1,0,1)),
                      margins=c(12,12),
                      col=palette)                  
    final.args=functionArgsList(default.args,additional.args)
    final.args=functionArgsList(final.args,par.heatmap)
    
    # ... finally plot   
    file.heatmap=paste("heatmap-mat_",name.mat,".pdf",sep="")
    default.pdf=list(file=file.heatmap,width=12,height=12)
    final.pdf=functionArgsList(default.pdf,par.pdf)
    do.call("pdf",final.pdf)
    do.call("heatmap.2",c(list(mat),final.args))
    dev.off()
  }
  setwd(current.dir)
}