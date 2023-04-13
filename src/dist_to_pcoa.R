dist_to_pcoa = function(data.dist,metadata,factor1,factor2=c(),
                       outputlabel=factor1,vec.coor=c(1,2,3)){
  # This function performs a pcoa of data.dist and represents it
  # considering as "factor1" for colors and, if given, "factor2"
  # for the center of mass. The function will create a pdf with
  # the plot in the directory in which the function is executed.
  # The PCoA coordinates to plot are determined by the vector vec.coor.
  # The data.dist names will be represented.
  #browser()
  coor.lab=paste0("Coor",paste0(vec.coor,collapse = ""))

  # ... Prepare a color palette
  library(RColorBrewer)
  n <- 60
  Kmax=length(unique(metadata[[factor1]]))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector = col_vector[1:(Kmax)]
  colorCodes = c("1"="pink2", "2"="green4", "3"="red1", "4"="gold","5"="deepskyblue1", "6"="violet") # Shannon
  col_vector=c(colorCodes,col_vector)
  # --- Principal coordinate analysis
  obs.pcoa=dudi.pco(as.dist(data.dist), scannf=F, nf=3)
  # --- Plot output
  id.els.exist=!is.na(metadata[[factor1]]) # there may be more elements in the metadata
  factor.plot=as.factor(metadata[[factor1]][id.els.exist])
  fileOut=paste0("Plot_PCoA_",coor.lab,"_",outputlabel,".pdf")
  pdf(file=fileOut, width = 10,heigh=10)
  s.class(obs.pcoa$li[,vec.coor], fac=factor.plot,
          cellipse=0,axesell = FALSE,cstar=1,
          grid=F,col=col_vector,clabel = 0)#, sub="Coordinates 1 and 2")
  text(obs.pcoa$li[,vec.coor],label=rownames(obs.pcoa$tab),
       adj=c(-.1,-.8),cex=0.5) # Add the names of the samples
  if(length(factor2) > 0){
    factor.plot2=as.factor(metadata[[factor2]][id.els.exist])
    s.class(obs.pcoa$li[,vec.coor], fac=factor.plot2, grid=F, 
          cstar=0,cpoint=0, cellipse=0,clabel=1,add.plot=TRUE) # Add the names of the
  }
  # Replot turning off all the properties of the graphic (dots, ellipses, 
  # etc) except the centers of mass given by factor2
  s.class(obs.pcoa$li[,vec.coor], fac=factor.plot,
            cellipse=0,axesell = FALSE,cstar=0,clabel=2,cpoint=0,
            grid=F,col=col_vector,add.plot=TRUE) 

  dev.off()
}
